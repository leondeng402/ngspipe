#! /usr/local/bin/python3.4
#####################################################################
#
#    ops.py:  the file contains post-annotation operations, 
#             including filtering, and ranking
#    author:  Liyong Deng
#    copyright (c) 2014 by Liyong Deng
#
#####################################################################
import os
import subprocess, shlex

import pandas as pd

import config

def apply_filter(dataframe, input): # return DataFrame
    # parse parameters maf, varantfunction, exonicfunction, 
    # nonsynonymous SNV pathogenicity parameter list.
    print('parsing filter file')
    maf = 0
    vflist=[]
    eflist=[]
    snv_pathogenic_parameters ={}
    
    with open(input, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            #print(line)
            if(line == '' or line.startswith('#')):
                continue
            formula = line.split('&')
        
            if(len(formula) == 1):
                words=formula[0].split('=')
                if(words[0] == 'maf'):
                    maf=float(words[1])
                elif(words[0] == 'VariantFunction'):
                    vflist=words[1].split(',')
                elif(words[0] == 'ExonicFunction'):
                    eflist=words[1].split(',')
                else:
                    print(words[0])
                    print('Error: not correct grammar for a filter file')
            elif(len(formula)== 2):
                words0 = formula[0].split('=')
                words1 = formula[1].split('=')
                words0[0]=words0[0].strip()
                words0[1]=words0[1].strip()        
                if(words0[0]=='ExonicFunction' and words0[1] == \
                   'nonsynonymous SNV'):
                    words1[0]=words1[0].strip()
                    words1[1]=words1[1].strip()
                    snv_pathogenic_parameters[words1[0]]=words1[1]
                
                else: 
                    print(formula, words0[0], words0[1], '2')
                    print('Error: not correct grammar for a filter file')
            else:
                print(formula, '>2')
                print('Error: not correct grammar for a filter file')
            
    # start to filter after gathering parameters
    print('\nmaf filtering')
    header = dataframe.columns.values
    nadf = dataframe[dataframe['PopFreqMax'] == '.']
    dataframe = dataframe[dataframe['PopFreqMax'] != '.']
    dataframe.is_copy = False
    dataframe[['PopFreqMax']] = dataframe[['PopFreqMax']].astype(float)
    dataframe = dataframe[dataframe['PopFreqMax'] <= maf]
    dataframe = dataframe.append(nadf)
    if(len(dataframe.index)==0):
        print('None DataFrame')
        return dataframe
    dataframe.to_csv('/home/leon/lab/ngspipe/maf_filtered.txt', index=False, sep='\t')
    print('Rows after maf filtering:', len(dataframe.index))
    print('\nVariantFunction filtering')
    if(len(vflist) > 0):
        df = pd.DataFrame(columns=header)
        for item in vflist:
            item = item.strip()
            #print(item)
            tempdf = dataframe[dataframe['VariantFunction']==item]
            df=df.append(tempdf)
        dataframe=df
    if(len(dataframe.index)==0):
        print('None DataFrame')
        return dataframe
    dataframe.to_csv('/home/leon/lab/ngspipe/var_filtered.txt', index=False, sep='\t')
    print('Rows after var filtering:', len(dataframe.index))
    
    print('\nExonicFunction filtering')
    #print(dataframe.columns.values)
    if(len(eflist) > 0):
        df = pd.DataFrame(columns=header)
        for item in eflist:
            item = item.strip()
            tempdf = dataframe[dataframe['ExonicFunction']==item]
            df=df.append(tempdf)
        dataframe=df
        
    if(len(dataframe.index)==0):
        print('None DataFrame')
        return dataframe 
    dataframe.to_csv('/home/leon/lab/ngspipe/exn_filtered.txt', index=False, \
                     sep='\t')

    print('Rows after exn filtering:', len(dataframe.index))
    
    print('\nnonsynonymous pathogenicity filtering')
    restdf=dataframe[dataframe['ExonicFunction']!='nonsynonymous SNV']
    ns_snv_df=dataframe[dataframe['ExonicFunction']=='nonsynonymous SNV']
    print('rest: ', len(restdf.index), '\tns_snv: ', len(ns_snv_df.index))
    
    df = pd.DataFrame(columns=header)
    if(len(snv_pathogenic_parameters ) > 0):
        # iterate through pathogenicity parameters to do the filtering              
        for key, value in snv_pathogenic_parameters.items():           
            key=key.strip()
            value=value.strip()
            #print(key, value)
            if('CADD' in key): 
                # CADD score filtering                            
                nnadf= ns_snv_df[ns_snv_df[key] != '.']                 
                nadf = ns_snv_df[ns_snv_df[key] == '.']  
                if(len(nnadf.index)>0):
                    nnadf.is_copy = False    
                    nnadf[[key]]=nnadf[[key]].astype(float)   
                    matchdf = nnadf[nnadf[key] >= float(value)]
                    nmatchdf = nnadf[nnadf[key] < float(value)]
                    df = df.append(matchdf, ignore_index=True) 
                    ns_snv_df = nadf.append(nmatchdf, ignore_index=True)                
            else:
                # sift and polyphen2                                 
                matchdf = ns_snv_df[ns_snv_df[key]==value]                
                nmatchdf = ns_snv_df[ns_snv_df[key]!=value]                
                df = df.append(matchdf, ignore_index=True)
                ns_snv_df = pd.DataFrame(nmatchdf, index=None)           
    print('done with filtering')      
    dataframe = df.append(restdf, ignore_index=True)
    print('Rows after all filtering:', len(dataframe.index))
    return dataframe


def rank(dataframe, input):
    print('Parsing ranking file')
    with open(input, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            #print(line)
            if(line == '' or line.startswith('#')):
                continue
def retrieve_genotype(input): # return DataFrame
    if(not input.endswith('.vcf')):
        return None
    with open(input, 'r')as f:
        results=[]
        for line in f.readlines():
            line = line.strip()            
            if(line.startswith('#CHROM')):
                line = line.replace('#', '')
                words=line.split('\t')
                
                header=words[:5]+words[9:]
            elif(not line.startswith('##')):
                words=line.split('\t')
                variant=words[:5]
                result=words[9:]
                row=[]
                for words in result:                    
                    word=words.split(':')
                    row.append(word[0])
                results.append(variant+row)
    #print(ids)
    #print(results)
    genotype_df = pd.DataFrame(results, columns=header)
    genotype_df[header] = genotype_df[header].astype(str)
    genotype_df.rename(columns={header[0]:'Chr', header[1]:'Start', \
                                header[3]:'Ref', header[4]:'Alt' }, \
                       inplace=True)
    return genotype_df
                
        
def output_results(output, variant_df, genotype_df, mode, bedf): 
    if(mode == 'mono'):
        print('output the results in a single file')
        df = pd.merge(variant_df, genotype_df, how='left', \
                      on=['Chr', 'Start', 'Ref', 'Alt'])
        df.to_csv(output, sep='\t', index=False)
    elif(mode == 'proband'):
        print('output the results of proband only ')
    elif(mode == 'trios'):
        print('output the results of trios and singleton')
    elif(mode == 'family'):
        print('output the results of family and singleton')
    elif(mode == 'all'):
        print('output the results in a single file, proband, trios and family')
    else:
        print('Error: ' + mode +' wrong mode')
    
    