#! /usr/local/bin/python3.4
#####################################################################
#
#    anv.py:  the file contains all the gene level functions related to 
#             annovar
#    author:  Liyong Deng
#    copyright (c) 2014 by Liyong Deng
#
#####################################################################
import os
import subprocess, shlex

import pandas as pd

import config



def fetch_refgene(inputfile, rflag): # return a DataFrame
    print('\n**************************************************************')
    print('*    Start to fetch geneinfo from refgene                    *')
    print('**************************************************************')
    cmd_str = 'annotate_variation.pl --geneanno -out ' + inputfile   \
            + ' -build ' + config.buildver + ' ' + inputfile + ' ' \
            + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    #  read in exonic file
    result = inputfile+config.gene_exn_suffix
    exn_df = pd.read_csv(result, header=None, names=config.refgene_exn_header,\
                         sep='\t', usecols=config.refgene_exn_use_header,\
                         low_memory=False)
    #  If a dataframe is empty, it contains all columns in names in pd.read_csv
    if (len(exn_df.index) == 0):
        exn_df = pd.DataFrame(columns=config.refgene_exn_use_header)
    if(rflag):
        os.remove(result)
    #  read in variant file
    result = inputfile+config.gene_var_suffix    
    var_df = pd.read_csv(result, header=None, names=config.refgene_var_header,\
                         sep='\t', usecols=config.refgene_var_use_header, \
                         low_memory=False)
    #  If a dataframe is empty, it contains all columns in names in pd.read_csv
    if (len(var_df.index) == 0):
        var_df = pd.DataFrame(columns=config.refgene_var_use_header)
    if(rflag):
        os.remove(result)
    #  merge the information from exonic and variant files    
    var_df = var_df[['exonic' not in x for x in var_df.variant_function]]    
    gene_df = exn_df.append(var_df)   
    gene_df[config.basic_header] = gene_df[config.basic_header].astype(str)
    return gene_df

def fetch_knowngene(inputfile, rflag): # return a DataFrame
    print('\n**************************************************************')
    print('*    Start to fetch geneinfo from knowngene                  *')
    print('**************************************************************')
    cmd_str = 'annotate_variation.pl --geneanno -out ' + inputfile   \
            + ' -build ' + config.buildver + ' ' + inputfile + ' ' \
            + config.dir_humandb + '  -dbtype ' + config.knowngene
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    #  read in exonic file
    result = inputfile+config.gene_exn_suffix
    exn_df = pd.read_csv(result, header=None, \
                         names=config.knowngene_exn_header,\
                         sep='\t', usecols=config.knowngene_exn_use_header,\
                         low_memory=False)
    #  If a dataframe is empty, it contains all columns in names in pd.read_csv
    if (len(exn_df.index) == 0):
        exn_df = pd.DataFrame(columns=config.knowngene_exn_use_header)
    if(rflag):
        os.remove(result)
    #  read in variant file
    result = inputfile+config.gene_var_suffix    
    var_df = pd.read_csv(result, header=None, \
                         names=config.knowngene_var_header, \
                         sep='\t', usecols=config.knowngene_var_use_header, \
                         low_memory=False)
    #  If a dataframe is empty, it contains all columns in names in pd.read_csv
    if (len(var_df.index) == 0):
        var_df = pd.DataFrame(columns=config.knowngene_var_use_header)
    if(rflag):
        os.remove(result)
        
    #  merge the information from exonic and variant files
    var_df = var_df[['exonic' not in x for x in var_df.k_variant_function]]   
    gene_df = exn_df.append(var_df)
    gene_df[config.basic_header] = gene_df[config.basic_header].astype(str)
    return gene_df

def fetch_ensgene(inputfile, rflag): # return a DataFrame
    print('\n**************************************************************')
    print('*    Start to fetch geneinfo from ensgene                    *')
    print('**************************************************************')
    cmd_str = 'annotate_variation.pl --geneanno -out ' + inputfile   \
            + ' -build ' + config.buildver + ' ' + inputfile + ' ' \
            + config.dir_humandb + ' -dbtype ' + config.ensgene
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    #  read in exonic file
    result = inputfile+config.gene_exn_suffix
    exn_df = pd.read_csv(result, header=None, \
                         names=config.ensgene_exn_header,\
                         sep='\t', usecols=config.ensgene_exn_use_header,\
                         low_memory=False)
    #  If a dataframe is empty, it contains all columns in names in pd.read_csv
    if (len(exn_df.index) == 0):
        exn_df = pd.DataFrame(columns=config.ensgene_exn_use_header)
    if(rflag):
        os.remove(result)
    #  read in variant file
    result = inputfile+config.gene_var_suffix    
    var_df = pd.read_csv(result, header=None, \
                         names=config.ensgene_var_header, \
                         sep='\t', usecols=config.ensgene_var_use_header, \
                         low_memory=False)
    #  If a dataframe is empty, it contains all columns in names in pd.read_csv
    if (len(var_df.index) == 0):
        var_df = pd.DataFrame(columns=config.ensgene_var_use_header)
    if(rflag):
        os.remove(result)
        
    #  merge the information from exonic and variant files
    var_df = var_df[['exonic' not in x for x in var_df.e_variant_function]]   
    gene_df = exn_df.append(var_df)
    gene_df[config.basic_header] = gene_df[config.basic_header].astype(str)
    return gene_df

def fetch_mergedgene(inputfile, rflag): # return a DataFrame
    #  refgene query, save known to gene_df and unknown to unknown_df
    refgene_df = fetch_refgene(inputfile, rflag)   
    refgene_unknown_df = refgene_df[refgene_df.variant_function=='unknown']
    gene_df = refgene_df[refgene_df.variant_function!='unknown']    
    #print(gene_df)
    unknown_df = refgene_unknown_df[config.basic_header]
    print('\nrefGene Number of unknown exonic entries:')
    print(len(unknown_df))
    # knowngene query, merge unknown_df with knowngene_df 
    knowngene_df = fetch_knowngene(inputfile, rflag)
    unknown_df = pd.merge(unknown_df, knowngene_df, how='left', \
                          on=config.basic_header)    
    unknown_df= unknown_df.rename(columns= \
                                  {config.knowngene_var_header[0]: \
                                   config.refgene_var_header[0], \
                                   config.knowngene_var_header[1]: \
                                   config.refgene_var_header[1]})    
    cols = unknown_df.columns.tolist()
    cols = cols[-2:] + cols[:-2]     
    unknown_df = unknown_df[cols]
    
    #  append unknows_df with known to gene_df
    gene_df = gene_df.append(unknown_df[unknown_df.variant_function \
                                        != 'unknown'])
    unknown_df=unknown_df[unknown_df.variant_function=='unknown']
    unknown_df = unknown_df[config.basic_header]
    print('\nknownGene Number of unknown exonic entries:')
    print(len(unknown_df))
    #  ensgene query to find the unknown after refgene and knowngene
    ensgene_df = fetch_ensgene(inputfile, rflag)
    unknown_df = pd.merge(unknown_df, ensgene_df, how='left', \
                          on=config.basic_header)    
    unknown_df= unknown_df.rename(columns= \
                                  {config.ensgene_var_header[0]: \
                                   config.refgene_var_header[0], \
                                   config.ensgene_var_header[1]: \
                                   config.refgene_var_header[1]})    
    cols = unknown_df.columns.tolist()
    cols = cols[-2:] + cols[:-2]
    unknown_df = unknown_df[cols]
    #print(unknown_df)
   
    gene_df = gene_df.append(unknown_df)
    #print(gene_df)
    print(len(gene_df.index))
    #gene_df = gene_df.append(unknown_df[unknown_df.variant_function \
       
    print('\nensGene Number of unknown exonic entries: ')
    print(len(unknown_df[unknown_df.variant_function == 'unknown']))
    gene_df[config.basic_header] = gene_df[config.basic_header].astype(str)
    
    return gene_df

#####################################################################
# expand_gene_entries
#####################################################################
def expand_gene_entries(dataframe, outfile):
    print('\n**************************************************************')
    print('*    Expanding gene entries                                  *')
    print('**************************************************************')
    ofile = open(outfile, 'w')
    k=0
    #write the header to outputfile
    for item in list(dataframe.columns.values):
        k=k+1
        if k==len(list(dataframe.columns.values)):
            ofile.write(str(item))
        else:
            ofile.write(str(item)+'\t')
    ofile.write('\n')
    
    j=0
    for index, row in dataframe.iterrows():
        #iterate the dataframe to find the row with multiple gene entries
        gene_piece = row['geneinfo'].strip().split(',')
        gene_piece= gene_piece[:-1]
        #print(gene_piece)
        if(len(gene_piece) <= 1):
            #  one gene entry: only copy it into outputfile
            k=0
            for item in row:
                k = k+1;
                if k==len(row):
                    ofile.write(str(item))
                else:
                    ofile.write(str(item)+'\t')
            
            ofile.write('\n')
            j=j+1
        
        else:   
            #multiple gene entries 
            for i in range(len(gene_piece)):            #print(i)
                row['geneinfo'] = gene_piece[i]            
                #temp_df.loc[j] = row
                k=0
                for item in row:
                    k = k+1;
                    if k==len(row):
                        ofile.write(str(item))
                    else:
                        ofile.write(str(item)+'\t')
                ofile.write('\n')
            j=j+1
     
    ofile.close()

