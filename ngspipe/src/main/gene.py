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

def merge_gene_entries(dataframe, outfile): # return a DataFrame
    #  refgene query, save known to gene_df and unknown to unknown_df
    headers = config.merged_custom_header
    if( not config.refgene_header[0] in headers):
        print('Error: refGene annotation is not in')
        exit()
    print('\n**************************************************************')
    print('*    merging gene entries from refgene, knowngene, & ensgene *')
    print('**************************************************************')
    ofstream = open(outfile, 'w')
    k=0
    #write the header to outputfile
    
    for item in headers:
        k=k+1
        if (k==len(headers)):
            ofstream.write(str(item))
        else:
            ofstream.write(str(item)+'\t')
    ofstream.write('\n')
    
    
    j=0
    for index, row in dataframe.iterrows():
        
        temp_row=row.copy(deep=True)
        temp_row = temp_row[headers]
        #temp_row['Gene'] ='a'
        #print(row[headers])
        #print(temp_row[headers])
        write_row(temp_row, ofstream)
        #print(row['VariantFunction'], row['VariantFunction.k'])
        #print(row['ExonicFunction'], row['ExonicFunction.k'])
        if(row['VariantFunction'] != row['VariantFunction.k'] \
           or row['ExonicFunction'] != row['ExonicFunction.k']) :
            #print('refGene & knownGene info is not the same')
            #print(temp_row['Chr'], temp_row['Start'], temp_row['End'], \
            #      temp_row['Ref'], temp_row['Alt'])
            temp_row['VariantFunction'] = row['VariantFunction.k']
            temp_row['Gene'] = row['Gene.k']
            temp_row['GeneDetail'] = row['GeneDetail.k']
            temp_row['ExonicFunction'] = row['ExonicFunction.k']
            temp_row['AAChange'] = row['AAChange.k']
            write_row(temp_row, ofstream)
        if(row['VariantFunction'] != row['VariantFunction.e'] \
           or row['ExonicFunction'] != row['ExonicFunction.e']):
            #print('refGene & ensGene info is not the same')
            #print(temp_row['Chr'], temp_row['Start'], temp_row['End'], \
            #      temp_row['Ref'], temp_row['Alt'])
            temp_row['VariantFunction'] = row['VariantFunction.e']
            temp_row['Gene'] = row['Gene.e']
            temp_row['GeneDetail'] = row['GeneDetail.e']
            temp_row['ExonicFunction'] = row['ExonicFunction.e']
            temp_row['AAChange'] = row['AAChange.e']
            write_row(temp_row, ofstream)
            
    ofstream.close()
    

#####################################################################
# expand_gene_entries
#####################################################################
def expand_gene_entries(dataframe, outfile):
    print('\n**************************************************************')
    print('*    Expanding gene entries                                  *')
    print('**************************************************************')
    ofstream = open(outfile, 'w')
    k=0
    #write the header to outputfile
    headers = list(dataframe.columns.values)
    for item in headers:
        k=k+1
        if k==len(headers):
            ofstream.write(str(item))
        else:
            ofstream.write(str(item)+'\t')
    ofstream.write('\n')
    
    j=0
    for index, row in dataframe.iterrows():
        #iterate the dataframe to find the row with multiple gene entries
        gene_piece = row['Gene'].strip().split(',')
        aa_change_piece = row['AAChange'].strip().split(',')
        
       
        #print('*************************')
        #print(row['Gene']) 
        #print('*************************')  
        
        for gene_item in gene_piece:
            for aa_item in aa_change_piece:
                #print('gene: ', gene_piece)
                row['Gene'] = gene_item
                row['AAChange']= aa_item
                write_row(row, ofstream)   
         
    ofstream.close()

def write_row(row, outstream):
    k=0
    for item in row:
        #print(type(item))
        i = str(item)
        outstream.write(i.strip())
        k=k+1
        if(k!=len(row)):
            outstream.write('\t')
    outstream.write('\n')
            
    
