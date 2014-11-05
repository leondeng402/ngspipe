#! /usr/local/bin/python3.4
#####################################################################
#
#    region.py:  the file contains all the region level functions  
#                 related to annovar
#    author:  Liyong Deng
#    copyright (c) 2014 by Liyong Deng
#
#####################################################################
import os
import subprocess, shlex

import pandas as pd

import config

def fetch_cytoband(inputfile, rflag): # return a DataFrame
    print('\n**************************************************************')
    print('*    Start to fetch cytoband info                            *')
    print('**************************************************************')
    md_str= 'annotate_variation.pl -regionanno  -out ' + inputfile \
            + ' -dbtype ' + config.cytoband +' -build ' + config.buildver \
            + ' ' + inputfile + ' ' + config.dir_humandb
    cmd = shlex.split(md_str)
    subprocess.call(cmd)
    result = inputfile + '.' + config.buildver + config.cytoband_suffix
    
    cyto_df =  pd.read_csv(result, header=None, names=config.cytoband_header, \
                          sep='\t', usecols=config.cytoband_use_header, \
                          low_memory=False)
    if(rflag):
        os.remove(inputfile + '.' + config.buildver + config.cytoband_suffix)
    #  If a dataframe is empty, it contains all columns in names in pd.read_csv
    if (len(cyto_df.index) == 0):
        cyto_df = pd.DataFrame(columns=config.cytoband_use_header)
    cyto_df[config.basic_header] = cyto_df[config.basic_header].astype(str)
    return cyto_df
        
def fetch_superdup(inputfile, rflag): # return a DataFrame
    print('\n**************************************************************')
    print('*    Start to fetch segmental duplications info              *')
    print('**************************************************************')
    md_str= 'annotate_variation.pl -regionanno  -out ' + inputfile \
            + ' -dbtype ' + config.superdup +' -build ' + config.buildver \
            + ' ' + inputfile + ' ' + config.dir_humandb
    cmd = shlex.split(md_str)
    subprocess.call(cmd)
    result = inputfile + '.' + config.buildver + config.superdup_suffix
    
    superdup_df =  pd.read_csv(result, header=None, \
                               names=config.superdup_header, \
                               sep='\t', usecols=config.superdup_use_header, \
                               low_memory=False)
    if(rflag):
        os.remove(inputfile + '.' + config.buildver + config.superdup_suffix)
    #  If a dataframe is empty, it contains all columns in names in pd.read_csv
    if (len(superdup_df.index) == 0):
        superdup_df = pd.DataFrame(columns=config.superdup_use_header)
    superdup_df[config.basic_header] = superdup_df[config.basic_header].astype(str)
    return superdup_df

def fetch_gwas(inputfile, rflag): # return a DataFrame
    print('\n**************************************************************')
    print('*    Start to fetch variants info in published gwas          *')
    print('**************************************************************')
    md_str= 'annotate_variation.pl -regionanno  -out ' + inputfile \
            + ' -dbtype ' + config.gwas +' -build ' + config.buildver \
            + ' ' + inputfile + ' ' + config.dir_humandb
    cmd = shlex.split(md_str)
    subprocess.call(cmd)
    result = inputfile + '.' + config.buildver + config.gwas_suffix
    
    gwas_df =  pd.read_csv(result, header=None, names=config.gwas_header, \
                          sep='\t', usecols=config.gwas_use_header, \
                          low_memory=False)
    if(rflag):
        os.remove(inputfile + '.' + config.buildver + config.gwas_suffix)
    #  If a dataframe is empty, it contains all columns in names in pd.read_csv
    if (len(gwas_df.index) == 0):
        gwas_df = pd.DataFrame(columns=config.gwas_use_header)
    gwas_df[config.basic_header] = gwas_df[config.basic_header].astype(str)
    return gwas_df