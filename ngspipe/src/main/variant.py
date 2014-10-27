#! /usr/local/bin/python3.4
#####################################################################
#
#    variant.py:  the file contains all the variant level functions  
#                 related to annovar
#    author:  Liyong Deng
#    copyright (c) 2014 by Liyong Deng
#
#####################################################################
import os
import subprocess, shlex

import pandas as pd

import config

def fetch_snpid(inputfile, rflag): # return a DataFrame
    print('\n**************************************************************')
    print('*    Start to fetch dbsnp ID                                 *')
    print('**************************************************************')
    md_str= 'annotate_variation.pl -filter  -out ' + inputfile + ' -dbtype ' \
        + config.snp +' -buildver ' + config.buildver + ' ' + inputfile \
        + ' ' + config.dir_humandb
    cmd = shlex.split(md_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.snp_drp_suffix
    snp_df =  pd.read_csv(result, header=None, names=config.snp_header, \
                          sep='\t', usecols=config.snp_use_header, \
                          low_memory=False)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.snp_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.snp_flt_suffix)
    #  If a dataframe is empty, it contains all columns in names in pd.read_csv
    if (len(snp_df.index) == 0):
        snp_df = pd.DataFrame(columns=config.snp_use_header)
    
    return snp_df

def fetch_hgmd(): # return DataFrame
    print('\n**************************************************************')
    print('*    Start to fetch hgmd                                     *')
    print('**************************************************************')
    hgmd_df = pd.read_csv(config.hgmd_file, header= None, 
                          names=config.hgmd_header, \
                          sep='!', index_col=False, 
                          usecols=config.hgmd_use_header, \
                          low_memory=False)

    hgmd_df.fillna('.')
    return hgmd_df
    
def fetch_g1k_all(inputfile, rflag): # return DataFrame    
    ### fetch g1k_all_allele_freq
    print('Fetch g1k_all allele frequency')
    cmd_str= 'annotate_variation.pl  -filter -out ' + inputfile   + ' -dbtype ' \
            + config.g1k + '_all -build ' + config.buildver + ' ' + inputfile \
            + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.g1k_all_drp_suffix
    g1k_all_df = pd.read_csv(result, header=None, names=config.g1k_all_header, \
                             sep='\t', usecols=config.g1k_all_use_header,\
                             low_memory=False)
    #g1k_all_df= g1k_all_df.drop('label', axis=1)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.g1k_all_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.g1k_all_flt_suffix)
    if (len(g1k_all_df.index) == 0):
        g1k_all_df = pd.DataFrame(columns=config.g1k_all_use_header)
    g1k_all_df[config.basic_header] = g1k_all_df[config.basic_header].astype(str)
    return g1k_all_df

def fetch_g1k_afr(inputfile, rflag): # return DataFrame    
    ### fetch g1k_afr_allele_freq
    print('Fetch g1k_afr allele frequency')
    cmd_str= 'annotate_variation.pl  -filter -out ' + inputfile   + ' -dbtype ' \
            + config.g1k + '_afr -build ' + config.buildver + ' ' + inputfile \
            + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.g1k_afr_drp_suffix
    g1k_afr_df = pd.read_csv(result, header=None, names=config.g1k_afr_header, \
                             sep='\t', usecols=config.g1k_afr_use_header,\
                             low_memory=False)
    #g1k_all_df= g1k_all_df.drop('label', axis=1)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.g1k_afr_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.g1k_afr_flt_suffix)
    if (len(g1k_afr_df.index) == 0):
        g1k_afr_df = pd.DataFrame(columns=config.g1k_afr_use_header)
    g1k_afr_df[config.basic_header] = g1k_afr_df[config.basic_header].astype(str)
    return g1k_afr_df

def fetch_g1k_amr(inputfile, rflag): # return Datamrame    
    ### fetch g1k_amr_allele_freq
    print('Fetch g1k_amr allele frequency')
    cmd_str= 'annotate_variation.pl  -filter -out ' + inputfile   + ' -dbtype ' \
            + config.g1k + '_amr -build ' + config.buildver + ' ' + inputfile \
            + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.g1k_amr_drp_suffix
    g1k_amr_df = pd.read_csv(result, header=None, names=config.g1k_amr_header, \
                             sep='\t', usecols=config.g1k_amr_use_header,\
                             low_memory=False)
    #g1k_all_df= g1k_all_df.drop('label', axis=1)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.g1k_amr_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.g1k_amr_flt_suffix)
    if (len(g1k_amr_df.index) == 0):
        g1k_amr_df = pd.DataFrame(columns=config.g1k_amr_use_header)
    g1k_amr_df[config.basic_header] = g1k_amr_df[config.basic_header].astype(str)
    return g1k_amr_df

def fetch_g1k_eas(inputfile, rflag): # return Dateasame    
    ### fetch g1k_eas_allele_freq
    print('Fetch g1k_eas allele frequency')
    cmd_str= 'annotate_variation.pl  -filter -out ' + inputfile   + ' -dbtype ' \
            + config.g1k + '_eas -build ' + config.buildver + ' ' + inputfile \
            + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.g1k_eas_drp_suffix
    g1k_eas_df = pd.read_csv(result, header=None, names=config.g1k_eas_header, \
                             sep='\t', usecols=config.g1k_eas_use_header,\
                             low_memory=False)
    #g1k_all_df= g1k_all_df.drop('label', axis=1)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.g1k_eas_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.g1k_eas_flt_suffix)
    if (len(g1k_eas_df.index) == 0):
        g1k_eas_df = pd.DataFrame(columns=config.g1k_eas_use_header)
    g1k_eas_df[config.basic_header] = g1k_eas_df[config.basic_header].astype(str)
    return g1k_eas_df

def fetch_g1k_eur(inputfile, rflag): # return Dateurame    
    ### fetch g1k_eur_allele_freq
    print('Fetch g1k_eur allele frequency')
    cmd_str= 'annotate_variation.pl  -filter -out ' + inputfile   + ' -dbtype ' \
            + config.g1k + '_eur -build ' + config.buildver + ' ' + inputfile \
            + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.g1k_eur_drp_suffix
    g1k_eur_df = pd.read_csv(result, header=None, names=config.g1k_eur_header, \
                             sep='\t', usecols=config.g1k_eur_use_header,\
                             low_memory=False)
    #g1k_all_df= g1k_all_df.drop('label', axis=1)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.g1k_eur_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.g1k_eur_flt_suffix)
    if (len(g1k_eur_df.index) == 0):
        g1k_eur_df = pd.DataFrame(columns=config.g1k_eur_use_header)
    g1k_eur_df[config.basic_header] = g1k_eur_df[config.basic_header].astype(str)
    return g1k_eur_df

def fetch_g1k_sas(inputfile, rflag): # return Datsasame    
    ### fetch g1k_sas_allele_freq
    print('Fetch g1k_sas allele frequency')
    cmd_str= 'annotate_variation.pl  -filter -out ' + inputfile   + ' -dbtype ' \
            + config.g1k + '_sas -build ' + config.buildver + ' ' + inputfile \
            + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.g1k_sas_drp_suffix
    g1k_sas_df = pd.read_csv(result, header=None, names=config.g1k_sas_header, \
                             sep='\t', usecols=config.g1k_sas_use_header,\
                             low_memory=False)
    #g1k_all_df= g1k_all_df.drop('label', axis=1)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.g1k_sas_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.g1k_sas_flt_suffix)
    if (len(g1k_sas_df.index) == 0):
        g1k_sas_df = pd.DataFrame(columns=config.g1k_sas_use_header)
    g1k_sas_df[config.basic_header] = g1k_sas_df[config.basic_header].astype(str)
    return g1k_sas_df

def fetch_esp_all(inputfile, rflag): # return Datsasame    
    ### fetch esp_all_allele_freq
    print('Fetch esp_all allele frequency')
    cmd_str= 'annotate_variation.pl  -filter -out ' + inputfile   + ' -dbtype ' \
            + config.esp + '_all -build ' + config.buildver + ' ' + inputfile \
            + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.esp_all_drp_suffix
    esp_all_df = pd.read_csv(result, header=None, names=config.esp_all_header, \
                             sep='\t', usecols=config.esp_all_use_header,\
                             low_memory=False)
    #esp_all_df= esp_all_df.drop('label', axis=1)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.esp_all_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.esp_all_flt_suffix)
    if (len(esp_all_df.index) == 0):
        esp_all_df = pd.DataFrame(columns=config.esp_all_use_header)
    esp_all_df[config.basic_header] = esp_all_df[config.basic_header].astype(str)
    return esp_all_df

def fetch_esp_aa(inputfile, rflag): # return Datsasame    
    ### fetch esp_aa_allele_freq
    print('Fetch esp_aa allele frequency')
    cmd_str= 'annotate_variation.pl  -filter -out ' + inputfile   + ' -dbtype ' \
            + config.esp + '_aa -build ' + config.buildver + ' ' + inputfile \
            + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.esp_aa_drp_suffix
    esp_aa_df = pd.read_csv(result, header=None, names=config.esp_aa_header, \
                             sep='\t', usecols=config.esp_aa_use_header,\
                             low_memory=False)
    #esp_aa_df= esp_aa_df.drop('label', axis=1)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.esp_aa_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.esp_aa_flt_suffix)
    if (len(esp_aa_df.index) == 0):
        esp_aa_df = pd.DataFrame(columns=config.esp_aa_use_header)
    esp_aa_df[config.basic_header] = esp_aa_df[config.basic_header].astype(str)
    return esp_aa_df

def fetch_esp_ea(inputfile, rflag): # return Datsasame    
    ### fetch esp_ea_allele_freq
    print('Fetch esp_ea allele frequency')
    cmd_str= 'annotate_variation.pl  -filter -out ' + inputfile   + ' -dbtype ' \
            + config.esp + '_ea -build ' + config.buildver + ' ' + inputfile \
            + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.esp_ea_drp_suffix
    esp_ea_df = pd.read_csv(result, header=None, names=config.esp_ea_header, \
                             sep='\t', usecols=config.esp_ea_use_header,\
                             low_memory=False)
    #esp_ea_df= esp_ea_df.drop('label', axis=1)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.esp_ea_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.esp_ea_flt_suffix)
    if (len(esp_ea_df.index) == 0):
        esp_ea_df = pd.DataFrame(columns=config.esp_ea_use_header)
    esp_ea_df[config.basic_header] = esp_ea_df[config.basic_header].astype(str)
    return esp_ea_df

def fetch_cg(inputfile, rflag): # return Datsasame    
    ### fetch cg_allele_freq
    print('\n**************************************************************') 
    print('*    Fetch cg allele frequency                               *')
    print('**************************************************************')

    cmd_str= 'annotate_variation.pl  -filter -out ' + inputfile   + ' -dbtype ' \
            + config.cg + ' -build ' + config.buildver + ' ' + inputfile \
            + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.cg_drp_suffix
    cg_df = pd.read_csv(result, header=None, names=config.cg_header, \
                             sep='\t', usecols=config.cg_use_header,\
                             low_memory=False)
    #cg_df= cg_df.drop('label', axis=1)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.cg_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.cg_flt_suffix)
    if (len(cg_df.index) == 0):
        cg_df = pd.DataFrame(columns=config.cg_use_header)
    cg_df[config.basic_header] = cg_df[config.basic_header].astype(str)
    return cg_df

def fetch_popfreq(inputfile, rflag): # return Datsasame    
    ### fetch popfreq_allele_freq
    print('\n**************************************************************') 
    print('*    Fetch popfreq allele frequency                          *')
    print('**************************************************************')
    cmd_str= 'annotate_variation.pl  -filter -out ' + inputfile   + ' -dbtype ' \
            + config.popfreq + ' -build ' + config.buildver + ' ' + inputfile \
            + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.popfreq_drp_suffix
    popfreq_df = pd.read_csv(result, header=None, names=config.popfreq_header, \
                             sep='\t', usecols=config.popfreq_use_header,\
                             low_memory=False)
    #popfreq_df= popfreq_df.drop('label', axis=1)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.popfreq_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.popfreq_flt_suffix)
    if (len(popfreq_df.index) == 0):
        popfreq_df = pd.DataFrame(columns=config.popfreq_use_header)
    popfreq_df[config.basic_header] = popfreq_df[config.basic_header].astype(str)
    return popfreq_df