#! /usr/local/bin/python3.4
#####################################################################
#
#    variant.py:  the file contains all the variant level functions  
#                 related to annovar
#    author:  Liyong Deng
#    copyright (c) 2014 by Liyong Deng in Dr. Chung's Lab
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
    hgmd_df[['Chr', 'Start', 'End']] = hgmd_df[['Chr', 'Start', 'End']].astype(str)
    return hgmd_df
    
def fetch_g1k(inputfile, race, vflag, rflag): # return DataFrame    
    ### fetch g1k_all_allele_freq
    print('\nFetch g1k_' + race + ' allele frequency')
    cmd_str= 'annotate_variation.pl  -filter -out ' + inputfile  \
             + ' -dbtype ' + config.g1k + '_' + race \
             +' -build ' + config.buildver + ' '\
             + inputfile  + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+'_' + race.upper() \
             + config.g1k_drp_suffix
    if(vflag):
        header = config.g1k_header
    else:
        header = config.g1k_header[:2] + config.txtinput_header
    #print(result)
    #print(header)   
    g1k_df = pd.read_csv(result, header=None, names=header, \
                             sep='\t', usecols=header[1:],\
                             low_memory=False)
    g1k_df = g1k_df.rename(columns={'Value': config.g1k+'_'+race})
    if(rflag):
        os.remove(inputfile+'.'+config.buildver + '_' + race.upper() \
                  + config.g1k_drp_suffix)
        os.remove(inputfile+'.'+config.buildver + '_' + race.upper() \
                  +config.g1k_flt_suffix)
    if (len(g1k_df.index) == 0):
        g1k_df = pd.DataFrame(columns=config.g1k_header[1:])
        g1k_df = g1k_df.rename(columns={'Value': 'g1k_'+race})
    if 'Comments' in g1k_df:
        g1k_df= g1k_df.drop('Comments', axis=1)
    g1k_df[config.basic_header] = g1k_df[config.basic_header].astype(str)
    #print(g1k_df.columns.values)
    return g1k_df

def fetch_g1k_afr(inputfile, rflag): # return DataFrame    
    ### fetch g1k_afr_allele_freq
    print('\nFetch g1k_afr allele frequency')
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
    print('\nFetch g1k_amr allele frequency')
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
    print('\nFetch g1k_eas allele frequency')
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
    print('\nFetch g1k_eur allele frequency')
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
    print('\nFetch g1k_sas allele frequency')
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
    print('\nFetch esp_all allele frequency')
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
    print('\nFetch esp_aa allele frequency')
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
    print('\nFetch esp_ea allele frequency')
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

def fetch_ljb(inputfile, rflag): # return Datsasame    
    ### fetch popfreq_allele_freq
    print('\n**************************************************************') 
    print('*    Fetch LJB pathogenicity scores                          *')
    print('**************************************************************')
    cmd_str= 'annotate_variation.pl  -filter -otherinfo -out ' + inputfile   \
             + ' -dbtype ' + config.ljb + ' -build ' + config.buildver + ' ' \
             + inputfile  + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.ljb_drp_suffix
    of = open(result+'.temp', 'w')
    with open(result, 'r') as f:
        data = f.read().replace('\t', ',')
        of.write(data)    
    of.close()
    os.remove(result)
    os.rename(result+'.temp', result)
    ljb_df = pd.read_csv(result, header=None, names=config.ljb_header, \
                             sep=',', usecols=config.ljb_use_header,\
                             low_memory=False)
    #ljb_df= ljb_df.drop('label', axis=1)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.ljb_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.ljb_flt_suffix)
    if (len(ljb_df.index) == 0):
        ljb_df = pd.DataFrame(columns=config.ljb_use_header)
    ljb_df[config.basic_header] = ljb_df[config.basic_header].astype(str)
    return ljb_df
    
    
def fetch_cadd(inputfile, vflag, rflag): # return Datsasame    
    ### fetch popfreq_allele_freq
    print('\n**************************************************************') 
    print('*    Fetch cadd pathogenicity scores                         *')
    print('**************************************************************')
    cmd_str= 'annotate_variation.pl  -filter -otherinfo -out ' + inputfile   \
             + ' -dbtype ' + config.cadd + ' -build ' + config.buildver + ' ' \
             + inputfile  + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.cadd_drp_suffix
    of = open(result+'.temp', 'w')
    with open(result, 'r') as f:
        data = f.read().replace('\t', ',')
        of.write(data)    
    of.close()
    os.remove(result)
    os.rename(result+'.temp', result)
    if(vflag):
        header = config.cadd_header
    else:
        header = config.cadd_header[:3] + config.txtinput_header
    print(header)
    print(header[1:])
    cadd_df = pd.read_csv(result, header=None, names=header, \
                             sep=',', usecols=header[1:],\
                             low_memory=False)
    if 'Comments' in cadd_df:
        cadd_df= cadd_df.drop('Comments', axis=1)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.cadd_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.cadd_flt_suffix)
    if (len(cadd_df.index) == 0):
        cadd_df = pd.DataFrame(columns=config.cadd_use_header)
    cadd_df[config.basic_header] = cadd_df[config.basic_header].astype(str)
    return cadd_df
    

def fetch_clinvar(inputfile, rflag): # return Datsasame    
    ### fetch clinvar_allele_freq
    print('\n**************************************************************') 
    print('*    Fetch clinvar info                                      *')
    print('**************************************************************')
    cmd_str= 'annotate_variation.pl  -filter -out ' + inputfile   + ' -dbtype ' \
            + config.clinvar + ' -build ' + config.buildver + ' ' + inputfile \
            + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.clinvar_drp_suffix
    clinvar_df = pd.read_csv(result, header=None, names=config.clinvar_header, \
                             sep='\t', usecols=config.clinvar_use_header,\
                             low_memory=False)
    #clinvar_df= clinvar_df.drop('label', axis=1)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.clinvar_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.clinvar_flt_suffix)
    if (len(clinvar_df.index) == 0):
        clinvar_df = pd.DataFrame(columns=config.clinvar_use_header)
    clinvar_df[config.basic_header] = clinvar_df[config.basic_header].astype(str)
    return clinvar_df
    
def fetch_cosmic(inputfile, rflag): # return Datsasame    
    ### fetch cosmic_allele_freq
    print('\n**************************************************************') 
    print('*    Fetch cosmic info                                      *')
    print('**************************************************************')
    cmd_str= 'annotate_variation.pl  -filter -out ' + inputfile   + ' -dbtype ' \
            + config.cosmic + ' -build ' + config.buildver + ' ' + inputfile \
            + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.cosmic_drp_suffix
    cosmic_df = pd.read_csv(result, header=None, names=config.cosmic_header, \
                             sep='\t', usecols=config.cosmic_use_header,\
                             low_memory=False)
    #cosmic_df= cosmic_df.drop('label', axis=1)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.cosmic_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.cosmic_flt_suffix)
    if (len(cosmic_df.index) == 0):
        cosmic_df = pd.DataFrame(columns=config.cosmic_use_header)
    cosmic_df[config.basic_header] = cosmic_df[config.basic_header].astype(str)
    return cosmic_df
      
    
    
def fetch_nci(inputfile, rflag): # return Datsasame    
    ### fetch nci_allele_freq
    print('\n**************************************************************') 
    print('*    Fetch nci info                                      *')
    print('**************************************************************')
    cmd_str= 'annotate_variation.pl  -filter -out ' + inputfile   + ' -dbtype ' \
            + config.nci + ' -build ' + config.buildver + ' ' + inputfile \
            + ' ' + config.dir_humandb
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    result = inputfile+'.'+config.buildver+config.nci_drp_suffix
    nci_df = pd.read_csv(result, header=None, names=config.nci_header, \
                             sep='\t', usecols=config.nci_use_header,\
                             low_memory=False)
    #nci_df= nci_df.drop('label', axis=1)
    if(rflag):
        os.remove(inputfile+'.'+config.buildver+config.nci_drp_suffix)
        os.remove(inputfile+'.'+config.buildver+config.nci_flt_suffix)
    if (len(nci_df.index) == 0):
        nci_df = pd.DataFrame(columns=config.nci_use_header)
    nci_df[config.basic_header] = nci_df[config.basic_header].astype(str)
    return nci_df
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    