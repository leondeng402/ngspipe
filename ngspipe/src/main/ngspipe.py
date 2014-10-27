#! /usr/local/bin/python3.4
#####################################################################
#
#    ngspipe.py: a program to annotate next generation sequencing 
#                data at variant,gene and region levels
#    author:     Liyong Deng
#    copyright (c) 2014 by Liyong Deng
#
#####################################################################

import os, time, subprocess, shlex
import argparse

import pandas as pd

import config, gene, variant, region


#####################################################################
#    handle command line arguments
#####################################################################
parser = argparse.ArgumentParser(description='Annotate Next Generation ' + \
                                 'Sequence file')
parser.add_argument('-i', nargs='?', help='input file')
parser.add_argument('-o', nargs='?', help='output directory')
#parser.add_argument('-n', nargs='?', default=5, \
#                    help='number of columns: range[5,6]')
parser.add_argument('-t', action='store_true', help='timing each operation')
parser.add_argument('-e', action='store_true', help='expanding gene entries')
parser.add_argument('-r', action='store_true', help='remove temporary files')
parser.add_argument('-m', action='store_true', help='merge gene info from ' \
                    +'refgene, knowngeen, and ensgene')
args=parser.parse_args()

#print('args.n:', args.n)
#n = int(args.n)
#if(n < 5 or n > 6):
#     print('Error Input n must be in range[5,6]: ' + args.n)
#     exit()

if(not os.path.isfile(args.i)):
    print('Error: File ' + args.i + ' does not exist.')
    exit()
if(not os.path.isdir(args.o)):
    print('Error: Directory ' + args.o + ' does not exist.')
    exit()
infile = args.i
outdir = args.o
if(args.t):
    start_time = time.time()
#####################################################################
#  Validate the input file:
#      1. Expect to take two kinds of files: .txt and .vcf
#      2. validate: chr, start, end, ref, alt
#      3. annovar can also validate it, too.
#      4. the number of columns in an input file
#      5. check whether there is a header
#####################################################################
print('\n**************************************************************')
print('*    Validate the input file                                 *')
print('**************************************************************')

if(infile.endswith('.vcf')):
    outfile = infile + '.anvf'    
    cmd_str = 'convert2annovar.pl -format vcf4 ' + infile + ' -outfile ' \
              + outfile + ' -allsample -withfreq'
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    infile = outfile
    df = pd.read_csv(infile, header=None, names=config.vcfinput_header, \
                     sep='\t', 
                     usecols=config.basic_header, low_memory=False)
else:
    ################################################################# 
    # need to validate the number of columns in an input file
    #################################################################
    fieldnum=0
    for line in open(infile):
        words = line.split('\t')
        fieldnum = len(words)
    if(fieldnum == 5):
        df = pd.read_csv(infile, header=None,  sep='\t', \
                             names=config.basic_header, low_memory=False)
    elif(fieldnum == 6):           
        df = pd.read_csv(infile, header=None,  sep='\t', \
                         names=config.txtinput_header, low_memory=False,\
                         usecols=config.basic_header)
    else:
        print('Error: columns number ' + n + '  is out of range[5,6]')
        #df = comment_df.drop('comments', 1)
    infile = infile + '.anvf'    
    # write cleaned data into infile
df.to_csv(infile, sep='\t', index=False, header=False, low_memory=False)


if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()
#####################################################################
#  fetch dbsnpid
#####################################################################
snp_df = variant.fetch_snpid(infile, args.r)
df = pd.merge(df, snp_df, how='left', on=config.basic_header)

if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()
#make 5 basic columns to be string to have less problem in join
df[config.basic_header] = df[config.basic_header].astype(str)
#####################################################################
#  fetch gene info from refGene, knownGene, and ensGene
#####################################################################
if(args.m):
    #################################################################
    #  have geneinfos, from refgene, knowngene, and ensgene,  merged into one
    #################################################################
    gene_df = gene.fetch_mergedgene(infile, args.r)
    df = pd.merge(df, gene_df, how='left', on=config.basic_header)
    
else:
    #################################################################
    #  three geneinfos from refgene, knowngene, and ensgene
    #################################################################
    gene_df = gene.fetch_refgene(infile, args.r)
    df = pd.merge(df, gene_df, how='left', on=config.basic_header)
                  

    gene_df = gene.fetch_knowngene(infile, args.r)
    df = pd.merge(df, gene_df, how='left', on=config.basic_header)

    gene_df = gene.fetch_ensgene(infile, args.r)
    df = pd.merge(df, gene_df, how='left', on=config.basic_header)


if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()

#####################################################################
#  fetch cytoband info
#####################################################################
cyto_df = region.fetch_cytoband(infile, args.r)
df = pd.merge(df, cyto_df, how='left', on=config.basic_header)


if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()

#####################################################################
#  fetch genomicSuperDups info
#####################################################################
superdup_df = region.fetch_superdup(infile, args.r)
df = pd.merge(df, superdup_df, how='left', on=config.basic_header)

if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()

#####################################################################
#  fetch gwas info
#####################################################################
gwas_df = region.fetch_gwas(infile, args.r)

df = pd.merge(df, gwas_df, how='left', on=config.basic_header)
#print(df.fillna('.'))

if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()

#####################################################################
#  fetch g1k allele frequency
#####################################################################
print('\n**************************************************************') 
print('*    Start to fetch g1k allele frequency: 6 sets             *')
print('**************************************************************')
#g1k_all
g1k_all_df = variant.fetch_g1k_all(infile, args.r)
df = pd.merge(df, g1k_all_df, how='left', on=config.basic_header)

#g1k_afr
g1k_afr_df = variant.fetch_g1k_afr(infile, args.r)
df = pd.merge(df, g1k_afr_df, how='left', on=config.basic_header)

#g1k_amr
g1k_amr_df = variant.fetch_g1k_amr(infile, args.r)
df = pd.merge(df, g1k_amr_df, how='left', on=config.basic_header)

#g1k_eas
g1k_eas_df = variant.fetch_g1k_eas(infile, args.r)
df = pd.merge(df, g1k_eas_df, how='left', on=config.basic_header)

#g1k_eur
g1k_eur_df = variant.fetch_g1k_eur(infile, args.r)
df = pd.merge(df, g1k_eur_df, how='left', on=config.basic_header)

#g1k_sas
g1k_sas_df = variant.fetch_g1k_sas(infile, args.r)
df = pd.merge(df, g1k_sas_df, how='left', on=config.basic_header)

if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()
#####################################################################
#  fetch esp allele frequency 3 sets
#####################################################################
print('\n**************************************************************') 
print('*    Start to fetch ESP allele frequency: 3 sets             *')
print('**************************************************************')
esp_all_df = variant.fetch_esp_all(infile, args.r)
df = pd.merge(df, esp_all_df, how='left', on=config.basic_header)

esp_aa_df = variant.fetch_esp_aa(infile, args.r)
df = pd.merge(df, esp_aa_df, how='left', on=config.basic_header)

esp_ea_df = variant.fetch_esp_ea(infile, args.r)
df = pd.merge(df, esp_ea_df, how='left', on=config.basic_header)

if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()
#####################################################################
#  fetch cg and popfreq allele frequency
#####################################################################
cg_df = variant.fetch_cg(infile, args.r)
df = pd.merge(df, cg_df, how='left', on=config.basic_header)

popfreq_df = variant.fetch_popfreq(infile, args.r)
df = pd.merge(df, popfreq_df, how='left', on=config.basic_header)

if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()
#####################################################################
#  fetch pathogenicity info from ljb
#####################################################################


if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()
#####################################################################
#  fetch cadd score
#####################################################################


if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()
#####################################################################
#  fetch region-based annotations from cytoband, genomicSuperDups, 
#  and gwas 
#####################################################################


if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()
#####################################################################
#  fetch hgmd info
#####################################################################
hgmd_df = variant.fetch_hgmd()
hgmd_df[['chr', 'start', 'end']] = hgmd_df[['chr', 'start', 'end']].astype(str)
df = df.merge(hgmd_df, how='left', on=['chr', 'start', 'end'])

if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()
#print(df.fillna('.').to_string())
#####################################################################
#  expand multiple entries in gene info
#####################################################################
df = df.fillna('.')
if(args.o.endswith('/')):
    file = args.o + '' + infile.split('/')[-1] + '.annotated'
else:
    file = args.o + '/' + infile.split('/')[-1] + '.annotated'
if(args.e):
    #expanding gene entries
    gene.expand_gene_entries(df, file)    
    if(args.t):
        end_time = time.time()
        print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
else:
    df.to_csv(file, sep='\t', index=False)




