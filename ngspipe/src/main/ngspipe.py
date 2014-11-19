#! /usr/local/bin/python3.4
#####################################################################
#
#    ngspipe.py: a program to annotate next generation sequencing 
#                data at variant,gene and region levels
#    author:     Liyong Deng
#    copyright (c) 2014 by Liyong Deng in Dr. Chung's Lab
#
#####################################################################

import os, time, subprocess, shlex
import argparse

import pandas as pd

import config, gene, variant, region, operation, ino


#####################################################################
#    handle command line arguments
#####################################################################
parser = argparse.ArgumentParser(description='Annotate Next Generation ' + \
                                 'Sequence file')
parser.add_argument('-c', action='store_true', help='consolidate gene info ' \
                    + 'from refgene, knowngene, and ensgene')
parser.add_argument('-d', action='store_true', help='delete temporary files')
parser.add_argument('-e', action='store_true', help='expanding gene entries')
parser.add_argument('-f', nargs='+', help='apply filter in the file')
parser.add_argument('-i', nargs='?', help='input file')
parser.add_argument('-o', nargs='?', help='output directory')
parser.add_argument('-r', nargs='+', help='output ranked list')
parser.add_argument('-t', action='store_true', help='timing each operation')
parser.add_argument('-m', nargs='?', help='proband/trios/family '+\
                     'information file')
parser.add_argument('--genetic_model', nargs='?', choices=['denovo', 'recessive', \
                                                  'xlinked', 'other'],\
                    default='denovo', help='choose the inheritance models')

args=parser.parse_args()

if(args.m != None):
    print('Validate the pedigree file')
    #validate args.m file exists with correct format
    
if(args.i == None or args.o == None):
    print('Error: -i and -o is required')
    exit()

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
vcf_flag = False
if(infile.endswith('.vcf')):
    vcf_flag = True
    path = os.path.dirname(infile)+'/'
    base = os.path.basename(infile)
    # fetch the header of chr, start, end, ref, alt and all sample ID
    with open(infile, 'r')as f:        
        for line in f.readlines():
            line = line.strip()    
            if(line.startswith('#CHROM')):            
                line = line.replace('#', '')
                words=line.split('\t')                
                sampleid=words[9:]
                break
    
    md_str = 'convert2annovar.pl  -format vcf4 -allsample -withfreq -include '\
             + infile + ' -outfile ' + path+base.replace('.vcf', '.temp')
    cmd = shlex.split(md_str)
    #subprocess.call(cmd)
    p = subprocess.Popen(cmd)
    p.wait()    # genotype_df as genotype dataframe
    geno_df = operation.retrieve_genotype(path+base.replace('.vcf', '.temp'), \
                                          sampleid)
    os.remove(path+base.replace('.vcf', '.temp'))
    infile = path+base.replace('.vcf', '.txt')
    
    geno_df.to_csv(infile, header=False, sep='\t', index=False, \
                   columns=config.basic_header)    
    geno_df.to_csv(infile.replace('.txt', '.geno'), sep='\t', index=False)
    fieldnum=len(config.basic_header)
else:
    ################################################################# 
    # need to validate the number of columns in an input file
    #################################################################
    fieldnum=0
    for line in open(infile):
        words = line.split('\t')
        fieldnum = len(words)
        break
    
    if(fieldnum not in range(5, 7)):    
        print('Error: columns number ' + str(fieldnum) + '  is out of range[5,6]')
        exit()

if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()
#####################################################################
#  table_annovar queries in one
#####################################################################
print('\n**************************************************************')
print('*    multiple database annotations                           *')
print('**************************************************************')
md_str= 'table_annovar.pl ' + infile + ' ' \
        + config.dir_humandb \
        + ' -buildver ' + config.buildver \
        + ' -protocol ' + ','.join(config.protocols) \
        + ' -operation ' + ','.join(config.operations)\
        + ' -nastring ' + config.nastring \
        + ' -csvout'
   
if(args.d):
    md_str=md_str + ' -remove'
#print(md_str)
cmd = shlex.split(md_str)
subprocess.call(cmd)
result = infile+'.'+config.buildver + config.table_anno_suffix
            
header=config.table_anno_header
df = pd.read_csv(result, header=None, names=header, sep=',', \
                 low_memory=False, skiprows=1)
df[config.basic_header] = df[config.basic_header].astype(str)
if(args.d):   
    os.remove(result)
if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()

#####################################################################
#  fetch 1000 genome project all
#####################################################################
for race in config.g1k_races:
    g1k_df = variant.fetch_g1k(infile, race, fieldnum, vcf_flag, args.d)
    #print(g1k_df)
    df = df.merge(g1k_df, how='left', on=config.basic_header)
   

if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()
#####################################################################
#  fetch cadd score
#####################################################################
cadd_df = variant.fetch_cadd(infile, fieldnum, vcf_flag, args.d)
#print(cadd_df)
df = pd.merge(df, cadd_df, how='left', on=config.basic_header)

if(args.d):
    os.remove(infile+'.log')
    
if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()
#####################################################################
#  fetch hgmd info
#####################################################################
hgmd_df = variant.fetch_hgmd()
df = df.merge(hgmd_df, how='left', on=['Chr', 'Start', 'End'])

if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()
#print(df.columns.values)

# keep only display columns
df=pd.DataFrame(df, columns=config.custom_header)
df=df.fillna('.')


if(args.o.endswith('/')):
    outdir = args.o  
    file = outdir + infile.split('/')[-1]  + '.annotated'
else:
    outdir = args.o + '/' 
    file = outdir + infile.split('/')[-1] + '.annotated'

#####################################################################
#  consolidate gene entries from different sources
#####################################################################
if(args.c):
    # merging gene entries
    df = gene.merge_gene_entries(df)
        
    df.to_csv(file.replace('.annotated', '.consilidated'), sep='\t', index=False)
    if(args.t):
        end_time = time.time()
        print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
#####################################################################
#  expand multiple entries in gene info
#####################################################################
if(args.e):
    # expanding gene entries
    df = gene.expand_gene_entries(df)    
    df.to_csv(file.replace('.annotated', '.expanded'), sep='\t', index=False)
    if(args.t):
        end_time = time.time()
        print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    

if(vcf_flag):
    ino.output_annotation(file, df, geno_df, 'singleton', args.m)
    ino.output_annotation(file, df, geno_df, 'trios', args.m)
    ino.output_annotation(file, df, geno_df, 'family', args.m)
    #df.to_csv(file, sep='\t', index=False)
else:
    df.to_csv(file, sep='\t', index=False)
    
#####################################################################
#  apply general filters
#####################################################################
print('\n**************************************************************')
print('*    start to apply filter(s)                                *')
print('**************************************************************')
if(args.f != None):
    for filter in args.f:
        df = operation.apply_filter(df, filter)

    df.to_csv(file.replace('.annotated', '.filtered'), sep='\t', index=False)

#####################################################################
#  apply inheritance filter
##################################################################### 
print('\n**************************************************************')
print('*    start to apply inheritance filter(s)                    *')
print('**************************************************************')
if(args.genetic_model=='denovo'):
    ino.output_denovo(file, df, geno_df, 'singleton', args.m)
    ino.output_denovo(file, df, geno_df, 'trios', args.m)
    ino.output_denovo(file, df, geno_df, 'family', args.m)
elif(args.genetic_model=='recessive'):
    ino.output_recessive(file, df, geno_df, 'singleton', args.m)
    ino.output_recessive(file, df, geno_df, 'trios', args.m)
    ino.output_recessive(file, df, geno_df, 'family', args.m)
elif(args.genetic_model=='xlinked'):
    ino.output_xlinked(file, df, geno_df, 'singleton', args.m)
    ino.output_xlinked(file, df, geno_df, 'trios', args.m)
    ino.output_xlinked(file, df, geno_df, 'family', args.m)
    
elif(args.genetic_model=='other'):
    pass
else:
    print("Error: genetic model is not implemented")


#####################################################################
#  rank variants/genes
##################################################################### 
print('\n**************************************************************')
print('*    start to rank variants/genes                            *')
print('**************************************************************') 
if(args.r != None):
    for ranker in args.r:
        df = operation.rank(df, ranker)  


