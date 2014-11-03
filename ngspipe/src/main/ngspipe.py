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
if(vcf_flag):
    md_str= md_str + ' -vcfinput'
    
if(args.r):
    md_str=md_str + ' -remove'
#print(md_str)
cmd = shlex.split(md_str)
subprocess.call(cmd)
result = infile+'.'+config.buildver + config.table_anno_suffix
            
#print(result)
header=config.table_anno_header
df = pd.read_csv(result, header=None, names=header, sep=',', \
                 usecols=header[:26]+header[32:],low_memory=False, skiprows=1)
#print(df.columns.values)

df.columns=header[:26]+header[32:]
if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()

#print(df)
#####################################################################
#  fetch 1000 genome project all
#####################################################################
for race in config.g1k_races:
    g1k_df = variant.fetch_g1k(infile, race, vcf_flag, args.r )
    df = df.merge(g1k_df, how='left', on=config.basic_header)
    
if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()
#####################################################################
#  fetch cadd score
#####################################################################
cadd_df = variant.fetch_cadd(infile, vcf_flag, args.r)
df = pd.merge(df, cadd_df, how='left', on=config.basic_header)


if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()
#####################################################################
#  fetch hgmd info
#####################################################################
hgmd_df = variant.fetch_hgmd()
print(hgmd_df.columns.values)
df = df.merge(hgmd_df, how='left', on=['Chr', 'Start', 'End'])

if(args.t):
    end_time = time.time()
    print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
    start_time = time.time()
print(df.columns.values)
#####################################################################
#  expand multiple entries in gene info
#####################################################################
#df = df.fillna('.')
if(args.o.endswith('/')):
    pass
    file = args.o + '' + infile.split('/')[-1] + '.annotated'
else:
    pass
    file = args.o + '/' + infile.split('/')[-1] + '.annotated'
if(args.e):
    pass
    #expanding gene entries
    gene.expand_gene_entries(df, file)    
    if(args.t):
        end_time = time.time()
        print(time.strftime('%H:%M:%S', time.gmtime(end_time-start_time)))
else:
    pass
df.fillna('.')
df.to_csv(file, sep='\t', index=False)




