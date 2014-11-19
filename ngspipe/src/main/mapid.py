#! /usr/local/bin/python3.4
#####################################################################
#
#    mapid.py: a program to map ids in vcf file
#    author:     Liyong Deng
#    copyright (c) 2014 by Liyong Deng in Dr. Chung's Lab
#
#####################################################################

import os, time
import argparse, re
import pandas as pd

parser = argparse.ArgumentParser(description='map the ids in vcf file')
parser.add_argument('-i', nargs='?', help='input file')
parser.add_argument('-o', nargs='?', help='output file name')
parser.add_argument('-f', action='store_true', \
                    help='force to overwrite existing file')
parser.add_argument('-m', nargs='+', help='apply id file with same id order')

args=parser.parse_args()

if(args.i == None or args.o == None):
    print('Error: -i and -o is required')
    exit()
if(not os.path.isfile(args.i)):
    print('Error: File ' + args.i + ' does not exist.')
    exit()
if(os.path.isfile(args.o) & (not args.f)):
    print('Error: File ' + args.o + '  already exists.')
    exit()
input = args.i    

if(not input.endswith('.vcf')):
        exit(0)
of = open(args.o, 'w')
with open(input, 'r')as f:
    results=[]
    for line in f.readlines():
        line = line.strip()            
        if(line.startswith('#CHROM')):            
            words=line.split('\t') 
            c = 0
            #print(len(words))             
            for word in words:
                if(c <= 8):
                    of.write(word)    
                else:
                    items = word.split('_')
                    _digit=re.compile('\d')
                    if(not bool(_digit.search(items[1]))):
                        of.write(items[1]+items[2])
                    else:
                        of.write(items[1])    
                if(c != (len(words) - 1)):
                    #print(c)
                    of.write('\t')  
                else:
                    #print(c)
                    of.write('\n') 
                                    
                c= c+1               
                
        else:
            of.write(line+'\n')
of.close()
            
                
