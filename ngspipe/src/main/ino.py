#! /usr/local/bin/python3.4
#####################################################################
#
#    ino.py:  the file input and output operations
#    author:  Liyong Deng
#    copyright (c) 2014 by Liyong Deng
#
#####################################################################
import os,re
import pandas as pd
import config


def output_annotation(output, dataframe, genotype_df, mode, pedf): 
    # parse pedigree information and build a pedigree dataframe
    print('parse the pedigree file')
    print(mode)
    if(mode != 'mono'):
        ped_df = pd.read_csv(pedf, header=0,  sep='\t', low_memory=False)
        
    # get the output directory and file
    dir = os.path.dirname(output)+'/annotated'
    if (not os.path.exists(dir)):
        os.mkdir(dir)
    print(dir)
    filename = os.path.basename(output)
    
    if(mode == 'mono'):
        print('output the results in a single file')
        df = pd.merge(variant_df, genotype_df, how='left', \
                      on=['Chr', 'Start', 'End', 'Ref', 'Alt'])
        df.to_csv(dir+'/'+filename, sep='\t', index=False)
        
    elif(mode == 'singleton'):
        print('output the results of singleton only ')
        probanddir=dir+'/singletons'
        if (not os.path.exists(probanddir)):
            os.makedirs(probanddir)
        proband= ped_df[ped_df['Person']=='Proband']['SampleID'].tolist()
        for sampleid in proband:
            sample_df = genotype_df[config.basic_header+[sampleid]]
            sample_df = sample_df[sample_df[sampleid]!='0/0']
            sample_df = sample_df[sample_df[sampleid]!='0|0']
            sample_df = sample_df[sample_df[sampleid]!='./.']
            #print(dataframe.columns.values)
            combined_df = pd.merge(sample_df, dataframe,  how='left', \
                           on=config.basic_header)
            familyid = ped_df[ped_df['SampleID']==sampleid]['FamilyID'].tolist()
            family_df = ped_df[ped_df['FamilyID']==familyid[0]]
            if(len(family_df.index) == 1):
                combined_df.to_csv(probanddir+'/'+sampleid+'.annotated', \
                                   sep='\t', index=False)
        # make a proband directory
        
    elif(mode == 'trios'):
        print('output the results of trios and singleton')
        triodir = dir+'/trios'
        if (not os.path.exists(triodir)):
            os.makedirs(triodir)
        familyid_list = ped_df[ped_df['Person']=='Proband']\
        ['FamilyID'].tolist()
        familyid_list = [id for id in familyid_list if id != 'Unknown']
        for familyid in familyid_list:
            #print(ped_df.columns.values)
            temp_df = ped_df[ped_df['FamilyID']==familyid]
            persons=temp_df['Person'].tolist()
            sampleids=[]
            if('Proband' in persons and 'Mother' in persons \
               and 'Father' in persons):
                print('This is a trio')
                trios=['Proband', 'Father', 'Mother']
                for person in trios:
                    item = temp_df[temp_df['Person']==person]\
                    ['SampleID'].tolist()                    
                    sampleids.append(item[0])
                print(sampleids)
                sample_df = genotype_df[config.basic_header+sampleids]
                # filter rows with all 0/0                  
                # iterate rows to filter out ./., 0/0 and 0|0
                retained_genotype=[]
                for index, row in sample_df.iterrows():
                    flag = 0
                    #print('index: ', index)
                    for id in sampleids:
                        sampleid_person = temp_df[temp_df['SampleID']==id]['Person'].tolist()                                                 
                        if(row[id] == './.' or row[id] == '0/0' 
                           or row[id] == '0|0'):
                            flag = 0
                            if(sampleid_person[0] == 'Proband'):
                                break
                        else:
                            flag=1
                        
                    if(flag > 0):
                        retained_genotype.append(row)
                header = sample_df.columns.values
                if(len(retained_genotype) > 0):
                    print('Trios ', sampleids, ' wirte genotypes to files')
                    sample_df = pd.DataFrame(retained_genotype, columns=header,\
                                         index=None)
                    combined_df = pd.merge(sample_df, dataframe, how='left', \
                                           on=config.basic_header)
                    combined_df.to_csv(triodir+'/'+familyid+'.annotated', \
                                       sep='\t', index=False)
                else:
                    print('Trios ', sampleids, ' have no genotypes left')
            
    elif(mode == 'family'):
        print('output the results of family and singleton')
        familydir = dir+'/families'
        if (not os.path.exists(familydir)):
            os.makedirs(familydir)
        familyid_list = ped_df[ped_df['Person']=='Proband']['FamilyID'].tolist()
        familyid_list = [id for id in familyid_list if id != 'Unknown']
        for familyid in familyid_list:
            persons = ped_df[ped_df['FamilyID']==familyid]['Person'].tolist()
            #print(persons)
            if(len(persons) >= 3):
                if( not ('Proband' in persons and 'Mother' in persons \
                         and 'Father' in persons)):
                    print('This is not  a trio')
                    print('This is a family')
            else:
                # either a family or single family member
                print('This is a family or singleton')
    else:
        print('Error: ' + mode +' wrong mode')
    
def output_denovo(output, dataframe, genotype_df, mode, pedf): 
    
    print('Start to write denovo filtering results ...')
    ## add denovo_maf
    print('\nmaf filtering')
    maf = config.denovo_maf
    header = dataframe.columns.values
    nadf = dataframe[dataframe['PopFreqMax'] == '.']
    dataframe = dataframe[dataframe['PopFreqMax'] != '.']
    dataframe.is_copy = False
    dataframe[['PopFreqMax']] = dataframe[['PopFreqMax']].astype(float)
    dataframe = dataframe[dataframe['PopFreqMax'] <= maf]
    dataframe = dataframe.append(nadf)
    if(len(dataframe.index)==0):
        print('None DataFrame')
        dataframe=None
    # parse the pedigree file
    print('parse the pedigree file')
    if(mode != 'mono'):
        ped_df = pd.read_csv(pedf, header=0,  sep='\t', low_memory=False)
        
    # get the output directory and file
    dir = os.path.dirname(output)+'/denovo'
    if (not os.path.exists(dir)):
        os.mkdir(dir)
    print(dir)
    filename = os.path.basename(output)
    
    if(mode == 'mono'):
        print('output the results in a single file')
        df = pd.merge(variant_df, genotype_df, how='left', \
                      on=['Chr', 'Start', 'End', 'Ref', 'Alt'])
        df.to_csv(dir+'/'+filename.replace('.annotated', '.denovo.txt'), sep='\t', index=False)
        
    elif(mode == 'singleton'):
        print('output the results of singleton only ')
        probanddir=dir+'/singletons'
        if (not os.path.exists(probanddir)):
            os.makedirs(probanddir)
        proband= ped_df[ped_df['Person']=='Proband']['SampleID'].tolist()
        for sampleid in proband:
            sample_df = genotype_df[config.basic_header+[sampleid]]
            sample_df = sample_df[sample_df[sampleid]!='0/0']
            sample_df = sample_df[sample_df[sampleid]!='0|0']
            sample_df = sample_df[sample_df[sampleid]!='./.']
            sample_df = sample_df[sample_df[sampleid]!='1/1']
            sample_df = sample_df[sample_df[sampleid]!='2/2']
            sample_df = sample_df[sample_df[sampleid]!='3/3']
            sample_df = sample_df[sample_df[sampleid]!='4/4']
            sample_df = sample_df[sample_df[sampleid]!='5/5']
            #print(dataframe.columns.values)
            combined_df = pd.merge(dataframe, sample_df, how='left', \
                           on=config.basic_header)
            familyid = ped_df[ped_df['SampleID']==sampleid]['FamilyID'].tolist()
            family_df = ped_df[ped_df['FamilyID']==familyid[0]]
            if(len(family_df.index) == 1):
                combined_df.to_csv(probanddir+'/'+sampleid+'.denovo.txt', \
                                   sep='\t', index=False)
        # make a proband directory
        
    elif(mode == 'trios'):
        print('output the results of trios')
        triodir = dir+'/trios'
        if (not os.path.exists(triodir)):
            os.makedirs(triodir)
        familyid_list = ped_df[ped_df['Person']=='Proband']\
        ['FamilyID'].tolist()
        familyid_list = [id for id in familyid_list if id != 'Unknown']
        for familyid in familyid_list:
            #print(ped_df.columns.values)
            temp_df = ped_df[ped_df['FamilyID']==familyid]
            persons=temp_df['Person'].tolist()
            sampleids=[]
            if('Proband' in persons and 'Mother' in persons \
               and 'Father' in persons):
                print('This is a trio')
                trios=['Proband', 'Father', 'Mother']
                for person in trios:
                    item = temp_df[temp_df['Person']==person]\
                    ['SampleID'].tolist()                    
                    sampleids.append(item[0])
                print(sampleids)
                sample_df = genotype_df[config.basic_header+sampleids]
                # filter rows with all 0/0                  
                # iterate rows to filter out ./., 0/0 and 0|0
                print('Filter homo wildtype')
                retained_genotype=[]
                for index, row in sample_df.iterrows():
                    flag = 0
                    #print('index: ', index)
                    father_genotypes=[]
                    mother_genotypes=[]
                    proband_genotypes=[]
                    for id in sampleids:
                        sampleid_person = temp_df[temp_df['SampleID']==id]\
                        ['Person'].tolist() 
                        # filter the rows if trios all have ./., 0/0 or 0|0
                        # Also, fiter the rows if proband genotype is ./., 
                        #0/0, 0|0                                                
                        if(row[id] == './.' or row[id] == '0/0' 
                           or row[id] == '0|0'):
                            flag = 0
                            if(sampleid_person[0] == 'Proband'):
                                break
                        else:
                            flag=1
                        
                        if(sampleid_person[0] == 'Father'):
                            
                            father_genotypes=re.split('/|', row[id])
                        elif(sampleid_person[0] == 'Mother'):
                            
                            mother_genotypes=re.split('/|', row[id])
                        elif(sampleid_person[0] == 'Proband'):
                           
                            proband_genotypes=re.split('/|', row[id])
                    # fitler the proand hets inherited from either mother
                    # or fatehr
                    parent_genotypes=father_genotypes + mother_genotypes
                    #print('proband_genotypes:', proband_genotypes)
                    #print('parent_genotypes:', parent_genotypes)
                    for one_genotype in proband_genotypes:
                        if(one_genotype != '0' \
                           and one_genotype in parent_genotypes):
                            #print('one_genotype', one_genotype)                            
                            flag = 0
                    # check the zygosity
                    if(len(proband_genotypes) > len(set(proband_genotypes))):
                        flag = 0
                    #print('flag:', flag)                        
                    if(flag == 1):
                        #print(sampleids)
                        retained_genotype.append(row)
                print(sampleids)
                header = sample_df.columns.values
                if(len(retained_genotype) > 0):
                    print('Trios ', sampleids, ' wirte genotypes to files')
                    sample_df = pd.DataFrame(retained_genotype, columns=header,\
                                         index=None)
                    combined_df = pd.merge(sample_df, dataframe, how='inner', \
                                           on=config.basic_header)
                    combined_df.to_csv(triodir+'/'+familyid+'.denovo.txt', \
                                       sep='\t', index=False)
                else:
                    print('Trios ', sampleids, ' have no genotypes left')
            
    elif(mode == 'family'):
        print('output the results of family and singleton')
        familydir = dir+'/families'
        if (not os.path.exists(familydir)):
            os.makedirs(familydir)
        familyid_list = ped_df[ped_df['Person']=='Proband']['FamilyID'].tolist()
        familyid_list = [id for id in familyid_list if id != 'Unknown']
        for familyid in familyid_list:
            persons = ped_df[ped_df['FamilyID']==familyid]['Person'].tolist()
            #print(persons)
            if(len(persons) >= 3):
                if( 'Proband' in persons and 'Mother' in persons \
                         and 'Father' in persons):
                    print('This is a trio in a family')
                else:
                    print('This is a family without a trio')
            elif(len(persons) == 2):               
                print('This is a family with two members')
            else:
                if('Proband' not in persons):
                    print('This is a family with single non-proband member')
    else:
        print('Error: ' + mode +' wrong mode')
    

def output_recessive(output, dataframe, genotype_df, mode, pedf): 
    print('Start to write recessive filtering results ...')
    ## add denovo_maf
    print('\nmaf filtering')
    maf = config.recessive_maf
    header = dataframe.columns.values
    nadf = dataframe[dataframe['PopFreqMax'] == '.']
    dataframe = dataframe[dataframe['PopFreqMax'] != '.']
    dataframe.is_copy = False
    dataframe[['PopFreqMax']] = dataframe[['PopFreqMax']].astype(float)
    dataframe = dataframe[dataframe['PopFreqMax'] <= maf]
    dataframe = dataframe.append(nadf)
    if(len(dataframe.index)==0):
        print('None DataFrame')
        dataframe=None
    # parse the pedigree file
    print('parse the pedigree file')
    if(mode != 'mono'):
        ped_df = pd.read_csv(pedf, header=0,  sep='\t', low_memory=False)
        
    # get the output directory and file
    dir = os.path.dirname(output)+'/recessive'
    if (not os.path.exists(dir)):
        os.mkdir(dir)
    print(dir)
    filename = os.path.basename(output)
    
    if(mode == 'mono'):
        print('output the results in a single file')
        df = pd.merge(variant_df, genotype_df, how='left', \
                      on=['Chr', 'Start', 'End', 'Ref', 'Alt'])
        df.to_csv(dir+'/'+filename.replace('.annotated', '.denovo.txt'), sep='\t', index=False)
        
    elif(mode == 'singleton'):
        print('output the results of singleton only ')
        probanddir=dir+'/singletons'
        if (not os.path.exists(probanddir)):
            os.makedirs(probanddir)
        proband= ped_df[ped_df['Person']=='Proband']['SampleID'].tolist()
        for sampleid in proband:
            sample_df = genotype_df[config.basic_header+[sampleid]]
            sample_df = sample_df[sample_df[sampleid]!='0/0']
            sample_df = sample_df[sample_df[sampleid]!='0|0']
            sample_df = sample_df[sample_df[sampleid]!='./.']
            #print(dataframe.columns.values)
            combined_df = pd.merge(dataframe, sample_df, how='left', \
                           on=config.basic_header)
            familyid = ped_df[ped_df['SampleID']==sampleid]['FamilyID'].tolist()
            family_df = ped_df[ped_df['FamilyID']==familyid[0]]
            if(len(family_df.index) == 1):
                combined_df.to_csv(probanddir+'/'+sampleid+'.denovo.txt', \
                                   sep='\t', index=False)
        # make a proband directory
        
    elif(mode == 'trios'):
        print('output the results of trios and singleton')
        triodir = dir+'/trios'
        if (not os.path.exists(triodir)):
            os.makedirs(triodir)
        familyid_list = ped_df[ped_df['Person']=='Proband']\
        ['FamilyID'].tolist()
        familyid_list = [id for id in familyid_list if id != 'Unknown']
        for familyid in familyid_list:
            #print(ped_df.columns.values)
            temp_df = ped_df[ped_df['FamilyID']==familyid]
            persons=temp_df['Person'].tolist()
            sampleids=[]
            if('Proband' in persons and 'Mother' in persons \
               and 'Father' in persons):
                print('This is a trio')
                trios=['Proband', 'Father', 'Mother']
                for person in trios:
                    item = temp_df[temp_df['Person']==person]\
                    ['SampleID'].tolist()                    
                    sampleids.append(item[0])
                print(sampleids)
                sample_df = genotype_df[config.basic_header+sampleids]
                # filter rows with all 0/0                  
                # iterate rows to filter out ./., 0/0 and 0|0
                retained_genotype=[]
                for index, row in sample_df.iterrows():
                    flag = 0
                    #print('index: ', index)
                    father_genotypes=[]
                    mother_genotypes=[]
                    proband_genotypes=[]
                    for id in sampleids:
                        sampleid_person = temp_df[temp_df['SampleID']==id]\
                        ['Person'].tolist() 
                        # filter the rows if trios all have ./., 0/0 or 0|0
                        # Also, fiter the rows if proband genotype is ./., 
                        #0/0, 0|0                                                
                        if(row[id] == './.' or row[id] == '0/0' 
                           or row[id] == '0|0'):
                            flag = 0
                            if(sampleid_person[0] == 'Proband'):
                                break
                        else:
                            flag=1
                        
                        if(sampleid_person[0] == 'Father'):                            
                            father_genotypes=re.split('/|', row[id])
                        elif(sampleid_person[0] == 'Mother'):
                            
                            mother_genotypes=re.split('/|', row[id])
                        elif(sampleid_person[0] == 'Proband'):
                           
                            proband_genotypes=re.split('/|', row[id])
                    # fitler the proand hets inherited from either mother
                    # or fatehr
                    parent_genotypes=father_genotypes + mother_genotypes
                    #print('proband_genotypes:', proband_genotypes)
                    #print('parent_genotypes:', parent_genotypes)
                    for one_genotype in proband_genotypes:
                        if(one_genotype != '0' \
                           and one_genotype not in father_genotypes):
                            #print('one_genotype', one_genotype)                            
                            flag = 0
                        if(one_genotype != '0' \
                           and one_genotype not in mother_genotypes):
                            #print('one_genotype', one_genotype)                            
                            flag = 0
                    # check the zygosity, did not deal with intrans
                    
                    if(len(proband_genotypes) != len(set(proband_genotypes))):
                        
                        flag = 0
                    #print('flag:', flag)                        
                    if(flag == 1):
                        #print(sampleids)
                        retained_genotype.append(row)
                header = sample_df.columns.values
                if(len(retained_genotype) > 0):
                    print('Trios ', sampleids, ' wirte genotypes to files')
                    sample_df = pd.DataFrame(retained_genotype, columns=header,\
                                         index=None)
                    combined_df = pd.merge(sample_df, dataframe, how='inner', \
                                           on=config.basic_header)
                    combined_df.to_csv(triodir+'/'+familyid+'.denovo.txt', \
                                       sep='\t', index=False)
                else:
                    print('Trios ', sampleids, ' have no genotypes left')
            
    elif(mode == 'family'):
        print('output the results of family and singleton')
        familydir = dir+'/families'
        if (not os.path.exists(familydir)):
            os.makedirs(familydir)
        familyid_list = ped_df[ped_df['Person']=='Proband']['FamilyID'].tolist()
        familyid_list = [id for id in familyid_list if id != 'Unknown']
        for familyid in familyid_list:
            persons = ped_df[ped_df['FamilyID']==familyid]['Person'].tolist()
            #print(persons)
            if(len(persons) >= 3):
                if( 'Proband' in persons and 'Mother' in persons \
                         and 'Father' in persons):
                    print('This is a trio in a family')
                else:
                    print('This is a family without a trio')
            elif(len(persons) == 2):               
                print('This is a family with two members')
            else:
                if('Proband' not in persons):
                    print('This is a family with single non-proband member')
    else:
        print('Error: ' + mode +' wrong mode')

def output_xlinked(output, dataframe, genotype_df, mode, pedf): 
    pass