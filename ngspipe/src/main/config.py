#! /usr/local/bin/python3.4
#####################################################################
#
#    config.py: the configuration file for database files and their
#               headers, and file names for the intermediate files
#    author:  Liyong Deng
#    copyright (c) 2014 by Liyong Deng in Dr. Chung's Lab
#
#####################################################################
buildver='hg19'
dir_humandb='/home/leon/biobin/annovar20141017/annovar/humandb'
basic_header = ['Chr', 'Start', 'End', 'Ref', 'Alt']
txtinput_header = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Comments'] 
vcfinput_header = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Zygosity', 'Qual', \
                  'Dep']
vcfoutput_header = basic_header + ['A', 'B', 'C'] 
#####################################################################
#  region-based: -regionanno
#####################################################################
cytoband = 'cytoBand'
superdup = 'genomicSuperDups'
gwas = 'gwasCatalog'

#####################################################################
#  gene level:'--geneanno -dbtype refGene'
#####################################################################
refgene = 'refGene'
knowngene='knownGene'
ensgene='ensGene'
refgene_header = ['VariantFunction', 'Gene', 'GeneDetail', 'ExonicFunction', \
                  'AAChange']
knowngene_header = ['VariantFunction.k', 'Gene.k', 'GeneDetail.k', \
                    'ExonicFunction.k', 'AAChange.k']
ensgene_header = ['VariantFunction.e', 'Gene.e', 'GeneDetail.e', \
                  'ExonicFunction.e', 'AAChange.e']

rvis_file='/home/leon/lab/rvis/RVIS_20141023.csv'
rvis_header=['Gene', 'ALL_0.01percent', 'ALL_0.1percent', \
             'ALL_1percent', 'PP2_aLL_0.1percent', 'Ea_0.1percent', \
             'Ea_1percent', 'Aa_0.1percent', 'Aa_1percent', 'OEratio']


####################################################################
#  variant leve: filter-based
####################################################################
snp='snp138'
#esp='esp6500si_all'
cg='cg69'
popfreq='popfreq_all'
ljb='ljb23_all'
clinvar='clinvar_20140929'
cosmic='cosmic70'
nci='nci60'

# g1k allele frequency
g1k='1000g2014sep'
g1k_drp_suffix='.sites.2014_09_dropped'
g1k_flt_suffix='.sites.2014_09_filtered'
g1k_races=['all', 'afr', 'amr', 'eas', 'eur', 'sas']
g1k_header=['Label', 'Value', 'Chr', 'Start', 'End', 'Ref', 'Alt']
g1k_use_header = [g1k.upper() + '_' + g1k_races[0].upper(), \
                  g1k.upper() + '_' + g1k_races[1].upper(), \
                  g1k.upper() + '_' + g1k_races[2].upper(), \
                  g1k.upper() + '_' + g1k_races[3].upper(), \
                  g1k.upper() + '_' + g1k_races[4].upper(), \
                  g1k.upper() + '_' + g1k_races[5].upper()]

# cadd raw score and phred score
cadd='caddgt10'
cadd_header = ['Label','CADD_raw', 'CADD_phred', 'Chr', 'Start', 'End', \
               'Ref', 'Alt']
cadd_use_header = cadd_header[1:]
cadd_drp_suffix='_'+cadd+'_dropped'
cadd_flt_suffix='_'+cadd+'_filtered'
# hgmd database configuration
hgmd_file= '/home/leon/lab/hgmd/hgmd_pro.allmut.txt'
hgmd_header = ['HGMD_Disease', 'HGMD_Gene', 'Chrom', 'HGMD_Genename', 'gdbid', 
               'Omimid', 'Amino', 'Deletion', 'Insertion', 'HGMD_Codon', 
               'CodonAff', 'Descr','Hgvs', 'HGMD_hgvsAll', 'dbsnp', 'Chr', 'Start', 
              'End', 'Tag', 'Author', 'HGMD_Magzine', 'Allname', 'Vol', 
              'Page', 'HGMD_PublishYear', 'HGMD_PMID', 'Reftag', 'Comments', 'ACC_NUM', 
              'New_date',  'Base']
hgmd_use_header = ['Chr', 'Start', 'End', 'HGMD_Disease', 'HGMD_Gene', \
                   'HGMD_Genename', 'HGMD_hgvsAll', 'HGMD_Codon', 'HGMD_PMID',\
                   'HGMD_Magzine', 'HGMD_PublishYear']
#####################################################################
#  table_annovar: protocols, operations
#####################################################################                     
protocols=[refgene, knowngene, ensgene, cytoband, superdup, gwas, snp,\
            cg, popfreq, ljb, clinvar, cosmic, nci]
operations=['g', 'g', 'g', 'r', 'r', 'r','f', 'f', 'f', 'f', 'f', 'f', \
            'f']
nastring='.'

# table_annovar contains a bug to include 1000G2012APR results
table_anno_header=['Chr', 'Start', 'End', 'Ref', 'Alt', \
                        'VariantFunction', 'Gene', 'GeneDetail', \
                        'ExonicFunction', 'AAChange', 'VariantFunction.k',\
                        'Gene.k', 'GeneDetail.k', 'ExonicFunction.k', \
                        'AAChange.k', 'VariantFunction.e', 'Gene.e', \
                        'GeneDetail.e',    'ExonicFunction.e', 'AAChange.e',\
                        'CytoBand', 'GenomicSuperDups', 'GwasCatalog', \
                        'dbsnp138',  'CG69', 'PopFreqMax', \
                        '1000G2012APR_ALL', '1000G2012APR_AFR', \
                        '1000G2012APR_AMR', '1000G2012APR_ASN', \
                        '1000G2012APR_EUR', 'ESP6500si_ALL', 'ESP6500si_AA',\
                        'ESP6500si_EA', 'CG46', 'LJB23_SIFT_score', \
                        'LJB23_SIFT_score_converted', 'LJB23_SIFT_pred', \
                        'LJB23_Polyphen2_HDIV_score', \
                        'LJB23_Polyphen2_HDIV_pred', \
                        'LJB23_Polyphen2_HVAR_score', \
                        'LJB23_Polyphen2_HVAR_pred', 'LJB23_LRT_score',\
                        'LJB23_LRT_score_converted', 'LJB23_LRT_pred',  \
                        'LJB23_MutationTaster_score', \
                        'LJB23_MutationTaster_score_converted', \
                        'LJB23_MutationTaster_pred', \
                        'LJB23_MutationAssessor_score', \
                        'LJB23_MutationAssessor_score_converted', \
                        'LJB23_MutationAssessor_pred', \
                        'LJB23_FATHMM_score', 'LJB23_FATHMM_score_converted', \
                        'LJB23_FATHMM_pred', 'LJB23_RadialSVM_score', \
                        'LJB23_RadialSVM_score_converted', \
                        'LJB23_RadialSVM_pred', 'LJB23_LR_score', \
                        'LJB23_LR_pred', 'LJB23_GERP++', 'LJB23_PhyloP', \
                        'LJB23_SiPhy', 'clinvar_20140929', 'cosmic70', 'nci60']

table_anno_suffix = '_multianno.csv'
# final display header
gene_header = ['VariantFunction', 'Gene', 'GeneDetail', 'ExonicFunction', \
               'AAChange', 'VariantFunction.k', 'Gene.k', 'GeneDetail.k', \
               'ExonicFunction.k', 'AAChange.k', 'VariantFunction.e', \
               'Gene.e', 'GeneDetail.e',    'ExonicFunction.e',    'AAChange.e']
region_header = ['CytoBand', 'GenomicSuperDups', 'GwasCatalog']

ljb_header = ['LJB23_SIFT_score', 'LJB23_SIFT_score_converted', \
              'LJB23_SIFT_pred', \
              'LJB23_Polyphen2_HDIV_score', 'LJB23_Polyphen2_HDIV_pred', \
              'LJB23_Polyphen2_HVAR_score', 'LJB23_Polyphen2_HVAR_pred', \
              'LJB23_LRT_score', 'LJB23_LRT_score_converted', \
              'LJB23_LRT_pred', \
              'LJB23_MutationTaster_score', \
              'LJB23_MutationTaster_score_converted',  \
              'LJB23_MutationTaster_pred', \
              'LJB23_MutationAssessor_score', \
              'LJB23_MutationAssessor_score_converted', \
              'LJB23_MutationAssessor_pred', \
              'LJB23_FATHMM_score', 'LJB23_FATHMM_score_converted', \
              'LJB23_FATHMM_pred', \
              'LJB23_RadialSVM_score', \
              'LJB23_RadialSVM_score_converted', \
              'LJB23_RadialSVM_pred', \
              'LJB23_LR_score', 'LJB23_LR_pred', \
              'LJB23_GERP++', 'LJB23_PhyloP', 'LJB23_SiPhy']

ljb_pred_header=['LJB23_SIFT_pred', 'LJB23_Polyphen2_HDIV_pred', \
                 'LJB23_Polyphen2_HVAR_pred', 'LJB23_LRT_pred', \
                 'LJB23_MutationTaster_pred', 'LJB23_MutationAssessor_pred', \
                 'LJB23_FATHMM_pred', 'LJB23_RadialSVM_pred', 'LJB23_LR_pred',\
                 'LJB23_GERP++', 'LJB23_PhyloP', 'LJB23_SiPhy']
allele_frequency_header = ['dbsnp138'] + g1k_use_header \
                           + ['ESP6500si_ALL', 'ESP6500si_AA', 'ESP6500si_EA'] \
                           +  ['CG69', 'PopFreqMax']
other_filter_header = ['clinvar_20140929', 'cosmic70', 'nci60']
custom_header = basic_header + gene_header + region_header \
                + allele_frequency_header + ljb_pred_header + cadd_header[1:3] \
                + other_filter_header + hgmd_use_header[3:] 
merged_custom_header = basic_header + refgene_header + region_header \
                + allele_frequency_header + ljb_pred_header + cadd_header[1:3] \
                + other_filter_header + hgmd_use_header[3:] 
                              
final_header = table_anno_header[:26] + table_anno_header[31:34] \
               + table_anno_header[35:] \
               + g1k_use_header + cadd_header[1:3] + hgmd_use_header[3:]
