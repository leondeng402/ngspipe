#! /usr/local/bin/python3.4
#####################################################################
#
#    config.py: the configuration file for database files and their
#               headers, and file names for the intermediate files
#    author:  Liyong Deng
#    copyright (c) 2014 by Liyong Deng
#
#####################################################################
buildver='hg19'
dir_humandb='/home/leon/biobin/annovar20141017/annovar/humandb'
basic_header = ['chr', 'start', 'end', 'ref', 'alt']
txtinput_header = ['chr', 'start', 'end', 'ref', 'alt', 'comments'] 
vcfinput_header = ['chr', 'start', 'end', 'ref', 'alt', 'zygosity', 'qual', \
                  'dep'] 
#####################################################################
#  region-based: -regionanno
#####################################################################
cytoband = 'cytoBand'
cytoband_suffix = '_cytoBand'
cytoband_header = ['label', 'cytoband', 'chr', 'start', 'end', 'ref', 'alt']
cytoband_use_header = cytoband_header[1:]
superdup = 'genomicSuperDups'
superdup_suffix = '_genomicSuperDups'
superdup_header = ['label', 'segmentaldups', 'chr', 'start', 'end', 'ref', \
                   'alt']
superdup_use_header = superdup_header[1:]
gwas = 'gwasCatalog'
gwas_suffix = '_gwasCatalog'
gwas_header = ['label', 'gwas_association', 'chr', 'start', 'end', 'ref', \
               'alt']
gwas_use_header = gwas_header[1:]
#####################################################################
#  gene level:'--geneanno -dbtype refGene'
#####################################################################
refgene = 'refGene'
refgene_exn_header = ['linenum', 'variant_function','geneinfo', 'chr', \
                      'start', 'end', 'ref', 'alt']
refgene_exn_use_header = refgene_exn_header[1:]
refgene_var_header = ['variant_function','geneinfo', 'chr', 'start', 'end', \
                      'ref', 'alt']
refgene_var_use_header =refgene_var_header
knowngene='knownGene'
knowngene_exn_header = ['linenum', 'k_variant_function','k_geneinfo', 'chr', \
                        'start', 'end', 'ref', 'alt']
knowngene_exn_use_header = knowngene_exn_header[1:]
knowngene_var_header = ['k_variant_function','k_geneinfo', 'chr', 'start', \
                        'end', 'ref', 'alt']
knowngene_var_use_header =knowngene_var_header

ensgene='ensGene'
ensgene_exn_header = ['linenum', 'e_variant_function','e_geneinfo', 'chr', \
                      'start', 'end', 'ref', 'alt']
ensgene_exn_use_header =ensgene_exn_header[1:]
ensgene_var_header = ['e_variant_function','e_geneinfo', 'chr', 'start', \
                      'end', 'ref', 'alt']
ensgene_var_use_header =ensgene_var_header


gene_var_suffix ='.variant_function'
gene_exn_suffix ='.exonic_variant_function'

####################################################################
#  variant leve: filter-based
####################################################################
snp='snp138'
snp_drp_suffix='_snp138_dropped'
snp_flt_suffix='_snp138_filtered'
snp_header= ['snpver', 'snpID', 'chr', 'start', 'end', 'ref', 'alt']
#snp_use_header = ['chr', 'start', 'end', 'ref', 'alt', 'snpID']
snp_use_header = snp_header[1:]

g1k='1000g2014sep'
g1k_afr_drp_suffix='_AFR.sites.2014_09_dropped'
g1k_afr_flt_suffix='_AFR.sites.2014_09_filtered'
g1k_afr_header = ['label','g1k_afr_allele_freq', 'chr', 'start', 'end', 'ref',\
                  'alt']
g1k_afr_use_header = g1k_afr_header[1:]

g1k_all_drp_suffix='_ALL.sites.2014_09_dropped'
g1k_all_flt_suffix='_ALL.sites.2014_09_filtered'
g1k_all_header = ['label','g1k_all_allele_freq', 'chr', 'start', 'end', 'ref',\
                  'alt']
g1k_all_use_header = g1k_all_header[1:]

g1k_amr_drp_suffix='_AMR.sites.2014_09_dropped'
g1k_amr_flt_suffix='_AMR.sites.2014_09_filtered'
g1k_amr_header = ['label','g1k_amr_allele_freq', 'chr', 'start', 'end', 'ref',\
                  'alt']
g1k_amr_use_header = g1k_amr_header[1:]

g1k_eas_drp_suffix='_EAS.sites.2014_09_dropped'
g1k_eas_flt_suffix='_EAS.sites.2014_09_filtered'
g1k_eas_header = ['label','g1k_eas_allele_freq', 'chr', 'start', 'end', 'ref',\
                  'alt']
g1k_eas_use_header = g1k_eas_header[1:]

g1k_eur_drp_suffix='_EUR.sites.2014_09_dropped'
g1k_eur_flt_suffix='_EUR.sites.2014_09_filtered'
g1k_eur_header = ['label','g1k_eur_allele_freq', 'chr', 'start', 'end', 'ref',\
                  'alt']
g1k_eur_use_header = g1k_eur_header[1:]

g1k_sas_drp_suffix='_SAS.sites.2014_09_dropped'
g1k_sas_flt_suffix='_SAS.sites.2014_09_filtered'
g1k_sas_header = ['label','g1k_sas_allele_freq', 'chr', 'start', 'end', 'ref',\
                  'alt']
g1k_sas_use_header = g1k_sas_header[1:]

esp='esp6500si'
esp_aa_drp_suffix='_esp6500si_aa_dropped'
esp_aa_flt_suffix='_esp6500si_aa_filtered'
esp_aa_header = ['label','esp_aa_allele_freq', 'chr', 'start', 'end', 'ref', \
              'alt']
esp_aa_use_header = esp_aa_header[1:]

esp_all_drp_suffix='_esp6500si_all_dropped'
esp_all_flt_suffix='_esp6500si_all_filtered'
esp_all_header = ['label','esp_all_allele_freq', 'chr', 'start', 'end', 'ref',\
              'alt']
esp_all_use_header = esp_all_header[1:]

esp_ea_drp_suffix='_esp6500si_ea_dropped'
esp_ea_flt_suffix='_esp6500si_ea_filtered'
esp_ea_header = ['label','esp_aa_allele_freq', 'chr', 'start', 'end', 'ref',\
              'alt']
esp_ea_use_header = esp_ea_header[1:]

cg='cg69'
cg_header = ['label','cg69_freq', 'chr', 'start', 'end', 'ref', 'alt']
cg_use_header = cg_header[1:]
popfreq='popfreq_all'
popfreq_header = ['label','popfreq_all', 'chr', 'start', 'end', 'ref', 'alt']
popfreq_use_header = popfreq_header[1:]

ljb='ljb23'
ljb_all_drp_suffix='_ljb23_all_dropped'
ljb_all_flt_suffix='_ljb23_all_filtered'
ljb_header = ['label','sift_score', 'sift_score_converted', 'sift_pred', 
              'pp2_hdiv_score', 'pp2_hdiv_pred',
              'pp2_hvar_score', 'pp2_hvar_pred',
              'lrt_score', 'lrt_score_converted', 'lrt_pred',
              'mt_score', 'mt_score_converted', 'mt_pred',
              'ma_score', 'ma_score_converted', 'ma_pred',
              'fathmm_score', 'fathmm_score_converted', 'fathmm_pred',
              'radialsvm_score', 'radialsvm_score_converted', 'radialsvm_pred',
              'lr_score', 'lr_pred', 'gerp++', 'phylop', 'siphy',               
              'chr', 'start', 'end', 'ref', 'alt']
ljb_use_header = ljb_header[1:]


clinvar='clianvar_20140211'
clinvar_header = ['label','clinvar', 'chr', 'start', 'end', 'ref', 'alt']
clinvar_use_header = clinvar_header[1:]
# cadd data directly from washington
cadd='caddgt10'
cadd_header = ['label','cadd_raw', 'cadd_phred', 'chr', 'start', 'end', \
               'ref', 'alt']
cadd_use_header = cadd_header[1:]
# hgmd database configuration
hgmd_file= '/home/leon/lab/hgmd/hgmd_pro.allmut.txt'
hgmd_header = ['disease', 'gene', 'chrom', 'genename', 'gdbid', 'omimid', \
              'amino', 'deletion', 'insertion', 'codon', 'codonAff', 'descr',\
              'hgvs', 'hgvsAll', 'dbsnp', 'chr', 'start', \
              'end', 'tag', 'author', 'fullname', 'allname', 'vol', \
              'page', 'year', 'pmid', 'reftag', 'comments', 'acc_num', \
              'new_date',  'base']
hgmd_use_header = ['chr', 'start', 'end', 'disease', 'gene', 'genename',\
                   'hgvsAll', 'codon']

