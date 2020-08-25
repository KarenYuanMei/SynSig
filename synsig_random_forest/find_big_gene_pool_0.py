#Goal: to find the genes in the big pool; find the genes in all of the features in features matrix
#final output file: big_pool_genes_index.csv; 11082 genes

import pandas as pd
import networkx as nx
import numpy as np

import ddot
from ddot import Ontology
import csv

import sys
import random

from collections import defaultdict

import matplotlib
matplotlib.use("TKAgg")
from matplotlib import pyplot as plt



#find the genes in all of the features:

def find_genes_in_no_brain_feature(filename):
	data=pd.read_csv("/Users/karenmei/Documents/Synapse_Ontology/normalized_no_brain_features/normalized_%s.csv"%filename, index_col=[0])
	genes=data['Norm_Symbol'].tolist()
	return genes

def find_genes_in_brain_feature(filename):
	data=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/normalized_features/normalized_%s.csv'%filename, index_col=[0])
	genes=data['Norm_Symbol'].tolist()
	return genes


#sum up all of the genes in both the no-brain-features and the brain features:
def find_all_feature_genes():
	all_genes=[]
	for item in feature_list:
		genes=find_genes_in_no_brain_feature(item)
		all_genes.append(genes)

	for item in brain_features:
		genes=find_genes_in_brain_feature(item)
		all_genes.append(genes)

	return all_genes
	

#flatten the list of gene lists into a single list:
def find_union_genes(all_genes):
	all_genes=[item for sublist in all_genes for item in sublist]
	union=list(set(all_genes))
	return union


#for each of the features, find the genes corresponding to that feature:
#output: dictionary of feature (key): gene_list (values)
def find_feature_genes():
	all_features=feature_list+brain_features
	genes_master=[]
	for item in all_features:
		print (item)
		if item in feature_list:
			genes_in_feature=find_genes_in_no_brain_feature(item)
			genes_master.append(genes_in_feature)
		if item in brain_features:
			genes_in_feature=find_genes_in_brain_feature(item)
			genes_master.append(genes_in_feature)

	tuples=list(zip(all_features, genes_master))
	feature_genes_dict=dict(tuples)
	return feature_genes_dict

#make a matrix of features with 1 or 0 indicating whether a gene is present in feature or not
def find_gene_feature_matrix(all_genes):
	features_gene_dict=find_feature_genes()
	
	union=find_union_genes(all_genes)

	counts_all_features=[]

	all_features=feature_list+brain_features
	for feature in all_features:
		counts=[]
		for gene in union:
		
			if gene in features_gene_dict[feature]:
				presence=1
			else: 
				presence=0
			counts.append(presence)
		counts_all_features.append(counts)

	df=pd.DataFrame({'Genes': union})

	for i in range(len(all_features)):
		df[all_features[i]]=counts_all_features[i]

	print (df)

	df.to_csv('all_genes_all_features.csv')
	return df



#find the pool of genes that are present in all features:
def find_big_pool_genes(df):
	#df=pd.read_csv('all_genes_all_features.csv', index_col=[0])
	df=df.set_index('Genes')

	df=df.drop(['HIP_Proteomics', 'DFC_Proteomics', 'V1C_Proteomics', 'AMY_Proteomics', 'MD_Proteomics', 'CBC_Proteomics', 'bioplex_source_feature', 'biogrid_source_feature', 'Phosphosite_mouse_no'], axis=1)

	feature_no=len(list(df.columns))

	df['sum']=df.sum(axis=1)
	
	df.loc['Total']= df.sum()

	df=df[df['sum']==feature_no]
	print (df)

	big_pool_genes=list(df.index)

	pool_genes_df=pd.DataFrame({'genes': big_pool_genes})

	pool_genes_df.to_csv('big_pool_genes_index.csv')
	return big_pool_genes

if __name__ == '__main__':

	#list of nonbrain features:

	feature_list=['colon_hpa_isoform_exp', 'ovary_hpa_isoform_exp', 'breast_hpa_isoform_exp', 'lung_hpa_isoform_exp', 'salivary gland_hpa_isoform_exp', 'seminal vesicle_hpa_isoform_exp', 'lymph node_hpa_isoform_exp', 'placenta_hpa_isoform_exp', 'kidney_hpa_isoform_exp', 'cervix, uterine_hpa_isoform_exp', 'adrenal gland_hpa_isoform_exp', 'thyroid gland_hpa_isoform_exp', 'stomach 1_hpa_isoform_exp', 'gallbladder_hpa_isoform_exp', 'duodenum_hpa_isoform_exp', 'fallopian tube_hpa_isoform_exp',
	'endometrium 1_hpa_isoform_exp', 'skin 1_hpa_isoform_exp', 'spleen_hpa_isoform_exp', 'gtex_no_brain_exp', 'appendix_hpa_isoform_exp', 'heart muscle_hpa_isoform_exp', 'small intestine_hpa_isoform_exp', 'epididymis_hpa_isoform_exp', 'testis_hpa_isoform_exp', 'liver_hpa_isoform_exp', 'esophagus_hpa_isoform_exp', 'urinary bladder_hpa_isoform_exp', 'skeletal muscle_hpa_isoform_exp', 'tonsil_hpa_isoform_exp', 'prostate_hpa_isoform_exp', 'parathyroid gland_hpa_isoform_exp',
	'adipose tissue_hpa_isoform_exp', 'smooth muscle_hpa_isoform_exp', 'rectum_hpa_isoform_exp', 'bone marrow_hpa_isoform_exp', 'mentha_source_feature', 'chr_no_source_feature', 'qPhos_site_number', 'Phosphosite_mouse_no', 'Phosphosite_hu_no', 'pFAM_domain_number', 'pFAM_domain', 'protein_mass', 'Ensembl_aa_length', 'Ensembl_isoform_no', 'trans_count', 'gc_content', 'trans_len', 'gene_length', 'exon_no', 'cds_length', 'bioplex_source_feature', 'biogrid_source_feature']

	#list of brain features:
	brain_features=['HIP_RNA', 'DFC_RNA', 'V1C_RNA', 'AMY_RNA', 'MD_RNA', 'STR_RNA', 'CBC_RNA', 'HIP_Proteomics', 'DFC_Proteomics', 'V1C_Proteomics', 'AMY_Proteomics', 'MD_Proteomics', 'CBC_Proteomics']

	#find genes in each feature
	all_genes=find_all_feature_genes()
	#make matrix of features and genes
	df=find_gene_feature_matrix(all_genes)

	df=pd.read_csv('all_genes_all_features.csv', index_col=[0])

	#find the genes that are in all features

	genes=find_big_pool_genes(df)

	print (len(genes))


