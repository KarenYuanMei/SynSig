#Goal: to find the genes that are expressed in the brain;
#final output file: brain_genes_index.csv

import pandas as pd
import networkx as nx
import numpy as np
import os
import ddot
from ddot import Ontology
import csv

import sys
import random

from collections import defaultdict

import matplotlib
matplotlib.use("TKAgg")
from matplotlib import pyplot as plt


def find_genes_in_feature(filename):
	data=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/normalized_features/normalized_%s.csv'%filename, index_col=[0])
	genes=data['Norm_Symbol'].tolist()
	return genes


def find_genes_in_feature_group(feature_list):
	all_genes=[]
	for item in feature_list:
		genes=find_genes_in_feature(item)
		all_genes.append(genes)

	all_genes=[item for sublist in all_genes for item in sublist]
	union=list(set(all_genes))
	return union

def find_all_index():
	data=pd.read_csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/big_pool_genes_index.csv')
	genes=data['genes'].tolist()
	return genes



if __name__ == '__main__':
	# nb_features=['colon_hpa_isoform_exp', 'ovary_hpa_isoform_exp', 'breast_hpa_isoform_exp', 'lung_hpa_isoform_exp', 'salivary gland_hpa_isoform_exp', 'seminal vesicle_hpa_isoform_exp', 'lymph node_hpa_isoform_exp', 'placenta_hpa_isoform_exp', 'kidney_hpa_isoform_exp', 'cervix, uterine_hpa_isoform_exp', 'adrenal gland_hpa_isoform_exp', 'thyroid gland_hpa_isoform_exp', 'stomach 1_hpa_isoform_exp', 'gallbladder_hpa_isoform_exp', 'duodenum_hpa_isoform_exp', 'fallopian tube_hpa_isoform_exp',
	# 'endometrium 1_hpa_isoform_exp', 'skin 1_hpa_isoform_exp', 'spleen_hpa_isoform_exp', 'gtex_no_brain_exp', 'appendix_hpa_isoform_exp', 'heart muscle_hpa_isoform_exp', 'small intestine_hpa_isoform_exp', 'epididymis_hpa_isoform_exp', 'testis_hpa_isoform_exp', 'liver_hpa_isoform_exp', 'esophagus_hpa_isoform_exp', 'urinary bladder_hpa_isoform_exp', 'skeletal muscle_hpa_isoform_exp', 'tonsil_hpa_isoform_exp', 'prostate_hpa_isoform_exp', 'parathyroid gland_hpa_isoform_exp',
	# 'adipose tissue_hpa_isoform_exp', 'smooth muscle_hpa_isoform_exp', 'rectum_hpa_isoform_exp', 'bone marrow_hpa_isoform_exp', 'mentha_source_feature', 'chr_no_source_feature', 'qPhos_site_number', 'Phosphosite_mouse_no', 'Phosphosite_hu_no', 'pFAM_domain_number', 'pFAM_domain', 'protein_mass', 'Ensembl_aa_length', 'Ensembl_isoform_no', 'trans_count', 'gc_content', 'trans_len', 'gene_length', 'exon_no', 'cds_length', 'bioplex_source_feature', 'biogrid_source_feature']


	#brain_features=['HIP_RNA', 'DFC_RNA', 'V1C_RNA', 'AMY_RNA', 'MD_RNA', 'STR_RNA', 'CBC_RNA']
	brain_features=['HIP_Proteomics', 'DFC_Proteomics', 'V1C_Proteomics', 'AMY_Proteomics', 'MD_Proteomics', 'CBC_Proteomics']

	brain_genes=find_genes_in_feature_group(brain_features)
	print (len(brain_genes))

	all_feature_genes=find_all_index()

	brain_genes_in_index=list(set(all_feature_genes)&set(brain_genes))

	print (len(brain_genes_in_index))
	
	brain_genes_df=pd.DataFrame({'genes': brain_genes_in_index})

	brain_genes_df.to_csv('brain_genes_index.csv')

	

	

