
#Goal: Define train-test gene pair objects for random forest to discover novel synapse genes

#do this on the server

import numpy as np
import pandas as pd
import csv
from scipy.stats.stats import pearsonr
from collections import defaultdict
from itertools import combinations, product
from scipy import spatial

#import networkx as nx
#import os
#import ddot
#from ddot import Ontology

import pickle

from define_gene_objects_run_rf_4 import Gene, PairOfGenes, define_features,find_input_features, load_feature, create_feature_value_dict, get_feature_value,create_GO_score_dict, create_gene_list, load_positives, find_pos_genes_in_training


def get_file_genes(filename):
	genes=pd.read_csv(filename)
	gene_names=genes['genes'].tolist()
	return gene_names

def get_all_training(pos_file, neg_file):
	pos=get_file_genes(pos_file)
	neg=get_file_genes(neg_file)
	training=list(set(pos+neg))
	return training

def create_gene_pair_objects(gene_pairs):
	gene_pair_objects=[]
	for item in gene_pairs:
		gene1=item[0]
		gene2=item[1]
		pair_objects=PairOfGenes(gene1, gene2, include_GO=False)
		gene_pair_objects.append(pair_objects)
	return gene_pair_objects


def find_data_genes(training_genes):
	#new_index=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/no_brain_genes_index.csv')

	new_index=pd.read_csv('big_pool_genes_index.csv')
	all_genes=new_index['genes'].tolist()
	data_genes=list(set(all_genes)-set(training_genes))
	return data_genes



def create_data_gene_list(gene_names,is_test_gene,feature_value_dict):
	#returns a list of Gene objects, corresponding to the names in gene_names
	#feature_value_dict is a dictionary containing all feature values for all genes

	feature_list = ['colon_hpa_isoform_exp', 'ovary_hpa_isoform_exp', 'breast_hpa_isoform_exp', 'lung_hpa_isoform_exp', 'salivary gland_hpa_isoform_exp', 'seminal vesicle_hpa_isoform_exp', 
		'lymph node_hpa_isoform_exp', 'placenta_hpa_isoform_exp', 'kidney_hpa_isoform_exp', 'cervix, uterine_hpa_isoform_exp', 'adrenal gland_hpa_isoform_exp', 'thyroid gland_hpa_isoform_exp', 
		'stomach 1_hpa_isoform_exp', 'gallbladder_hpa_isoform_exp', 'duodenum_hpa_isoform_exp', 'fallopian tube_hpa_isoform_exp','endometrium 1_hpa_isoform_exp', 'skin 1_hpa_isoform_exp', 
		'spleen_hpa_isoform_exp', 'gtex_rna_tissue_expression', 'appendix_hpa_isoform_exp', 'heart muscle_hpa_isoform_exp', 'small intestine_hpa_isoform_exp', 'epididymis_hpa_isoform_exp', 'testis_hpa_isoform_exp', 
		'liver_hpa_isoform_exp', 'esophagus_hpa_isoform_exp', 'urinary bladder_hpa_isoform_exp', 'skeletal muscle_hpa_isoform_exp', 'tonsil_hpa_isoform_exp', 'prostate_hpa_isoform_exp', 
		'parathyroid gland_hpa_isoform_exp','adipose tissue_hpa_isoform_exp', 'smooth muscle_hpa_isoform_exp', 'rectum_hpa_isoform_exp', 'bone marrow_hpa_isoform_exp', 'mentha_source_feature', 
		'chr_no_source_feature', 'qPhos_site_number','Phosphosite_hu_no', 'pFAM_domain_number', 'pFAM_domain', 'protein_mass', 'Ensembl_aa_length', 'Ensembl_isoform_no', 
		'trans_count', 'gc_content', 'trans_len', 'gene_length', 'exon_no', 'cds_length']

	brain_features=['HIP_RNA', 'DFC_RNA', 'V1C_RNA', 'AMY_RNA', 'MD_RNA', 'STR_RNA', 'CBC_RNA']

	feature_list=feature_list+brain_features

	gene_list = []
	for name in gene_names:
		new_gene = Gene(name, is_test_gene)
		for feature_name in feature_list:
			feature_value = get_feature_value(name,feature_name,feature_value_dict)
			new_gene.create_feature(feature_name, feature_value)
		gene_list.append(new_gene)

	return gene_list

def create_data_sets(data_genes, training_genes):
	positives=load_positives()
	print ('one')
	positive_training_genes=find_pos_genes_in_training(training_genes, positives)
	print ('two')
	feature_value_dict=create_feature_value_dict(positive_training_genes)
	positive_training_objects=create_data_gene_list(positive_training_genes, False, feature_value_dict)
	print ('DONE')

	feature_value_dict = create_feature_value_dict(data_genes)
	data_gene_objects = create_data_gene_list(data_genes,False,feature_value_dict)
	print (len(data_gene_objects))
	print (len(positive_training_objects))
	
	train_data_pairs=product(positive_training_objects, data_gene_objects)
	return train_data_pairs

def create_data_pair_objects(gene_pairs):
	gene_pair_objects=[]
	for item in gene_pairs:
		gene1=item[0]
		gene2=item[1]
		pair_objects=PairOfGenes(gene1, gene2, include_GO=False)
		gene_pair_objects.append(pair_objects)
	return gene_pair_objects


if __name__ == '__main__':

	pos_file='synapse_positives.csv'
	neg_file='synapse_negatives.csv'


	#make object pairs of pos_training and new gene objects:
	training=get_all_training(pos_file, neg_file)

	data_genes=find_data_genes(training)
	
	train_data_pairs=create_data_sets(data_genes, training)
	print ('three')
	train_data_pair_objects=create_data_pair_objects(train_data_pairs)
	print (len(list(train_data_pair_objects)))

	print ('DONE')

	with open('train_data_pair_objects.pkl', 'wb') as output:
		pickle.dump(train_data_pair_objects, output, pickle.HIGHEST_PROTOCOL)


	