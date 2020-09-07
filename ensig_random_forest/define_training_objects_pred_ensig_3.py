import numpy as np
import pandas as pd
import csv
from scipy.stats.stats import pearsonr
from collections import defaultdict
from itertools import combinations, combinations_with_replacement
from itertools import product
from scipy import spatial

#import networkx as nx
#import os
#import ddot
#from ddot import Ontology

import pickle


#from sklearn.ensemble import RandomForestRegressor

#from sklearn import metrics
#from sklearn.metrics import explained_variance_score, mean_absolute_error, r2_score
#from scipy.stats.stats import pearsonr, spearmanr
#import pylab
#from sklearn.datasets import make_regression
from define_ensig_objects_run_rf_0 import define_features,Gene,PairOfGenes,find_input_genes,find_input_features,load_feature, create_feature_value_dict, get_feature_value, find_pos_genes_in_training


def get_training_gene_names(training_file):
	#function that returns a list of training gene names (strings)
	genes=pd.read_csv(training_file)
	gene_names=genes['genes'].tolist()
	return gene_names

def get_test_gene_names(test_file):
	#function that returns a list of training gene names (strings)
	genes=pd.read_csv(test_file)
	gene_names=genes['genes'].tolist()
	return gene_names

def get_file_genes(filename):
	genes=pd.read_csv(filename)
	gene_names=genes['genes'].tolist()
	return gene_names

def get_all_training(pos_file, neg_file):
	pos=get_file_genes(pos_file)
	neg=get_file_genes(neg_file)
	training=list(set(pos+neg))
	return training

	
def create_GO_score_dict():
	#df=pd.read_csv('GO_training_score_matrix_for_big_pool_genes.csv', index_col=[0])
	df=pd.read_csv('../synsig_random_forest/GO_training_score_matrix_for_big_pool_genes.csv', index_col=[0])
	
	
	idx=list(df.index)
	cols=list(df.columns)

	all_dict=[]
	for gene1 in idx:
		gene1_scores=[]
		for gene2 in cols:
			score=df.loc[gene1, gene2]
			gene1_scores.append(score)
		gene2_dict=dict(zip(cols, gene1_scores))
		all_dict.append(gene2_dict)
	
	master_dict=dict(zip(idx, all_dict))
	#print (master_dict['STX4'])
	return master_dict



def create_gene_list(gene_names,is_test_gene,feature_value_dict):
	#returns a list of Gene objects, corresponding to the names in gene_names
	#feature_value_dict is a dictionary containing all feature values for all genes

	#feature_list = ['Bayes_PSD', 'HPA_IHC', 'trans_count', 'gc_content', 'trans_len', 'gene_length', 'exon_no', 'cds_length','mentha_source_feature', 'bioplex_source_feature', 'biogrid_source_feature', 'CBC_Proteomics', 'STR_Proteomics', 'MD_Proteomics', 'AMY_Proteomics', 'V1C_Proteomics', 'DFC_Proteomics', 'HIP_Proteomics', 'CBC_RNA', 'STR_RNA', 'MD_RNA', 'AMY_RNA', 'V1C_RNA', 'DFC_RNA', 'HIP_RNA', 'mathieson_halflife', 'fornasiero_halflife', 'fornasiero_aa', 'hpa_isoform_exp', 'chr_no_source_feature']
	feature_list = ['colon_hpa_isoform_exp', 'ovary_hpa_isoform_exp', 'breast_hpa_isoform_exp', 'lung_hpa_isoform_exp', 'salivary gland_hpa_isoform_exp', 'seminal vesicle_hpa_isoform_exp', 
		'lymph node_hpa_isoform_exp', 'placenta_hpa_isoform_exp', 'kidney_hpa_isoform_exp', 'cervix, uterine_hpa_isoform_exp', 'adrenal gland_hpa_isoform_exp', 'thyroid gland_hpa_isoform_exp', 
		'stomach 1_hpa_isoform_exp', 'gallbladder_hpa_isoform_exp', 'duodenum_hpa_isoform_exp', 'fallopian tube_hpa_isoform_exp','endometrium 1_hpa_isoform_exp', 'skin 1_hpa_isoform_exp', 
		'spleen_hpa_isoform_exp', 'gtex_no_brain_exp', 'appendix_hpa_isoform_exp', 'heart muscle_hpa_isoform_exp', 'small intestine_hpa_isoform_exp', 'epididymis_hpa_isoform_exp', 'testis_hpa_isoform_exp', 
		'liver_hpa_isoform_exp', 'esophagus_hpa_isoform_exp', 'urinary bladder_hpa_isoform_exp', 'skeletal muscle_hpa_isoform_exp', 'tonsil_hpa_isoform_exp', 'prostate_hpa_isoform_exp', 
		'parathyroid gland_hpa_isoform_exp','adipose tissue_hpa_isoform_exp', 'smooth muscle_hpa_isoform_exp', 'rectum_hpa_isoform_exp', 'bone marrow_hpa_isoform_exp',
		'chr_no_source_feature', 'qPhos_site_number','Phosphosite_hu_no', 'pFAM_domain_number', 'pFAM_domain', 'protein_mass', 'Ensembl_aa_length', 'Ensembl_isoform_no', 
		'trans_count', 'gc_content', 'trans_len', 'gene_length', 'exon_no', 'cds_length']

	gene_list = []
	for name in gene_names:
		new_gene = Gene(name, is_test_gene)
		for feature_name in feature_list:
			feature_value = get_feature_value(name,feature_name,feature_value_dict)
			new_gene.create_feature(feature_name, feature_value)
		gene_list.append(new_gene)

	GO_score_dict=create_GO_score_dict()

	for gene1 in gene_list:
		gene1_name =gene1.name
		go_scores=GO_score_dict[gene1_name]
		gene1.create_GO_scores(go_scores)

	return gene_list

def load_positives():
	#positives=pd.read_csv('synapse_positives.csv')
	positives=pd.read_csv('../synsig_random_forest/synapse_positives.csv')
	
	positives=positives['genes'].tolist()
	return positives


def create_training_sets(pos_file, neg_file):
	training_gene_names=get_all_training(pos_file, neg_file)

	feature_value_dict = create_feature_value_dict(training_gene_names)
	training_gene_objects = create_gene_list(training_gene_names,False,feature_value_dict)
	print ('number of training objects', len(training_gene_objects))
	training_pairs=combinations(training_gene_objects,2)
	return training_pairs


def create_gene_pair_objects(gene_pairs):
	gene_pair_objects=[]
	for item in gene_pairs:
		gene1=item[0]
		gene2=item[1]
		pair_objects=PairOfGenes(gene1, gene2)
		gene_pair_objects.append(pair_objects)
	return gene_pair_objects

if __name__ == "__main__":
	pos_file='../synsig_random_forest/synapse_positives.csv'
	neg_file='../synsig_random_forest/synapse_negatives.csv'

	#create objects pairs of training genes (200 synapse pos and 200 synapse neg):
	training_pairs=create_training_sets(pos_file, neg_file)

	training_gene_pair_objects=create_gene_pair_objects(training_pairs)
	print ('training', len(list(training_gene_pair_objects)))

	with open('training_gene_pair_objects.pkl', 'wb') as output:
		pickle.dump(training_gene_pair_objects, output, pickle.HIGHEST_PROTOCOL)

	print ('DONE')

	