#Goal: to generate gene objects and gene_pair objects for all genes so that we can identify synapse genes from non-brain features. 

#do this on the server

import numpy as np
import pandas as pd
import csv
from scipy.stats.stats import pearsonr
from collections import defaultdict
from itertools import combinations, product
from scipy import spatial

import networkx as nx
import os
import ddot
from ddot import Ontology

import pickle

from define_ensig_objects_run_rf_0 import define_features,Gene,PairOfGenes,find_input_genes,find_input_features,load_feature, create_feature_value_dict, get_feature_value, find_pos_genes_in_training


def get_file_genes(filename):
	genes=pd.read_csv('../synsig_random_forest/%s'%filename)
	gene_names=genes['genes'].tolist()
	return gene_names

def get_all_training(pos_file, neg_file):
	pos=get_file_genes(pos_file)
	neg=get_file_genes(neg_file)
	training=list(set(pos+neg))
	return training

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

def find_input_genes(training_genes, test_genes):
	input_genes=list(set(training_genes+test_genes))
	return input_genes
	
def create_GO_score_dict():
	#df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/GO_scores/GO_training_score_matrix_no_brain_features.csv', index_col=[0])
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
	feature_list=define_features()

	gene_list = []
	for name in gene_names:
		new_gene = Gene(name, is_test_gene)
		for feature_name in feature_list:
			feature_value = get_feature_value(name,feature_name,feature_value_dict)
			new_gene.create_feature(feature_name, feature_value)
		gene_list.append(new_gene)

	return gene_list

def load_positives():
	#positives=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_9/input_genes/synapse_positives.csv')
	positives=pd.read_csv('../synsig_random_forest/synapse_positives.csv')
	
	positives=positives['genes'].tolist()
	return positives


def create_training_and_test_sets(training_file, test_file):
	training_gene_names = get_training_gene_names(training_file)
	test_gene_names = get_test_gene_names(test_file)
	input_genes=find_input_genes(training_gene_names, test_gene_names)
	print (input_genes)
	feature_value_dict = create_feature_value_dict(input_genes)
	training_gene_objects = create_gene_list(training_gene_names,False,feature_value_dict)
	print ('number of training objects', len(training_gene_objects))
	training_pairs=combinations(training_gene_objects,2)
	test_gene_objects=create_gene_list(test_gene_names, True, feature_value_dict)
	print ('number of test objects', len(test_gene_objects))
	positives=load_positives()
	positive_training_genes=find_pos_genes_in_training(training_gene_names, positives)
	print ('positive_training_genes', len(positive_training_genes))
	positive_training_objects=create_gene_list(positive_training_genes, False, feature_value_dict)
	train_test_pairs=product(positive_training_objects, test_gene_objects)

	#test_pairs=combinations(test_gene_objects, 2)
	return training_pairs, train_test_pairs

def create_gene_pair_objects(gene_pairs):
	gene_pair_objects=[]
	for item in gene_pairs:
		gene1=item[0]
		gene2=item[1]
		pair_objects=PairOfGenes(gene1, gene2, include_GO=False)
		gene_pair_objects.append(pair_objects)
	return gene_pair_objects


def load_positives_negatives():
	#positives=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_9/input_genes/synapse_positives.csv')
	positive_file='../synsig_random_forest/synapse_positives.csv'
	
	#positives=pd.read_csv('synapse_positives.csv')
	positives=positives['genes'].tolist()
	#negatives=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/synapse_negatives.csv')
	positive_file='../synsig_random_forest/synapse_negatives.csv'
	
	#negatives=pd.read_csv('synapse_negatives.csv')
	negatives=negatives['genes'].tolist()
	return positives, negatives

def find_data_genes(training_genes):
	#new_index=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/no_brain_genes_index.csv')

	new_index=pd.read_csv('../synsig_random_forest/big_pool_genes_index.csv')
	all_genes=new_index['genes'].tolist()
	data_genes=list(set(all_genes)-set(training_genes))
	return data_genes



def create_data_gene_list(gene_names,is_test_gene,feature_value_dict):
	#returns a list of Gene objects, corresponding to the names in gene_names
	#feature_value_dict is a dictionary containing all feature values for all genes

	feature_list=define_features()

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

def find_pos_genes_in_training(training_genes, positives):
	input_genes=[]
	for item in training_genes:
		if item in positives:
			input_genes.append(item)
	input_genes=list(set(input_genes))
	#print (len(input_genes))
	return input_genes


if __name__ == '__main__':
	pos_file='../synsig_random_forest/synapse_positives.csv'
	neg_file='../synsig_random_forest/synapse_negatives.csv'

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



	

	