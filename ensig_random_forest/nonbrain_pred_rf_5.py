#Goal: use all 200 synapse positives and negatives as training examples for running the random forest to find the synapse scores for all genes; 


import numpy as np
import pandas as pd
import csv
from scipy.stats.stats import pearsonr
from collections import defaultdict
from itertools import combinations, combinations_with_replacement
from itertools import product
from scipy import spatial

import networkx as nx
import os
import ddot
from ddot import Ontology

import pickle

from sklearn.ensemble import RandomForestRegressor

from sklearn import metrics
from sklearn.metrics import explained_variance_score, mean_absolute_error, r2_score
from scipy.stats.stats import pearsonr, spearmanr
import pylab
from sklearn.datasets import make_regression

from define_training_objects_pred_synapse_5 import PairOfGenes


#load the pairs of training genes (200 synapse positives and 200 negatives)
#load the pairs of training-new genes (training to non-training genes)

def load_objects():

	with open('training_gene_pair_objects.pkl', 'rb') as input:
	    training_gene_pair_objects = pickle.load(input)

	with open('train_data_pair_objects.pkl', 'rb') as input:
		train_data_pair_objects = pickle.load(input)

	return training_gene_pair_objects, train_data_pair_objects


#for the training gene pairs, find the feature values and the reference GO scores:
def find_train_array(pair_objects):
	feature_list = ['colon_hpa_isoform_exp', 'ovary_hpa_isoform_exp', 'breast_hpa_isoform_exp', 'lung_hpa_isoform_exp', 'salivary gland_hpa_isoform_exp', 'seminal vesicle_hpa_isoform_exp', 
		'lymph node_hpa_isoform_exp', 'placenta_hpa_isoform_exp', 'kidney_hpa_isoform_exp', 'cervix, uterine_hpa_isoform_exp', 'adrenal gland_hpa_isoform_exp', 'thyroid gland_hpa_isoform_exp', 
		'stomach 1_hpa_isoform_exp', 'gallbladder_hpa_isoform_exp', 'duodenum_hpa_isoform_exp', 'fallopian tube_hpa_isoform_exp','endometrium 1_hpa_isoform_exp', 'skin 1_hpa_isoform_exp', 
		'spleen_hpa_isoform_exp', 'gtex_no_brain_exp', 'appendix_hpa_isoform_exp', 'heart muscle_hpa_isoform_exp', 'small intestine_hpa_isoform_exp', 'epididymis_hpa_isoform_exp', 'testis_hpa_isoform_exp', 
		'liver_hpa_isoform_exp', 'esophagus_hpa_isoform_exp', 'urinary bladder_hpa_isoform_exp', 'skeletal muscle_hpa_isoform_exp', 'tonsil_hpa_isoform_exp', 'prostate_hpa_isoform_exp', 
		'parathyroid gland_hpa_isoform_exp','adipose tissue_hpa_isoform_exp', 'smooth muscle_hpa_isoform_exp', 'rectum_hpa_isoform_exp', 'bone marrow_hpa_isoform_exp',
		'chr_no_source_feature', 'qPhos_site_number','Phosphosite_hu_no', 'pFAM_domain_number', 'pFAM_domain', 'protein_mass', 'Ensembl_aa_length', 'Ensembl_isoform_no', 
		'trans_count', 'gc_content', 'trans_len', 'gene_length', 'exon_no', 'cds_length']

	feature_array=[]
	score_array=[]
	for item in pair_objects:
		pair_GO_score=item.GO_score
		score_array.append(pair_GO_score)
		pair_feature_array=[]
		for feature_name in feature_list:
			pair_feature_values=item.__dict__[feature_name]
			pair_feature_array.append(pair_feature_values)
		feature_array.append(pair_feature_array)
	feature_array=np.array(feature_array)
	score_array=np.array(score_array)
	return feature_array, score_array

#for the training to non-training gene pairs, find the feature values, and the identity of the first gene (gene1) of the pair, and the identity of the second gene (gene2) of the pair
def find_data_array(pair_objects):
	feature_list = ['colon_hpa_isoform_exp', 'ovary_hpa_isoform_exp', 'breast_hpa_isoform_exp', 'lung_hpa_isoform_exp', 'salivary gland_hpa_isoform_exp', 'seminal vesicle_hpa_isoform_exp', 
		'lymph node_hpa_isoform_exp', 'placenta_hpa_isoform_exp', 'kidney_hpa_isoform_exp', 'cervix, uterine_hpa_isoform_exp', 'adrenal gland_hpa_isoform_exp', 'thyroid gland_hpa_isoform_exp', 
		'stomach 1_hpa_isoform_exp', 'gallbladder_hpa_isoform_exp', 'duodenum_hpa_isoform_exp', 'fallopian tube_hpa_isoform_exp','endometrium 1_hpa_isoform_exp', 'skin 1_hpa_isoform_exp', 
		'spleen_hpa_isoform_exp', 'gtex_no_brain_exp', 'appendix_hpa_isoform_exp', 'heart muscle_hpa_isoform_exp', 'small intestine_hpa_isoform_exp', 'epididymis_hpa_isoform_exp', 'testis_hpa_isoform_exp', 
		'liver_hpa_isoform_exp', 'esophagus_hpa_isoform_exp', 'urinary bladder_hpa_isoform_exp', 'skeletal muscle_hpa_isoform_exp', 'tonsil_hpa_isoform_exp', 'prostate_hpa_isoform_exp', 
		'parathyroid gland_hpa_isoform_exp','adipose tissue_hpa_isoform_exp', 'smooth muscle_hpa_isoform_exp', 'rectum_hpa_isoform_exp', 'bone marrow_hpa_isoform_exp',
		'chr_no_source_feature', 'qPhos_site_number','Phosphosite_hu_no', 'pFAM_domain_number', 'pFAM_domain', 'protein_mass', 'Ensembl_aa_length', 'Ensembl_isoform_no', 
		'trans_count', 'gc_content', 'trans_len', 'gene_length', 'exon_no', 'cds_length']

	feature_array=[]
	gene1_all=[]
	gene2_all=[]
	for item in pair_objects:
		gene1=item.gene1_name
		gene1_all.append(gene1)
		gene2=item.gene2_name
		gene2_all.append(gene2)
		pair_feature_array=[]
		for feature_name in feature_list:
			pair_feature_values=item.__dict__[feature_name]
			pair_feature_array.append(pair_feature_values)
		feature_array.append(pair_feature_array)
	feature_array=np.array(feature_array)
	return feature_array, gene1_all, gene2_all

# def find_score_array(pair_objects):
# 	score_array=[]
# 	for item in pair_objects:
# 		pair_GO_score=item.GO_score
# 		score_array.append(pair_GO_score)
# 	score_array=np.array(score_array)
# 	return score_array


#run the random forest to train on all of the training gene pairs and use it to predict the training-new gene pairs:
def run_rf():
	training_gene_pair_objects, train_data_pair_objects=load_objects()
	X_train, y_train=find_train_array(training_gene_pair_objects)
	print (X_train.shape)

	forest = RandomForestRegressor(n_estimators=100, oob_score=True, random_state=0)
	#forest = RandomForestRegressor(200)
	forest.fit(X_train, y_train)

	#------actual new genes-------------------------------------------------------------------------------------------------------------

	data_test, data_gene1, data_gene2=find_data_array(train_data_pair_objects)
	print (data_test.shape)
	print (data_test)

	#print('nan', np.isnan(data_test).any())
	nan_idx=np.where(np.isnan(data_test))
	data_test[nan_idx]=0
	#print(np.where(np.isnan(data_test)))
	#print('infinity', np.isfinite(data_test).all())

	data_fit=forest.predict(data_test)

	print (len(data_fit))

	df=pd.DataFrame({'ypredict':data_fit})

	df['Gene1']=data_gene1
	df['Gene2']=data_gene2

	df = df[['Gene1', 'Gene2', 'ypredict']]
	print (df)
	df.to_csv('nonbrain_all_gene_predictions.csv')

#---when need to run file, use the following command-------------
if __name__ == '__main__':
	run_rf()