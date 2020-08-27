import numpy as np
import pandas as pd
import csv
from scipy.stats.stats import pearsonr
from collections import defaultdict
from itertools import combinations, combinations_with_replacement
from itertools import product
from scipy import spatial

import os
import ddot
from ddot import Ontology

import pickle

from sklearn.ensemble import RandomForestRegressor

from sklearn import metrics
from sklearn.metrics import explained_variance_score, mean_absolute_error, r2_score
from scipy.stats.stats import pearsonr, spearmanr
import pylab
#from sklearn.datasets import make_regression

from define_training_synapse_objects_8 import PairOfGenes


def load_objects():

	with open('training_gene_pair_objects.pkl', 'rb') as input:
	    training_gene_pair_objects = pickle.load(input)
	#     for item in training_gene_pair_objects:
	#     	print (item.GO_score)

	with open('train_data_pair_objects.pkl', 'rb') as input:
		train_data_pair_objects = pickle.load(input)

	return training_gene_pair_objects, train_data_pair_objects

#print (len(list(train_data_pair_objects)))
#     for item in test_gene_pair_objects:
#     	print (item.GO_score)

def find_train_array(pair_objects):
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

def find_data_array(pair_objects):
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

def find_score_array(pair_objects):
	score_array=[]
	for item in pair_objects:
		pair_GO_score=item.GO_score
		score_array.append(pair_GO_score)
	score_array=np.array(score_array)
	return score_array

def run_rf():
	training_gene_pair_objects, train_data_pair_objects=load_objects()
	X_train, y_train=find_train_array(training_gene_pair_objects)
	print (X_train.shape)

	#X_test, y_test=find_train_array(train_test_gene_pair_objects)
	#print (X_test.shape)

	forest = RandomForestRegressor(n_estimators=100, oob_score=True)
	#forest = RandomForestRegressor(200)
	forest.fit(X_train, y_train)

	#yfit=forest.predict(X_test)
	#print (pearsonr(y_test, yfit))

	#----add the training gene predictions------------------------------
	#data_test, data_gene1, data_gene2=find_data_array(training_gene_pair_objects)
	#data_true=find_score_array(training_gene_pair_objects)
	#nan_idx=np.where(np.isnan(data_test))
	#data_test[nan_idx]=0
	#data_fit=forest.predict(data_test)

	#df=pd.DataFrame({'ytrue':data_true, 'ypredict':data_fit})

	#df['Gene1']=data_gene1
	#df['Gene2']=data_gene2

	#df = df[['Gene1', 'Gene2', 'ytrue', 'ypredict']]
	#print (df)
	#df.to_csv('training_gene_pair_predictions.csv')

	#------actual new genes-------------------------------------------------------------------------------------------------------------

	data_test, data_gene1, data_gene2=find_data_array(train_data_pair_objects)
	print (data_test.shape)
	print (data_test)

	print('nan', np.isnan(data_test).any())
	nan_idx=np.where(np.isnan(data_test))
	data_test[nan_idx]=0
	print(np.where(np.isnan(data_test)))
	print('infinity', np.isfinite(data_test).all())

	data_fit=forest.predict(data_test)

	print (len(data_fit))

	df=pd.DataFrame({'ypredict':data_fit})

	df['Gene1']=data_gene1
	df['Gene2']=data_gene2

	df = df[['Gene1', 'Gene2', 'ypredict']]
	print (df)
	df.to_csv('all_gene_predictions.csv')
#---when need to run file, use the following command-------------
run_rf()