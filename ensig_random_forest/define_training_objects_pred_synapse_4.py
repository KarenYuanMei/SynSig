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
import ddot
from ddot import Ontology

import pickle


#from sklearn.ensemble import RandomForestRegressor

#from sklearn import metrics
#from sklearn.metrics import explained_variance_score, mean_absolute_error, r2_score
#from scipy.stats.stats import pearsonr, spearmanr
#import pylab
#from sklearn.datasets import make_regression


class Gene:
	def __init__(self,name,is_test_gene):
		assert isinstance(name,str), "You screwed up! The name you gave is not a string!"
		assert isinstance(is_test_gene,bool), "You screwed up! is_test_gene needs to be True or False!"

		self.name = name
		self.is_test_gene = is_test_gene

	def create_feature(self,feature_name,feature_value):
		self.__dict__[feature_name] = feature_value


	def create_GO_scores(self, GO_scores):
		#GO_scores is a dictionary with th gene2_names as key and the go_score as the value
		self.go_scores = GO_scores

class PairOfGenes:
	def __init__(self,gene1,gene2):
		assert isinstance(gene1,Gene), "You screwed up! gene1 needs to be a Gene!"
		assert isinstance(gene2,Gene), "You screwed up! gene2 needs to be a Gene!"

		self.gene1_name = gene1.name
		self.gene2_name = gene2.name

		feature_list = ['colon_hpa_isoform_exp', 'ovary_hpa_isoform_exp', 'breast_hpa_isoform_exp', 'lung_hpa_isoform_exp', 'salivary gland_hpa_isoform_exp', 'seminal vesicle_hpa_isoform_exp', 
		'lymph node_hpa_isoform_exp', 'placenta_hpa_isoform_exp', 'kidney_hpa_isoform_exp', 'cervix, uterine_hpa_isoform_exp', 'adrenal gland_hpa_isoform_exp', 'thyroid gland_hpa_isoform_exp', 
		'stomach 1_hpa_isoform_exp', 'gallbladder_hpa_isoform_exp', 'duodenum_hpa_isoform_exp', 'fallopian tube_hpa_isoform_exp','endometrium 1_hpa_isoform_exp', 'skin 1_hpa_isoform_exp', 
		'spleen_hpa_isoform_exp', 'gtex_no_brain_exp', 'appendix_hpa_isoform_exp', 'heart muscle_hpa_isoform_exp', 'small intestine_hpa_isoform_exp', 'epididymis_hpa_isoform_exp', 'testis_hpa_isoform_exp', 
		'liver_hpa_isoform_exp', 'esophagus_hpa_isoform_exp', 'urinary bladder_hpa_isoform_exp', 'skeletal muscle_hpa_isoform_exp', 'tonsil_hpa_isoform_exp', 'prostate_hpa_isoform_exp', 
		'parathyroid gland_hpa_isoform_exp','adipose tissue_hpa_isoform_exp', 'smooth muscle_hpa_isoform_exp', 'rectum_hpa_isoform_exp', 'bone marrow_hpa_isoform_exp',
		'chr_no_source_feature', 'qPhos_site_number','Phosphosite_hu_no', 'pFAM_domain_number', 'pFAM_domain', 'protein_mass', 'Ensembl_aa_length', 'Ensembl_isoform_no', 
		'trans_count', 'gc_content', 'trans_len', 'gene_length', 'exon_no', 'cds_length']
#find the genes in all of the features:

		for feature_name in feature_list:
			gene1_feature = gene1.__dict__[feature_name]
			gene2_feature = gene2.__dict__[feature_name]
			if self.check_for_missing_feature(gene1_feature, gene2_feature) == "Missing":
				self.__dict__[feature_name] = 0
			else:
				self.combine_features(gene1_feature,gene2_feature,feature_name)

		self.create_pair_GO_score(gene1, gene2)

	def check_for_missing_feature(self, gene1_feature,gene2_feature):
		#check if numpy array:
		#then check if all values are 0
		#check if not numpy array:
		#then check if values are 0
		if type(gene1_feature) is np.ndarray or type(gene2_feature) is np.ndarray:
			if np.all(gene1_feature==0) == True or np.all(gene2_feature==0)==True:
				return "Missing"
		else:
			if gene1_feature ==0 or gene2_feature == 0:
				return "Missing"
			
	def combine_features(self,gene1_feature,gene2_feature,feature_name):
	
		pearson_features = ['colon_hpa_isoform_exp', 'ovary_hpa_isoform_exp', 'breast_hpa_isoform_exp', 'lung_hpa_isoform_exp', 'salivary gland_hpa_isoform_exp', 'seminal vesicle_hpa_isoform_exp',
		'lymph node_hpa_isoform_exp', 'placenta_hpa_isoform_exp', 'kidney_hpa_isoform_exp', 'cervix, uterine_hpa_isoform_exp', 'adrenal gland_hpa_isoform_exp', 'thyroid gland_hpa_isoform_exp', 
		'stomach 1_hpa_isoform_exp', 'gallbladder_hpa_isoform_exp', 'duodenum_hpa_isoform_exp', 'fallopian tube_hpa_isoform_exp','endometrium 1_hpa_isoform_exp', 'skin 1_hpa_isoform_exp', 
		'spleen_hpa_isoform_exp', 'gtex_no_brain_exp', 'appendix_hpa_isoform_exp', 'heart muscle_hpa_isoform_exp', 'small intestine_hpa_isoform_exp', 'epididymis_hpa_isoform_exp', 'testis_hpa_isoform_exp',
		'liver_hpa_isoform_exp', 'esophagus_hpa_isoform_exp', 'urinary bladder_hpa_isoform_exp', 'skeletal muscle_hpa_isoform_exp', 'tonsil_hpa_isoform_exp', 'prostate_hpa_isoform_exp', 
		'parathyroid gland_hpa_isoform_exp','adipose tissue_hpa_isoform_exp', 'smooth muscle_hpa_isoform_exp', 'rectum_hpa_isoform_exp', 'bone marrow_hpa_isoform_exp', 
		'HIP_RNA', 'DFC_RNA', 'V1C_RNA', 'AMY_RNA', 'MD_RNA', 'STR_RNA', 'CBC_RNA']
		subtraction_features=['Phosphosite_hu_no', 'qPhos_site_number', 'Ensembl_isoform_no', 'Ensembl_aa_length', 'pFAM_domain_number', 'protein_mass', "trans_count", 'gc_content', 'trans_len', 'gene_length', 'exon_no', 'cds_length']
		jaccard_features=['pFAM_domain', 'mentha_source_feature']

		if feature_name in pearson_features:
			self.__dict__[feature_name] = pearsonr(gene1_feature,gene2_feature)[0]

		elif feature_name in subtraction_features:
			self.__dict__[feature_name]=abs(gene1_feature-gene2_feature)[0]

		elif feature_name in jaccard_features:
			no1=set(gene1_feature)
			no2=set(gene2_feature)
			self.__dict__[feature_name]=len(no1.intersection(no2)) / len(no1.union(no2))


		elif feature_name=="chr_no_source_feature":
			if gene1_feature == gene2_feature:
				result=1
				self.chr_no_source_feature=result
			else:
				result=0
				self.chr_no_source_feature=result

	def create_pair_GO_score(self, gene1, gene2,):
		gene2_name=gene2.name
		GO_score=gene1.go_scores[gene2_name]
		self.GO_score = GO_score

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

def get_file_genes(filename):
	genes=pd.read_csv(filename)
	gene_names=genes['genes'].tolist()
	return gene_names

def get_all_training(pos_file, neg_file):
	pos=get_file_genes(pos_file)
	neg=get_file_genes(neg_file)
	training=list(set(pos+neg))
	return training

def find_input_features(filename, input_genes):
	string_files=['pFAM_domain', 'mentha_source_feature','biogrid_source_feature', 'bioplex_source_feature', 'chr_no_source_feature']
	if filename not in string_files:
		#df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/normalized_no_brain_features/normalized_%s.csv'%filename)

		df=pd.read_csv('../normalized_features/normalized_%s.csv'%filename)

		symbol=df['Norm_Symbol']
		df.drop(labels=['Norm_Symbol', 'Genes'], axis=1,inplace = True)
		df.insert(0, 'Genes', symbol)
		#print (df)
		df=df.set_index('Genes')
		df=df.loc[input_genes]

	else:
		#df = pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/normalized_no_brain_features/normalized_%s.csv'%filename,converters={"Interactors": lambda x: x.strip("[]").split(", ")})
		df = pd.read_csv('../normalized_features/normalized_%s.csv'%filename,converters={"Interactors": lambda x: x.strip("[]").split(", ")})
		symbol=df['Norm_Symbol']
		df.drop(labels=['Norm_Symbol', 'Genes'], axis=1,inplace = True)
		df.insert(0, 'Genes', symbol)

		df=df.set_index('Genes')
		df=df.loc[input_genes]

	print ('DF', df)
		
	return df

def load_feature(filename, input_genes):
	feature=find_input_features(filename, input_genes)
	string_files=['pFAM_domain', 'mentha_source_feature','biogrid_source_feature', 'bioplex_source_feature', 'chr_no_source_feature']
	#feature=pd.read_csv(filename, index_col=[0])
	#fill in missing values with 0
	feature=feature.fillna(0)
	idx=list(feature.index)
	#values=feature.values
	if filename in string_files:
		#print ('TRUE')
		values=feature[feature.columns[0]]
		values=values.tolist()
		feature_dict=dict(list(zip(idx, values)))

	else:
		#in case of duplicate values, do a defaultdict
		values=feature.values
		idx_values=list(zip(idx, values))
		d=defaultdict(list)
		for idx, value in idx_values:
			d[idx].append(value)
		new_values=[]
		keys=[]
		for key in d:
			feature_values=d[key]
			keys.append(key)
			new_value=np.mean(feature_values, axis=0)
			new_values.append(new_value)
		feature_dict=dict(zip(keys, new_values))
		#print (feature_dict['KCNMA1'])
	return feature_dict



#def create_feature_value_dict():
	#returns a dictionary containing all feature values for all genes
	#format is: dict[feature_name] returns another dictionary. This dictionary d can be used as follows: d[gene_name] = feature value for gene_name

def create_feature_value_dict(input_genes):
	#source_feature_files=['Bayes_PSD', 'HPA_IHC', 'trans_count', 'gc_content', 'trans_len', 'gene_length', 'exon_no', 'cds_length','mentha_source_feature', 'bioplex_source_feature', 'biogrid_source_feature', 'CBC_Proteomics', 'STR_Proteomics', 'MD_Proteomics', 'AMY_Proteomics', 'V1C_Proteomics', 'DFC_Proteomics', 'HIP_Proteomics', 'CBC_RNA', 'STR_RNA', 'MD_RNA', 'AMY_RNA', 'V1C_RNA', 'DFC_RNA', 'HIP_RNA', 'mathieson_halflife', 'fornasiero_halflife', 'fornasiero_aa', 'hpa_isoform_exp', 'chr_no_source_feature']
	feature_list = ['colon_hpa_isoform_exp', 'ovary_hpa_isoform_exp', 'breast_hpa_isoform_exp', 'lung_hpa_isoform_exp', 'salivary gland_hpa_isoform_exp', 'seminal vesicle_hpa_isoform_exp', 
		'lymph node_hpa_isoform_exp', 'placenta_hpa_isoform_exp', 'kidney_hpa_isoform_exp', 'cervix, uterine_hpa_isoform_exp', 'adrenal gland_hpa_isoform_exp', 'thyroid gland_hpa_isoform_exp', 
		'stomach 1_hpa_isoform_exp', 'gallbladder_hpa_isoform_exp', 'duodenum_hpa_isoform_exp', 'fallopian tube_hpa_isoform_exp','endometrium 1_hpa_isoform_exp', 'skin 1_hpa_isoform_exp', 
		'spleen_hpa_isoform_exp', 'gtex_no_brain_exp', 'appendix_hpa_isoform_exp', 'heart muscle_hpa_isoform_exp', 'small intestine_hpa_isoform_exp', 'epididymis_hpa_isoform_exp', 'testis_hpa_isoform_exp', 
		'liver_hpa_isoform_exp', 'esophagus_hpa_isoform_exp', 'urinary bladder_hpa_isoform_exp', 'skeletal muscle_hpa_isoform_exp', 'tonsil_hpa_isoform_exp', 'prostate_hpa_isoform_exp', 
		'parathyroid gland_hpa_isoform_exp','adipose tissue_hpa_isoform_exp', 'smooth muscle_hpa_isoform_exp', 'rectum_hpa_isoform_exp', 'bone marrow_hpa_isoform_exp',
		'chr_no_source_feature', 'qPhos_site_number','Phosphosite_hu_no', 'pFAM_domain_number', 'pFAM_domain', 'protein_mass', 'Ensembl_aa_length', 'Ensembl_isoform_no', 
		'trans_count', 'gc_content', 'trans_len', 'gene_length', 'exon_no', 'cds_length']
	source_feature_files=feature_list
#find the genes in all of the features:

	all_feature_values=[]
	items=[]
	for item in source_feature_files:
		feature_values=load_feature(item, input_genes)
		items.append(item)
		all_feature_values.append(feature_values)

	feature_dict=dict(zip(items, all_feature_values))
	return feature_dict

def get_feature_value(gene_name, feature_name, feature_value_dict):
	#returns the feature value for a specific gene
	feature_values_for_all_genes = feature_value_dict[feature_name]
	feature_value = feature_values_for_all_genes[gene_name]
	return feature_value
	
def create_GO_score_dict():
	#df=pd.read_csv('GO_training_score_matrix_for_big_pool_genes.csv', index_col=[0])
	df=pd.read_csv('../brain_RNA_big_gene_pool_pipeline/GO_training_score_matrix_for_big_pool_genes.csv', index_col=[0])
	
	
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
	positives=pd.read_csv('../brain_RNA_big_gene_pool_pipeline/synapse_positives.csv')
	
	positives=positives['genes'].tolist()
	return positives

def find_pos_genes_in_training(training_genes, positives):
	input_genes=[]
	for item in training_genes:
		if item in positives:
			input_genes.append(item)
	input_genes=list(set(input_genes))
	#print (len(input_genes))
	return input_genes

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
	pos_file='../brain_RNA_big_gene_pool_pipeline/synapse_positives.csv'
	neg_file='../brain_RNA_big_gene_pool_pipeline/synapse_negatives.csv'

	#create objects pairs of training genes (200 synapse pos and 200 synapse neg):
	training_pairs=create_training_sets(pos_file, neg_file)

	training_gene_pair_objects=create_gene_pair_objects(training_pairs)
	print ('training', len(list(training_gene_pair_objects)))

	with open('training_gene_pair_objects.pkl', 'wb') as output:
		pickle.dump(training_gene_pair_objects, output, pickle.HIGHEST_PROTOCOL)

	print ('DONE')

	