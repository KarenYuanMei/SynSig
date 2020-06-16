#Goal: analyze tissue expression in the synapse; #analyze the difference between synapse positive vs. negative genes in terms of this feature:
#source: HPA:
#RNA GTEx tissue gene data: https://www.proteinatlas.org/about/download; number 6
#Transcript expression levels summarized per gene in 36 tissues based on RNA-seq. The tab-separated file includes Ensembl gene identifier ("Gene"), analysed sample ("Tissue"), transcripts per million ("TPM"), protein-transcripts per million ("pTPM") and normalized expression ("NX"). The data was obtained from GTEx and is based on The Human Protein Atlas version 19 and Ensembl version 92.38.
import csv
import math
import numpy as np

from scipy.stats.stats import pearsonr  

from mlxtend.evaluate import permutation_test

import pandas as pd

import mygene

from collections import Counter
from collections import defaultdict

import matplotlib
matplotlib.use("TKAgg")
from matplotlib import pyplot as plt

#import seaborn as sns; sns.set()
from scipy import stats

from numpy.random import seed 
#from numpy.random import randn 

import ddot
from ddot import Ontology

import random
from statsmodels.sandbox.stats.multicomp import multipletests


#plt.style.use('seaborn-bright')

plt.style.use('seaborn-deep')
#plt.style.use('seaborn-deep')
matplotlib.rcParams.update({'font.size': 22})


#source datafile for this feature:
filename='/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_7/source_features/gtex_rna_tissue_expression.csv'

#load the positive synapse genes (training)

def load_synapse_positives():
	#synapse_genes=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_8/input_genes/synapse_positives.csv')
	synapse_genes=pd.read_csv('synapse_positives.csv')
	synapse_genes=synapse_genes['genes'].tolist()
	return synapse_genes

#load the negative synapse genes (training)
def load_synapse_negatives():
	#synapse_negatives=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_8/input_genes/synapse_negatives.csv')
	synapse_genes=pd.read_csv('synapse_negatives.csv')
	synapse_negatives=synapse_genes['genes'].tolist()
	return synapse_negatives

#for the expanded network: for the novel random forest derived network, repeat the above calculations: =========================================================================
def load_expanded_positives():
	#novel_synapse_genes=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_8/synapse_genes_above_5.3.csv', usecols=['genes'])
	novel_synapse_genes=pd.read_csv('pred_genes_above_4.7.csv', usecols=['genes'])
	novel_synapse_genes=novel_synapse_genes['genes'].tolist()
	return novel_synapse_genes

def load_expanded_negatives():
	#negative_genes=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Analyze_Synapse_Features/expanded_negative_list.csv')
	negative_genes=pd.read_csv('synsig_negative_list.csv')
	negatives=negative_genes['Genes'].tolist()
	return negatives

def load_expanded_positives_negatives():
	training_synapse=load_synapse_positives()
	training_negative=load_synapse_negatives()
	training_genes=list(set(training_synapse+training_negative))

	synapse_genes=load_expanded_positives()
	synapse_genes=list(set(synapse_genes)-set(training_genes))

	synapse_negatives=load_expanded_negatives()
	synapse_negatives=list(set(synapse_negatives)-set(training_genes))
	return synapse_genes, synapse_negatives

#load the source dataset
def load_tissue_expression(tissue):
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_7/source_features/gtex_rna_tissue_expression.csv', usecols=['Genes', tissue])
	df=df.set_index('Genes')
	print (df)
	return df

#find the list of all tissues in the data
def find_tissue_list():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_7/source_features/gtex_rna_tissue_expression.csv')
	df=df.set_index('Genes')
	cols=list(df.columns)
	return cols

#find all the brain tissues
def find_brain_tissue_list():
	brain_tissues=['midbrain', 'pituitary gland', 'hippocampal formation', 'cerebellum', 'hypothalamus', 'basal ganglia', 'cerebral cortex', 'amygdala']
	return brain_tissues

#find the data values for the positive and negative synapse genes
def find_positive_negative(tissue):
	tissue_df=load_tissue_expression(tissue)
	synapse_genes=load_synapse_positives()
	synapse_negatives=load_synapse_negatives()
	positive_df=tissue_df.loc[synapse_genes]
	positive_df=positive_df.dropna()
	positives=positive_df[tissue].tolist()
	positive_values=np.array(positives)

	negative_df=tissue_df.loc[synapse_negatives]
	negative_df=negative_df.dropna()
	negatives=negative_df[tissue].tolist()
	negative_values=np.array(negatives)
	return positive_values, negative_values


#find the data values for the positive and negative synapse genes
def find_expanded_positive_negative(tissue, expanded_positives, expanded_negatives):
	tissue_df=load_tissue_expression(tissue)
	synapse_genes=expanded_positives
	synapse_negatives=expanded_negatives
	positive_df=tissue_df.loc[synapse_genes]
	positive_df=positive_df.dropna()
	positives=positive_df[tissue].tolist()
	positive_values=np.array(positives)

	negative_df=tissue_df.loc[synapse_negatives]
	negative_df=negative_df.dropna()
	negatives=negative_df[tissue].tolist()
	negative_values=np.array(negatives)
	return positive_values, negative_values

#plot the distributions of the synapse vs. non-synapse proteins
def plot_distributions(positives, negatives, tissue):

	bins=np.histogram(np.hstack((positives,negatives)), bins=40)[1] #get the bin edges
	
	plt.hist(positives, bins, alpha=0.5, edgecolor='black', linewidth=0.5, color='darkred')
	plt.hist(negatives, bins, alpha=0.5, edgecolor='black', linewidth=0.5, color='lightblue')

	plt.xlabel('mRNA Expression (TPM)', fontweight='bold')
	plt.ylabel('Frequency', fontweight = 'bold')
	#plt.xscale('log')
	plt.yscale('log')
	
	#plt.legend(labels=['Training Synapse', 'Training Negatives'])
	plt.legend(labels=['SynSig Genes', 'Non-SynSig Genes'])
	#plt.grid(b=None)
	plt.grid(False)
	plt.show()
	plt.close()


def analyze_by_tissue(tissue_list, positive_genes, negative_genes):
	fc_list=[]
	sem_list=[]
	p_list=[]
	for item in tissue_list:
		print (item)
		col_list=['Genes', item]
		data=pd.read_csv(filename, usecols=col_list)
		data=data.set_index('Genes')
		#print (data)
		data['Mean']=data.mean(axis=1)

		positive_df=data.loc[positive_genes]
		#print ('positive', positive_df)
		positive_df=positive_df.dropna()
		positives=positive_df['Mean'].tolist()
		#print (positives)
		positive_values=np.array(positives)
		#print ('positives', len(positive_genes))
		#print ('positive', positive_df)
		#print (len(positive_values))

		negative_df=data.loc[negative_genes]
		#print ('negative', len(negative_genes), negative_genes[:5])
		#print ('negative', negative_df)
		negative_df=negative_df.dropna()


		negatives=negative_df['Mean'].tolist()
		#print (len(negatives))
		negative_values=np.array(negatives)
		#print (negative_values[:5])
		
		fc=np.mean(positive_values)/np.mean(negative_values)

		norm_pos=positive_values/np.mean(negative_values)
		pos_sem=stats.sem(norm_pos)

		p_value=permutation_test(positive_values, negative_values, method='approximate', num_rounds=10000, seed=0)

		#print (len(positive_values), len(negative_values))

		#plot_distributions(positive_values, negative_values, item)

		fc_list.append(fc)
		p_list.append(p_value)
		sem_list.append(pos_sem)


	tissue_list=[x.capitalize() for x in tissue_list]

	df=pd.DataFrame({'Tissues': tissue_list, 'Fold_Change': fc_list, 'Significance': p_list, 'SEM': sem_list})
	return df


def analyze_tissues_by_group(tissue_list, positives, negatives, group):
	df=analyze_by_tissue(tissue_list, positives, negatives)
	final_df = df.rename({'Fold_Change': '%s_FC'%group, 'Significance': '%s_Significance'%group}, axis=1)
	return final_df 

def correct_fdr(df, col_name):
		col_list=np.array(df[col_name].tolist())
		corr_list=multipletests(col_list, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
		return corr_list

def analyze_training_tissues(brain_tissue_list, non_brain_tissues):
	positive_genes=load_synapse_positives()
	negative_genes=load_synapse_negatives()

	brain_training=analyze_tissues_by_group(brain_tissue_list, positive_genes, negative_genes, 'Training')
	corr_brain_train=correct_fdr(brain_training, 'Training_Significance')
	
	brain_training['Corr_Significance']=corr_brain_train[1]
	print (brain_training)
	brain_training.to_csv('brain_training.csv')

	nonbrain_training=analyze_tissues_by_group(non_brain_tissues, positive_genes, negative_genes, 'Training')
	corr_nonbrain_train=correct_fdr(nonbrain_training, 'Training_Significance')
	
	nonbrain_training['Corr_Significance']=corr_nonbrain_train[1]
	print (nonbrain_training)
	#print (nonbrain_df)
	nonbrain_training.to_csv('nonbrain_training.csv')
	return brain_training, nonbrain_training

def analyze_expanded_tissues(brain_tissue_list, non_brain_tissues):
	#positive_genes, negative_genes=load_expanded_positives_negatives()

	positive_genes=load_expanded_positives()
	negative_genes=load_expanded_negatives()
	print (len(negative_genes), negative_genes[:5])

	brain_synsig=analyze_tissues_by_group(brain_tissue_list, positive_genes, negative_genes, 'SynSig')
	corr_brain_synsig=correct_fdr(brain_synsig, 'SynSig_Significance')
	
	brain_synsig['Corr_Significance']=corr_brain_synsig[1]
	print (brain_synsig)
	brain_synsig.to_csv('brain_synsig.csv')

	nonbrain_synsig=analyze_tissues_by_group(non_brain_tissues, positive_genes, negative_genes, 'SynSig')
	corr_nonbrain_synsig=correct_fdr(nonbrain_synsig, 'SynSig_Significance')
	
	nonbrain_synsig['Corr_Significance']=corr_nonbrain_synsig[1]
	print (nonbrain_synsig)
	nonbrain_synsig.to_csv('nonbrain_synsig.csv')
	#print (nonbrain_df)
	return brain_synsig, nonbrain_synsig

		
if __name__ == "__main__":
	#analyze_by_tissue(brain_tissue_list, 'All_Brain_Tissues')
	brain_tissue_list=find_brain_tissue_list()
	print (brain_tissue_list)

	tissue_list=find_tissue_list()

	non_brain_tissues=list(set(tissue_list)-set(brain_tissue_list))

	#brain_training,  nonbrain_training=analyze_training_tissues(brain_tissue_list, non_brain_tissues)

	brain_synsig, nonbrain_synsig=analyze_expanded_tissues(brain_tissue_list, non_brain_tissues)

	