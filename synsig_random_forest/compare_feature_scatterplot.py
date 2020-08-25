import csv
import numpy as np
import math
import pandas as pd
from mlxtend.evaluate import permutation_test

import matplotlib
matplotlib.use("TKAgg")
from matplotlib import pyplot as plt

from scipy import stats

from numpy.random import seed 
from numpy.random import randn 
from scipy.stats import mannwhitneyu


import ddot
from ddot import Ontology
import random
plt.style.use('seaborn-deep')
matplotlib.rcParams.update({'font.size': 22})

from statsmodels.sandbox.stats.multicomp import multipletests
#import seaborn as sns; sns.set()

#load the source dataset
def load_data(filename):
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_7/source_features/%s.csv'%filename)
	df=df.set_index('Genes')
	#print (df)
	return df

#find the title of the data:
def load_data_col_title(filename):
	df=load_data(filename)
	cols=list(df.columns)
	return cols[0]

def load_expanded_positives():
	#novel_synapse_genes=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_8/synapse_genes_above_5.3.csv', usecols=['genes'])
	novel_synapse_genes=pd.read_csv('pred_genes_above_4.7.csv', usecols=['genes'])
	novel_synapse_genes=novel_synapse_genes['genes'].tolist()
	return novel_synapse_genes

def scale_features(v):
	mean=np.mean(v)
	std=np.std(v)
	scaled=(v-mean)/std
	return scaled

#find the data values for the positive and negative synapse genes
def find_positive_negative(filename1, filename2, synapse_genes):
	#synapse_genes=load_synapse_positives()
	#synapse_negatives=load_synapse_negatives()
	data_df=load_data(filename1)
	col_title=load_data_col_title(filename1)

	positive_df=data_df.loc[synapse_genes]
	#positive_df=positive_df.dropna()
	positives=positive_df[col_title].tolist()
	positive_values1=np.array(positives)

	data_df=load_data(filename2)
	col_title=load_data_col_title(filename2)

	positive_df=data_df.loc[synapse_genes]
	#positive_df=positive_df.dropna()
	positives=positive_df[col_title].tolist()
	positive_values2=np.array(positives)

	df=pd.DataFrame({'Genes': synapse_genes, '%s'%filename1: positive_values1, '%s'%filename2: positive_values2})
	df.to_csv('%s_%s_corr.csv'%(filename1, filename2))
	return positive_values1, positive_values2, df

def make_scatterplot(v1, v2, name1, name2):
	#v1 = v1 + 0.1 * np.random.rand(len(v1)) -0.05
	#random.seed(0)
	#v2 = v2+ 0.5 * np.random.rand(len(v2))
	plt.scatter(v1, v2, s=10, alpha=0.5, edgecolors='darkblue')
	plt.xticks(rotation=45)
	#plt.title('Feature Correlation')
	plt.yscale('log')
	plt.xscale('log')
	plt.xlabel(name1)
	plt.ylabel(name2)
	plt.xlim(left=1000)
	#plt.ylim(bottom=1)
	plt.grid(False)
	plt.annotate('CNTNAP2', color='red', xy=(v1[506], v2[506]), xytext=(v1[506]+1, v2[506]-2),
             arrowprops=dict(facecolor='darkgray', lw=1, arrowstyle='->'))

	plt.show()
	print (stats.spearmanr(v1, v2))
	print (stats.pearsonr(v1, v2))
	print ('log corr', stats.pearsonr(np.log10(v1), np.log10(v2)))

	#plt.annotate('IL1RAPL2', (v1[299], v2[299]))
	#print (v1[299], v2[299])
	

def sn_scatterplot(name1, name2, df):
	sns.stripplot(x=name1, y=name2, data=df, jitter=0.4)
	#sns.despine()
	plt.yscale('log')
	plt.xscale('log')
	plt.show()

synsig=load_expanded_positives()

v1, v2, df=find_positive_negative("gene_length", "Ensembl_isoform_no", synsig)
print (len(v1), len(v2))

make_scatterplot(v1, v2, 'Gene Length', 'Isoform Number')
#sn_scatterplot("gene_length", "Ensembl_isoform_no", df)


v1, v2, df=find_positive_negative("gene_length", "pFAM_domain_number", synsig)
print (len(v1), len(v2))

#v1=scale_features(v1)
#v2=scale_features(v2)

make_scatterplot(v1, v2, 'Gene Length', 'pFAM_domain_number')


'ENSEMBL_aa_length'

v1, v2, df=find_positive_negative("gene_length", "ENSEMBL_aa_length", synsig)
print (len(v1), len(v2))

#v1=scale_features(v1)
#v2=scale_features(v2)

make_scatterplot(v1, v2, 'Gene Length', 'ENSEMBL_aa_length')

v1, v2, df=find_positive_negative("protein_mass", "ENSEMBL_aa_length", synsig)
print (len(v1), len(v2))

#v1=scale_features(v1)
#v2=scale_features(v2)

make_scatterplot(v1, v2, 'protein_mass', 'ENSEMBL_aa_length')