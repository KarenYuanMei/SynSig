#Goal: analyze number of transcripts in the human genes as reported by Ensembl to determine if there are differences between synapse and non-synapse genes in this feature

#source: Ensembl

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

#load the source dataset
def load_data(filename):
	df=pd.read_csv(filename)
	df=df.set_index('Genes')
	#print (df)
	return df

#find the title of the data:
def load_data_col_title(filename):
	df=load_data(filename)
	cols=list(df.columns)
	return cols[0]

def load_synapse_positives_negatives():
	synapse_genes=load_synapse_positives()
	synapse_negatives=load_synapse_negatives()
	return synapse_genes, synapse_negatives


#find the data values for the positive and negative synapse genes
def find_positive_negative(filename, synapse_genes, synapse_negatives):
	#synapse_genes=load_synapse_positives()
	#synapse_negatives=load_synapse_negatives()
	data_df=load_data(filename)
	col_title=load_data_col_title(filename)

	positive_df=data_df.loc[synapse_genes]
	positive_df=positive_df.dropna()
	positives=positive_df[col_title].tolist()
	positive_values=np.array(positives)

	negative_df=data_df.loc[synapse_negatives]
	negative_df=negative_df.dropna()
	negatives=negative_df[col_title].tolist()
	negative_values=np.array(negatives)
	return positive_values, negative_values

def find_permutation(positives, negatives):
	p_value=permutation_test(positives, negatives, method='approximate', num_rounds=10000, seed=0)
	#print (p_value)
	return p_value
	#plot_distributions(positive_values, negative_values, title, col_title)
	

#for the expanded network: for the novel random forest derived network, repeat the above calculations: =========================================================================
def load_expanded_positives():
	#novel_synapse_genes=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_8/synapse_genes_above_5.3.csv', usecols=['genes'])
	novel_synapse_genes=pd.read_csv('pred_genes_above_4.7.csv', usecols=['genes'])
	novel_synapse_genes=novel_synapse_genes['genes'].tolist()
	return novel_synapse_genes

#generate an expanded list of negative synapse genes and keep it consistent across analyses:
def generate_expanded_negatives():

	negatives=pd.read_csv('negative_pool.csv')
	negatives=negatives['genes'].tolist()
	random.seed(4)
	random.shuffle(negatives)

	positives=load_expanded_positives()
	negatives=list(set(negatives)-set(positives))
	negatives.sort()
	print ('negatives', len(negatives))

	training_negative=load_synapse_negatives()
	training_positive=load_synapse_positives()
	negatives=list(set(negatives)-set(training_negative)-set(training_positive))
	negatives.sort()

	negatives=negatives[:len(positives)]
	print (negatives[:5])
	negatives_df=pd.DataFrame({'Genes': negatives})
	negatives_df.to_csv('synsig_negative_list.csv')
	return negatives_df

def load_expanded_negatives():
	#negative_genes=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Analyze_Synapse_Features/expanded_negative_list.csv')
	negative_genes=pd.read_csv('synsig_negative_list.csv')
	negatives=negative_genes['Genes'].tolist()
	return negatives

def load_expanded_positives_negatives():
	training_synapse=load_synapse_positives()
	training_negative=load_synapse_negatives()

	synapse_genes=load_expanded_positives()
	synapse_genes=list(set(synapse_genes)-set(training_synapse))

	synapse_negatives=load_expanded_negatives()
	synapse_negatives=list(set(synapse_negatives)-set(training_negative))
	return synapse_genes, synapse_negatives

def make_dataframe_R(synapse_genes, synapse_negatives, title, feature_name):
	positive_values, negative_values=find_positive_negative(synapse_genes, synapse_negatives)
	pos=pd.DataFrame({'Type': 'Positive', 'Values': positive_values})
	neg=pd.DataFrame({'Type': 'Negative', 'Values': negative_values})
	final=pd.concat([pos, neg], axis=0)
	print (final)
	#final.to_csv('%s_%s_feature_histogram_R.csv'%(title, feature_name))
	return final

#plot the distributions of the synapse vs. non-synapse proteins
def plot_training_distributions(positives, negatives, col_title):

	bins=np.histogram(np.hstack((positives,negatives)), bins=40)[1] #get the bin edges
	
	plt.hist(positives, bins, alpha=0.5, edgecolor='black', linewidth=0.5)
	plt.hist(negatives, bins, alpha=0.5, edgecolor='black', linewidth=0.5)

	#plt.ylabel('Non-Synapse Genes in Brain')
	plt.xlabel('%s'%col_title, fontweight='bold')
	plt.ylabel('Frequency', fontweight = 'bold')
	#plt.xscale('log')
	plt.yscale('log')
	plt.legend(labels=['Training Positives', 'Training Negatives'])
	#plt.show()
	#plt.close()
	plt.grid(b=False)
	#plt.savefig(title, bbox_inches='tight')
	plt.show()
	plt.close()

def plot_expanded_distributions(positives, negatives, col_title):
	bins=np.histogram(np.hstack((positives,negatives)), bins=40)[1] #get the bin edges
	
	plt.hist(positives, bins, alpha=0.5, edgecolor='black', linewidth=0.5, color='navy')
	plt.hist(negatives, bins, alpha=0.5, edgecolor='black', linewidth=0.5, color='honeydew' )

	#plt.ylabel('Non-Synapse Genes in Brain')
	plt.xlabel('%s'%col_title, fontweight='bold')
	plt.ylabel('Frequency', fontweight = 'bold')
	#plt.xscale('log')
	plt.yscale('log')
	plt.legend(labels=['SynSig Positives', 'SynSig Negatives'])
	#plt.show()
	#plt.close()
	plt.grid(b=False)
	#plt.savefig(title, bbox_inches='tight')
	plt.show()
	plt.close()


#combine the above functions to find the significance of the original synapse genes and non-synapse genes:
def test_significance(filename, synapse_genes, negative_genes):
	#positive_values, negative_values=find_positive_negative(filename, synapse_genes, synapse_negatives)
	positive_values, negative_values=find_positive_negative(filename, synapse_genes, negative_genes)
	p_value=find_permutation(positive_values, negative_values)

	pos_mean=np.mean(positive_values)
	neg_mean=np.mean(negative_values)
	norm_pos=positive_values/neg_mean
	pos_sem=stats.sem(norm_pos)
	return p_value, pos_sem, pos_mean/neg_mean

def compare_features():
	expanded_negatives=generate_expanded_negatives()
	features=['ENSEMBL_aa_length', "cds_length", "exon_no", "gc_content", "gene_length", "Ensembl_isoform_no", "trans_count", "pFAM_domain_number", "Phosphosite_hu_no", "protein_mass"]
	training_p=[]
	training_ratio=[]
	exp_p=[]
	exp_ratio=[]
	training_pos_sem_list=[]
	exp_pos_sem_list=[]
	for item in features:
		filename='/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_7/source_features/%s.csv'%item

		data_type=item
		print (item)

		col_title=load_data_col_title(filename)
		new_title = col_title.replace("_"," ") 

		#group_name='Training'
		synapse_genes, negative_genes=load_synapse_positives_negatives()
		p, pos_sem, ratio=test_significance(filename, synapse_genes, negative_genes)
		training_p.append(p)
		training_ratio.append(ratio)
		training_pos_sem_list.append(pos_sem)
		#make_dataframe_R(synapse_genes, synapse_negatives, 'training', item)
		positive_values, negative_values=find_positive_negative(filename, synapse_genes, negative_genes)
		#plot_training_distributions(positive_values, negative_values, new_title)

		#group_name='SynSig'
		synapse_genes, negative_genes=load_expanded_positives_negatives()
		p, pos_sem, ratio=test_significance(filename, synapse_genes, negative_genes)
		exp_p.append(p)
		exp_ratio.append(ratio)
		exp_pos_sem_list.append(pos_sem)
		#make_dataframe_R(positives, negatives, 'expanded', item)
		positive_values, negative_values=find_positive_negative(filename, synapse_genes, negative_genes)
		#plot_expanded_distributions(positive_values, negative_values, new_title)

	final=pd.DataFrame({'Features': features, 'Training_P': training_p, 'Training_Ratio': training_ratio, 'Training_SEM': training_pos_sem_list, 'Exp_P': exp_p, 'Exp_Ratio': exp_ratio, 'Exp_SEM': exp_pos_sem_list})
	final.to_csv('compare_train_synsig_features.csv')
	return final


def correct_fdr(df, col_name):
		col_list=np.array(df[col_name].tolist())

		corr_list=multipletests(col_list, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

		return corr_list

def corr_df(df):
	df=df.set_index('Features')
	df=df.drop(['trans_count'])
	df=df.reset_index()

	trainingP_corr=correct_fdr(df, 'Training_P')
	expP_corr=correct_fdr(df, 'Exp_P')

	features=['Amino Acid Length', 'CDS Length', 'Exon Number', 'GC Content', 'Gene Length', 'Isoform Number', 'pFAM Domain Number', 'Phosphorylation Site Number', 'Protein Mass']

	trainingFC=np.array(df['Exp_Ratio'].tolist())

	expFC=np.array(df['Training_Ratio'].tolist())

	df['Features']=features

	print (df)

	df.to_csv('compare_train_synsig_features.csv')
	return df

	
if __name__ == "__main__":
	df=compare_features()
	corr_df(df)