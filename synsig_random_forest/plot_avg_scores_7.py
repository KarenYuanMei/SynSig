#Goal: to plot the similarity score distributions of synapse-synapse and synapse-nonsynapse genes; plots mean with standard deviation

import numpy as np
#from igraph import *
import pandas as pd
import sys
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

from mlxtend.evaluate import permutation_test

#-----load your ontology onto HiView-------------------------------------------------------------------------
#code for uploading to HiView taken from DDOT package: https://github.com/michaelkyu/ddot/blob/master/examples/Tutorial.ipynb

import matplotlib
matplotlib.use("TKAgg")
from matplotlib import pyplot as plt
matplotlib.rcParams.update({'font.size': 22})


#import seaborn as sns; sns.set()

import pylab

#construct the random forest so that when doing the 5X cross validation, the model is not seeing 20% of the genes, not just rows--------------------
import random
import pickle
#from define_gene_objects_rf_5 import PairOfGenes

from sklearn.metrics import auc
from itertools import combinations
from collections import defaultdict

from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score



def find_test_genes(test_genes_file):
	test_genes=pd.read_csv(test_genes_file)
	test_genes=test_genes['genes'].tolist()
	print (len(test_genes))
	test_genes=list(set(test_genes))
	print ('should be 80', len(test_genes))
	return test_genes

def find_training_genes(training_genes_file):
	training_genes=pd.read_csv(training_genes_file)
	training_genes=training_genes['genes'].tolist()
	training_genes=list(set(training_genes))
	print ('should be 320', len(training_genes))
	return training_genes

def find_pos_genes_in_training(training_genes, positives):
	input_genes=[]
	for item in training_genes:
		if item in positives:
			input_genes.append(item)
	input_genes=list(set(input_genes))
	#print (len(input_genes))
	return input_genes


#convert each predicted score df into a different df with test genes as index, and all of the positive training genes as columns, so that it's easy to compute the average score per test gene to all 160 positive training genes
def make_avg_score_df(i):

	data=pd.read_csv('ypredict_ytest_%s.csv'%i)


	test_genes_file='test_genes_%s.csv'%i
	training_genes_file='training_genes_%s.csv'%i

	positive_filename='synapse_positives.csv'

	positive_genes=pd.read_csv(positive_filename)
	positives=positive_genes['genes'].tolist()

	negative_filename='synapse_negatives.csv'

	negative_genes=pd.read_csv(negative_filename)
	negatives=negative_genes['genes'].tolist()

	test_genes=find_test_genes(test_genes_file)
	training_genes=find_training_genes(training_genes_file)

	training_positives=find_pos_genes_in_training(training_genes, positives)

	print ('training_pos', len(training_positives))

	gene1=data['Gene1'].tolist()

	gene2=data['Gene2'].tolist()

	df=data[['Gene2', 'Gene1', 'ytest', 'ypredict']]

	table = pd.pivot_table(df, values='ypredict', index=['Gene2'], columns=['Gene1'], aggfunc=np.sum)
	print (table)
	table.to_csv('table.csv')

	#print (len(list(set(table.index)&set(overlap))))
	table['mean']=table.mean(axis=1)

	pos_overlap=list(set(positives)&set(test_genes))
	print ('pos_overlap', len(pos_overlap))
	print (pos_overlap[:5])


	pos_pos_df=table.loc[pos_overlap]
	print (pos_pos_df)

	pos_pos_df['group']=1

	neg_overlap=list(set(negatives)&set(test_genes))
	print ('neg_overlap', len(neg_overlap))
	print (neg_overlap[:5])

	neg_neg_df=table.loc[neg_overlap]

	neg_neg_df['group']=0

	final=pd.concat([pos_pos_df, neg_neg_df], axis=0)

	print (final)
	#final.to_csv('final.csv')

	return final

def find_distributions(gene_list1, gene_list2):
	all_scores=[]
	for gene1 in gene_list1:
		#print (gene1)
		gene2_array=master_dict[gene1]
		scores=[]
		for gene2 in gene2_array:
			if gene2 in gene_list2:
				score=gene2_array[gene2]
				scores.append(score)
		all_scores.append(scores)

	score_np=np.array(all_scores)
	score_array=np.mean(score_np, axis=1)
	#print (score_array)
	return score_array


def plot_distributions(positive_array, mixed_array):

	bins=np.histogram(np.hstack((positive_array,mixed_array)), bins=40)[1] #get the bin edges
	
	
	plt.hist(positive_array, bins, color='darkblue', label='Synapse to Synapse', alpha=0.5, edgecolor='black', linewidth=0.8)
	plt.hist(mixed_array, bins, color='grey', label='Negative to Synapse', alpha=0.5, edgecolor='black', linewidth=0.8)


	plt.xlabel('Average Predicted Semantic Similarity Per Gene', fontweight='bold')
	plt.ylabel('Frequency', fontweight = 'bold')
	#plt.title('Average Predicted Semantic Similarity', fontweight = 'bold')
	#plt.xlim(2,6)
	#plt.ylim(0,8)
	plt.grid(False)
	plt.legend()
	plt.show()
	plt.close()

def plot_sim_score_dist():
	all_pos=[]
	all_neg=[]
	for i in range(5):
		final=make_avg_score_df(i)
		print (final)
		pos=final[final['group']==1]
		print (pos)
		pos_mean=pos['mean'].tolist()
		all_pos.append(pos_mean)
		

		neg=final[final['group']==0]
		print (neg)
		neg_mean=neg['mean'].tolist()
		all_neg.append(neg_mean)

	pos_list = [item for sublist in all_pos for item in sublist]
	neg_list = [item for sublist in all_neg for item in sublist]

	p_value = permutation_test(pos_list, neg_list,
	                           method='approximate',
	                           num_rounds=10000,
	                           seed=0)
	print(p_value)

	plot_distributions(pos_list, neg_list)

if __name__ == '__main__':
	plot_sim_score_dist()
