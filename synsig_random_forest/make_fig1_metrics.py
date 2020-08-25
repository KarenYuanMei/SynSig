#Goal: to plot the ROC curve with python scikit learn; plots mean with standard deviation; fig 2c

import numpy as np
from igraph import *
import pandas as pd
import sys
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

#-----load your ontology onto HiView-------------------------------------------------------------------------
#code for uploading to HiView taken from DDOT package: https://github.com/michaelkyu/ddot/blob/master/examples/Tutorial.ipynb


import networkx as nx

import matplotlib
matplotlib.use("TKAgg")
from matplotlib import pyplot as plt

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
matplotlib.rcParams.update({'font.size': 22})



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

	#data=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/random_forest/ypredict_ytest_%s.csv'%i, index_col=[0])

	#data=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/random_forest/no_ppi_ypredict_ytest_%s.csv'%i, index_col=[0])

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
	#table.to_csv('table.csv')

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

#this is the MASTER function that combines all of the above functions to plot an average ROC curve:
def plot_mean_ROC():

	#tprs = []
	#mean_fpr = np.linspace(0, 1, 80)

	#auc_list=[]

	for i in range(5):
		final=make_avg_score_df(i)
		probs=final['mean'].tolist()
		y=final['group'].tolist()

		fpr, tpr, thresholds = roc_curve(y, probs)

		#true_false=list(zip(tpr, fpr))

		#ROC=list(zip(thresholds, tpr, fpr))

		#ROC_df=pd.DataFrame({'Threshold': thresholds, "True_Positives": tpr, "False_Positives": fpr})

		#ROC_df.to_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/random_forest/ROC_df_%s.csv'%i)

		#print (ROC)
		auc = roc_auc_score(y, probs)
		#auc_list.append(auc)
		print('AUC: %.3f' % auc)


		ns_probs = [0 for _ in range(len(y))]
		ns_auc = roc_auc_score(y, ns_probs)
		ns_fpr, ns_tpr, _ = roc_curve(y, ns_probs)

		lr_auc = roc_auc_score(y, probs)
		fpr, tpr, _ = roc_curve(y, probs)
		#tprs.append(np.interp(mean_fpr, fpr, tpr))
		#tprs[-1][0] = 0.0

		#plt.plot(fpr, tpr, marker='.', label='Random Forest Classifier')
		#auc_mean=np.mean(auc_list)
		#print (auc_list)
		#print ('mean', np.mean(auc_list))
		plt.plot(ns_fpr, ns_tpr, linestyle='--', color='k')
		#mean_tpr = np.mean(tprs, axis=0)
		#mean_tpr[-1] = 1.0
		#mean_auc = auc(mean_fpr, mean_tpr)
		#std_auc = np.std(aucs)
		plt.plot(fpr, tpr, linewidth=3, color='navy')

		# std_tpr = np.std(tprs, axis=0)
		# tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
		# tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
		# plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
		#                  label=r'$\pm$ 1 std. dev.')
		# #plt.title('Avg ROC Curve for Predicting Synapse Genes \n 5-Fold Cross-Validation', fontweight = 'bold')

		plt.xlabel('False Positive Rate', fontweight='bold')
		plt.ylabel('True Positive Rate', fontweight='bold')
		plt.grid(False)
		# show the legend
		#plt.legend()
			# show the plot
	plt.show()

if __name__ == '__main__':
	plot_mean_ROC()