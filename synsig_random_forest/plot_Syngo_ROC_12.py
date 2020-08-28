#Goal: to plot ROC of predicted gene scores with SynGO; generate Supp Fig 2a

import pandas as pd
#import networkx as nx
import numpy as np
#import os
import ddot
from ddot import Ontology
import csv

import sys
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

#sys.path.append("/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_8/")
#from scipy.stats import hypergeom


import matplotlib
#matplotlib.use("TKAgg")
#print(matplotlib.get_backend())
from matplotlib import pyplot as plt

#from matplotlib_venn import venn3, venn3_circles
#from matplotlib_venn import venn2, venn2_circles
#import venn

from sklearn.metrics import auc

from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
plt.rcParams.update({'font.size': 22})
plt.rcParams["font.family"] = "Helvetica"

def load_syngo_genes():
	syngo=Ontology.from_table('../prev_databases/SynGO_BP.txt')
	syngo_bp_genes=syngo.genes
	syngo=Ontology.from_table('../prev_databases/SynGO_CC.txt')
	syngo_cc_genes=syngo.genes
	syngo_genes=list(set(syngo_bp_genes+syngo_cc_genes))
	return syngo_genes
	
def load_training_genes():
	filename='synapse_positives.csv'
	df=pd.read_csv(filename, index_col=[0])
	pos=df['genes'].tolist()

	filename='synapse_negatives.csv'
	df=pd.read_csv(filename, index_col=[0])
	neg=df['genes'].tolist()

	genes=list(set(pos+neg))
	return genes

#Load the files with the avg predicted scores for each gene:
def load_predicted_df():
	df=pd.read_csv('brain_RNA_big_pool_novel_synapse_genes_avg_scores.csv', index_col=[0])
	#print ('pred', df)
	return df

def find_true_y():
	full=load_syngo_genes()
	training_genes=load_training_genes()
	full=list(set(full)-set(training_genes))

	df=load_predicted_df()
	avg_scores=df['avg_scores'].tolist()
	pred_genes=df['genes'].tolist()

	y_list=[]
	for item in pred_genes:
		if item in full:
			group=1
		else:
			group=0
		y_list.append(group)

	final=pd.DataFrame({'genes': pred_genes, 'avg_scores': avg_scores , 'union': y_list})
	return final


def plot_ROC(df):
	probs=df['avg_scores'].tolist()
	y=df['union'].tolist()

	fpr, tpr, thresholds = roc_curve(y, probs)

	true_false=list(zip(tpr, fpr))

	ROC=list(zip(thresholds, tpr, fpr))

	ROC_df=pd.DataFrame({'Threshold': thresholds, "True_Positives": tpr, "False_Positives": fpr})

	ROC_df.to_csv('ROC_on_syngo.csv')

	#print (ROC)
	auc = roc_auc_score(y, probs)
	print('AUC: %.3f' % auc)

	idx_list=ROC_df[ROC_df['Threshold'] > 4.7].index.tolist()
	rnd_idx=idx_list[-1]
	x_axis=fpr[rnd_idx]
	y_axis=tpr[rnd_idx]
	print (x_axis, y_axis)

	#import matplotlib.pyplot as plt
	#plt.title('ROC Curve for Novel Synapse Predictions Using Experiments')
	f=plt.figure()
	
	plt.plot(fpr, tpr, 'g', label = 'AUC = %0.2f' %auc, linewidth=3)
	plt.legend(loc = 'lower right')
	plt.plot([0, 1], [0, 1],'--', color='black')
	plt.plot([0, x_axis], [y_axis, y_axis],'--', color='darkgray')
	plt.plot([x_axis, x_axis], [0, y_axis],'--', color='darkgray')
	
	#rnd_idx=517
	plt.annotate('SynSig\nThreshold', color='black', xy=(fpr[rnd_idx], tpr[rnd_idx]), xytext=(fpr[rnd_idx]+0.05, tpr[rnd_idx]),
             arrowprops=dict(facecolor='darkgray', lw=2, arrowstyle='->'),)
	plt.xlim([0, 1])
	plt.ylim([0, 1])
	plt.ylabel('True Positive Rate')
	plt.xlabel('False Positive Rate')
	plt.show()
	f.savefig("ROC_Syngo_plot.pdf", bbox_inches='tight')


def plot_adult_ROC():
	
	final=find_true_y()

	#print (union_df)
	plot_ROC(final)


if __name__ == '__main__':
	plot_adult_ROC()

	
