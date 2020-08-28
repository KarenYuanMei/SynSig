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

import tkinter
import matplotlib
matplotlib.use("TKAgg")
print(matplotlib.get_backend())
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
	syngo=Ontology.from_table('../features/SynGO_BP.txt')
	syngo_bp_genes=syngo.genes
	syngo=Ontology.from_table('../features/SynGO_CC.txt')
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


def plot_adult_ROC():
	
	final=find_true_y()

	#print (union_df)
	plot_ROC(final)

def plot_pr(df):
	from sklearn.datasets import make_classification
	from sklearn.linear_model import LogisticRegression
	from sklearn.model_selection import train_test_split
	from sklearn.metrics import precision_recall_curve
	from sklearn.metrics import f1_score
	from sklearn.metrics import auc
	from matplotlib import pyplot
	probs=df['avg_scores'].tolist()
	y=df['union'].tolist()

	lr_precision, lr_recall, thresholds = precision_recall_curve(y, probs)
	print (len(thresholds), len(lr_precision), len(lr_recall))
	#print (len(probs))
	print (thresholds)
	print (lr_precision)
	print (lr_recall)
	lr_auc = auc(lr_recall, lr_precision)


	#pr_curve=list(zip(thresholds, lr_precision, lr_recall))

	pr_curve=pd.DataFrame({"Precision": lr_precision, "Recall": lr_recall})

	pr_curve.to_csv('PR_on_syngo.csv')

	# summarize scores
	print('Logistic: auc=%.3f' % (lr_auc))
	# plot the precision-recall curves
	no_skill = len(df[df['union']==1]) / len(y)
	pyplot.plot([0, 1], [no_skill, no_skill], linestyle='--', label='No Skill')
	pyplot.plot(lr_recall, lr_precision, marker='.', label='Logistic')
	# axis labels
	pyplot.xlabel('Recall')
	pyplot.ylabel('Precision')
	# show the legend
	pyplot.legend()
	# show the plot
	pyplot.show()

def plot_adult_pr():
	final=find_true_y()

	#print (union_df)
	plot_pr(final)

if __name__ == '__main__':
	plot_adult_ROC()
	plot_adult_pr()
	
