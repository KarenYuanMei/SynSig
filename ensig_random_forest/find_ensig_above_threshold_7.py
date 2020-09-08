#Goal: 
#.     1) Find the synapse genes above a threshold

import numpy as np
import pandas as pd
import csv

import ddot
from ddot import Ontology

import matplotlib
#matplotlib.use("TKAgg")
from matplotlib import pyplot as plt

from sklearn.metrics import auc

from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
plt.rcParams.update({'font.size': 22})
plt.rcParams["font.family"] = "Arial"

#Load the files with the avg predicted scores for each gene:
def load_predicted_df():
	df=pd.read_csv('nonbrain_RNA_big_pool_novel_synapse_genes_avg_scores.csv', index_col=[0])
	#print ('pred', df)
	return df


def find_synapse_above_threshold(df, threshold):
	#print ('df', df)
	df=df[df['avg_scores']>threshold]
	#print ('thresholded', df)

	novel_genes=df['genes'].tolist()
	df=pd.DataFrame({'genes': novel_genes})

	df.to_csv('nonbrain_pred_genes_above_%s.csv'%threshold)
	return novel_genes



def find_network_above_threshold(threshold):
	pred=load_predicted_df()

	novel_genes=find_synapse_above_threshold(pred, threshold)
	print ('novel', len(novel_genes))

	#pred_val_final=pred_val.set_index('Genes')
	#print (novel_genes)

	#thres_sum=pred_val_final.loc[novel_genes]

	#print (thres_sum[thres_sum['sum']>0])
	return novel_genes

if __name__ == '__main__':
	#find_frequency_val_df()
	#plot_lit_ROC()
	find_network_above_threshold(4.67)

