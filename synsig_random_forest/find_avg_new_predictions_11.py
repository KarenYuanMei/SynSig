import numpy as np
from igraph import *
import pandas as pd
import sys
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

#-----load your ontology onto HiView-------------------------------------------------------------------------
#code for uploading to HiView taken from DDOT package: https://github.com/michaelkyu/ddot/blob/master/examples/Tutorial.ipynb


#construct the random forest so that when doing the 5X cross validation, the model is not seeing 20% of the genes, not just rows--------------------
import random
import pickle

def get_file_genes(filename):
	genes=pd.read_csv(filename)
	gene_names=genes['genes'].tolist()
	return gene_names

def get_all_training(pos_file, neg_file):
	pos=get_file_genes(pos_file)
	neg=get_file_genes(neg_file)
	training=list(set(pos+neg))
	return training

def find_data_genes(training_genes):
	#new_index=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/no_brain_genes_index.csv')

	new_index=pd.read_csv('big_pool_genes_index.csv')
	all_genes=new_index['genes'].tolist()
	data_genes=list(set(all_genes)-set(training_genes))
	return data_genes


def find_avg_scores():
	pos_file='synapse_positives.csv'
	neg_file='synapse_negatives.csv'

	#make object pairs of pos_training and new gene objects:
	training=get_all_training(pos_file, neg_file)
	data_genes=find_data_genes(training)


	print ('data_genes', len(data_genes))

	pred_filename='all_gene_predictions.csv'
	pred=pd.read_csv(pred_filename, index_col=[0])
	gene1=pred['Gene1'].tolist()
	gene2=pred['Gene2'].tolist()

	avg_scores=[]
	novel_genes=[]
	for gene in data_genes:
		df1=pred.loc[pred['Gene1'] == gene]
		
		df2=pred.loc[pred['Gene2'] == gene]
		df=df1.append(df2)
		scores=df['ypredict'].tolist()
		print ('length of scores', len(scores))
		scores_np=np.array(scores)
		avg_score=np.mean(scores_np)
		
		avg_scores.append(avg_score)
		novel_genes.append(gene)

	print ('novel_genes', len(novel_genes), 'all_average_scores', len(avg_scores))

	df=pd.DataFrame({'genes': novel_genes, 'avg_scores': avg_scores})
	df.to_csv('brain_RNA_big_pool_novel_synapse_genes_avg_scores.csv')
	return df

if __name__ == '__main__':
	find_avg_scores()

