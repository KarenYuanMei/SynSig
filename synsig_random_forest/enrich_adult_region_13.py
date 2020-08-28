#Goal: #1: analyze the new experimental results and compare the enrichment of SynGO and Synch
		#2: find the previously unvalidated fetal specific genes that are now validated by experiments
		#3: find the previously unvalidated adult specific genes that are now validated by experiments
		#4: draw the intersection between fetal and adult and predicted genes for figure
import pandas as pd
#import networkx as nx
import numpy as np
#import os
import ddot
from ddot import Ontology
#from collections import defaultdict
#import collections
import csv

#import sys
#import os
#os.environ['KMP_DUPLICATE_LIB_OK']='True'

#sys.path.append("/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_8/")
from scipy.stats import hypergeom

import matplotlib
#matplotlib.use("TKAgg")
#print(matplotlib.get_backend())
from matplotlib import pyplot as plt

from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles
import venn

def load_training_genes():
	filename='synapse_positives.csv'
	df=pd.read_csv(filename, index_col=[0])
	pos=df['genes'].tolist()

	filename='synapse_negatives.csv'
	df=pd.read_csv(filename, index_col=[0])
	neg=df['genes'].tolist()

	genes=list(set(pos+neg))
	return genes

def load_pred_genes():
	pred_file='pred_genes_above_4.7.csv'
	pred_df=pd.read_csv(pred_file, index_col=[0])
	pred_genes=pred_df['genes'].tolist()
	pred_genes=[x.upper() for x in pred_genes]
	return pred_genes

def find_pred_no_training():
	pred_genes=load_pred_genes()
	training_genes=load_training_genes()
	pred_no_training=list(set(pred_genes)-set(training_genes))
	return pred_no_training


def load_syngo_genes():
	syngo=Ontology.from_table('../prev_databases/SynGO_BP.txt')
	syngo_bp_genes=syngo.genes
	syngo=Ontology.from_table('../prev_databases/SynGO_CC.txt')
	syngo_cc_genes=syngo.genes
	syngo_genes=list(set(syngo_bp_genes+syngo_cc_genes))
	training=load_training_genes()
	syngo_genes=list(set(syngo_genes)-set(training))
	return syngo_genes


def find_hypergeometric(genes, pred_no_training):

	overlap=list(set(genes)&set(pred_no_training))
	M=10683
	#M=20000
	N=len(genes)
	n=len(pred_no_training)
	x=len(overlap)
	#print ('x', x)

	#print ('M', M, 'N', N, 'n', n, 'x', x)
	#x=190
	pval = hypergeom.sf(x-1, M, n, N)

	rv = hypergeom(M, n, N)
	distr = np.arange(0, n+1)
	#print (x)
	prob = rv.pmf(distr)

	maximum=np.max(prob)
	result = np.where(prob == maximum)
	#print (result)
	#result=result.tolist()
	result=result[0]
	#print (result)
	fold=x/result
	fold=fold.tolist()
	print ('Fold Enrichment', fold)
	print ('hypergeometric p-value', pval)
	return fold


def load_data_genes_no_training(data_file):
	genes=load_data_genes(data_file)
	training=load_training_genes()
	final=list(set(genes)-set(training))
	return final

def load_adult_ctx():
	df=pd.read_csv('../experimental_validation/weijun_ctx_uniprot.csv', sep='\t')
	#print (df)
	genes=df['To'].tolist()
	training=load_training_genes()
	genes=[x.upper() for x in genes]
	genes=list(set(genes)-set(training))
	return genes

def load_adult_str():
	df=pd.read_csv('../experimental_validation/weijun_str_uniprot.csv', sep='\t')
	#print (df)
	genes=df['To'].tolist()
	training=load_training_genes()
	genes=[x.upper() for x in genes]
	genes=list(set(genes)-set(training))
	return genes

def find_synsysnet():
	df=pd.read_csv('../prev_databases/SynSysNet_genes.csv')
	#print (df)
	genes=df['gene_name'].tolist()
	#print (len(genes))
	training=load_training_genes()
	genes=list(set(genes)-set(training))
	return genes

def find_synDB():
	df=pd.read_csv('../prev_databases/SynDB_Master.csv')
	#print (df)
	genes=df['Symbol'].tolist()
	training=load_training_genes()
	genes=list(set(genes)-set(training))
	return genes

def find_GO_synapse():
	df=pd.read_csv('../prev_databases/GO_Synapse.csv')
	#print (df)
	genes=df['genes'].tolist()
	training=load_training_genes()
	genes=list(set(genes)-set(training))
	return genes

if __name__ == '__main__':
	
	pred=load_pred_genes()

	print ('cortex')
	ctx=load_adult_ctx()
	find_hypergeometric(ctx, pred)

	print ('striatum')
	stria=load_adult_str()
	find_hypergeometric(stria, pred)

	print ('overlap')
	overlap=list(set(ctx)&set(stria))
	find_hypergeometric(overlap, pred)


	print ('syngo')
	syngo=load_syngo_genes()
	find_hypergeometric(pred, syngo)
	#find_hypergeometric(stria, syngo)






