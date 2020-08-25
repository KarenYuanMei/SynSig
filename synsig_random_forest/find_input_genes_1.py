#Goal: to select random genes for synapse and non-synapse categories:

#select positive synapse examples by using the intersection of synsysnet, syndb, and go-synapse


import pandas as pd
import networkx as nx
import numpy as np

import ddot
from ddot import Ontology
import csv

import sys
import random


def GO_Focus(Term):
	ndex_server ='http://public.ndexbio.org' 
	ndex_user, ndex_pass = 'ym2', 'Synapse'
	go_human = Ontology.from_ndex('http://public.ndexbio.org/v2/network/16565851-4bfc-11e9-9f06-0ac135e8bacf')
	go_genes=go_human.genes
	subont= go_human.focus(Term)
	subont_genes=subont.genes
	return go_genes, subont_genes

def find_synDB_genes():
	synDB=pd.read_csv('/Users/karenmei/Documents/BrainHierarchyDataSource/SynDB_Master.csv')
	synDB_genes=synDB['Symbol'].tolist()
	return synDB_genes

def find_synsysnet_genes():
	synsysnet=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Synapse_Genes/SynSysNet_genes.csv')
	synsysnet_genes=synsysnet['gene_name'].tolist()
	return synsysnet_genes

def find_pool_genes():
	union=pd.read_csv('big_pool_genes_index.csv', index_col=[0])
	union=union['genes'].tolist()
	return union

#find positive pool:
def find_positives(subont_genes, synDB_genes, synsysnet_genes, pool_genes):
	pos=list(set(subont_genes)&set(synDB_genes)&set(synsysnet_genes))
	print (len(pos))
	overlap=list(set(pos)&set(pool_genes))
	print (len(overlap))
	df=pd.DataFrame({'genes': overlap})
	df.to_csv('positive_pool.csv')
	return overlap

#find negative pool:----------------------------------------------------------------------------------------------------
def find_negatives(go_genes, subont_genes, synDB_genes, synsysnet_genes, pool_genes):
	negatives=list(set(pool_genes)-set(subont_genes)-set(synDB_genes)-set(synsysnet_genes))
	overlap=list(set(negatives)&set(pool_genes))
	print (len(overlap))
	df=pd.DataFrame({'genes': overlap})
	df.to_csv('negative_pool.csv')
	return overlap

#find positive and negative pools for input examples:
def find_pos_neg_pools():
	go_genes, subont_genes=GO_Focus('GO:0045202')
	#Extract all of the genes in the synapse in SynDB---------------------------------
	synDB_genes=find_synDB_genes()
	#Extract the genes in SynSysNet---------------------------------------------
	synsysnet_genes=find_synsysnet_genes()
	pool_genes=find_pool_genes()
	positives=find_positives(subont_genes, synDB_genes, synsysnet_genes, pool_genes)
	negatives=find_negatives(go_genes, subont_genes, synDB_genes, synsysnet_genes, pool_genes)
	return positives, negatives


#find the synapse_positives, negatives and holdout genes:
def find_pos_neg_training(gene_pool, gene_list, name):
	#gene_pool=list(df.index)

	overlap=list(set(gene_pool)&set(gene_list))
	overlap.sort()
	print (len(overlap))

	random.seed(4)

	random.shuffle(overlap)
	genes=overlap[:200]
	
	print (genes[:5])

	df=pd.DataFrame(genes, columns=['genes'])

	df.to_csv('synapse_%s.csv'%name)

	
	return genes


if __name__ == '__main__':

	#find the positive and negative pools:

	positives, negatives=find_pos_neg_pools()
	
	#find the gene pool that has information in all features
	pool_genes=find_pool_genes()

	#find the positive, negative training, test and hold-out genes

	pos=find_pos_neg_training(pool_genes, positives, 'positives')
	neg=find_pos_neg_training(pool_genes, negatives, 'negatives')
	
	print (set(pos)&set(neg))
	
