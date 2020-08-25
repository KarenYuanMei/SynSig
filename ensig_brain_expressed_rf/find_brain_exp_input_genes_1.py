#Goal: to select random genes for synapse and non-synapse categories:

#select positive synapse examples by using the intersection of synsysnet, syndb, and go-synapse


import pandas as pd
import networkx as nx
import numpy as np
import os
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
	print (subont)
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

def find_brain_genes():
	union=pd.read_csv('brain_genes_index.csv', index_col=[0])
	union=union['genes'].tolist()
	print ('brain_genes', len(union))
	return union


def find_positives(subont_genes, synDB_genes, synsysnet_genes, brain_genes):

	pos=list(set(subont_genes)&set(synDB_genes)&set(synsysnet_genes))

	print (len(pos))


	overlap=list(set(pos)&set(brain_genes))
	overlap.sort()

	print (len(overlap))
	print (overlap[:5])


	df=pd.DataFrame({'genes': overlap})
	df.to_csv('nb_brain_positive_pool.csv')

	return overlap

#Negative synapse examples----------------------------------------------------------------------------------------------------


def find_negatives(go_genes, subont_genes, synDB_genes, synsysnet_genes, brain_genes):

	negatives=list(set(go_genes)-set(subont_genes)-set(synDB_genes)-set(synsysnet_genes))

	overlap=list(set(negatives)&set(brain_genes))
	overlap.sort()

	print (len(overlap))
	print (overlap[:5])

	df=pd.DataFrame({'genes': overlap})

	df.to_csv('nb_brain_negative_pool.csv')

	return overlap


def find_pos_neg_pools():
	go_genes, subont_genes=GO_Focus('GO:0045202')


	#Extract all of the genes in the synapse in SynDB---------------------------------
	synDB_genes=find_synDB_genes()
	#Extract the genes in SynSysNet---------------------------------------------

	synsysnet_genes=find_synsysnet_genes()

	brain_genes=find_brain_genes()

	positives=find_positives(subont_genes, synDB_genes, synsysnet_genes, brain_genes)
	negatives=find_negatives(go_genes, subont_genes, synDB_genes, synsysnet_genes, brain_genes)
	return positives, negatives




def find_pos_neg_training(gene_list, name):
	random.seed(0)
	random.shuffle(gene_list)
	genes=gene_list[:200]

	print (gene_list[:5])
	
	df=pd.DataFrame(genes, columns=['genes'])

	df.to_csv('nb_brain_synapse_%s.csv'%name)

	
	return genes


if __name__ == '__main__':
	

	positives, negatives=find_pos_neg_pools()
	
	pos=find_pos_neg_training(positives, 'positives')
	neg=find_pos_neg_training(negatives, 'negatives')
	