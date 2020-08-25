#GOAL: find the GO semantics with synapse branch embedded in GO, with both positive and negative synapse examples

#make sure to do this on the server, because this takes a lot of memory

import pandas as pd
import networkx as nx
import numpy as np
import os
import ddot
from ddot import Ontology
import csv

#There are different GO score matrices:
#1) All_GO_score_matrix_for_synapse.csv: 1177 synapse genes by 1177 synapse genes from all GO
#2) GO_score_matrix.csv: 200 positive and 200 negative examples from all GO
#3) GO_synapse_score_matrix.csv: 1177 synapse genes scores from just GO_synapse
#4) syngo_score_matrix.csv: genes from syngo_cc collapsed into semantic similarity scores


def get_positive_gene_names(positive_filename):
	#function that returns a list of training gene names (strings)
	genes=pd.read_csv(positive_filename)
	gene_names=genes['genes'].tolist()
	return gene_names

def get_negative_gene_names(negative_filename):
	#function that returns a list of training gene names (strings)
	genes=pd.read_csv(negative_filename)
	gene_names=genes['genes'].tolist()
	return gene_names

def find_GO_score_matrix():
	ndex_server ='http://public.ndexbio.org' 
	ndex_user, ndex_pass = 'ym2', 'Synapse'
	go_human = Ontology.from_ndex('http://public.ndexbio.org/v2/network/16565851-4bfc-11e9-9f06-0ac135e8bacf')

	print (go_human)
	sim, genes=go_human.flatten()

	sim_df = pd.DataFrame(sim, index=genes, columns=genes)
	return sim_df

def find_input_gene_GO_scores(positive_filename, negative_filename):
	positive_gene_names=get_positive_gene_names(positive_filename)
	negative_gene_names=get_negative_gene_names(negative_filename)
	input_gene_names=list(set(positive_gene_names+negative_gene_names))
	GO_sim=find_GO_score_matrix()
	GO_score_matrix=GO_sim.loc[input_gene_names, input_gene_names]
	GO_score_matrix.to_csv('GO_training_score_matrix_for_big_pool_genes.csv')
	return GO_score_matrix


positive_filename='synapse_positives.csv'
negative_filename='synapse_negatives.csv'

find_input_gene_GO_scores(positive_filename, negative_filename)
