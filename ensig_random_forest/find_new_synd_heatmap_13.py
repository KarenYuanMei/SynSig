#Goal: find the systems that syndromic autism genes enrich for (Fig 7f)

import pandas as pd
#import networkx as nx
import numpy as np

import csv

#import sys
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

#sys.path.append("/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_8/")
from scipy.stats import hypergeom

import ddot
from ddot import Ontology

from mlxtend.evaluate import permutation_test
import random
from statsmodels.sandbox.stats.multicomp import multipletests


def load_new_synd():
	df=pd.read_csv('nb_val_new_syndromic.csv')
	#print (df)
	genes=df['Genes'].tolist()
	return genes

def load_pheno_hpo():
	hpo=Ontology.from_table('/Users/karenmei/Documents/Synapse_Ontology/HPO/making_HPO/HPO_parent_child.txt')
	#print (hpo)
	hpo=hpo.propagate(direction='forward', gene_term=True, term_term=False)
	hpo=hpo.focus(branches=['Phenotypic abnormality'])
	return hpo

def find_gene_phenotypes(genes):
	hpo=load_pheno_hpo()
	phe_types=hpo.parent_2_child['Phenotypic abnormality']

	all_gene_phe=[]
	for item in phe_types:
		gene_phe=[]
		for gene in genes:
			subont=hpo.focus(branches=item)
			if gene in subont.genes:
				entry=1
			else:
				entry=0
			gene_phe.append(entry)
		all_gene_phe.append(gene_phe)

	mat=np.array(all_gene_phe)
	df=pd.DataFrame(mat, index=phe_types, columns=genes)
	df=df.T
	df.to_csv('synapse_phenotype_heatmap.csv')
	return df
genes=load_new_synd()
find_gene_phenotypes(genes)