#Goal: compare syndromic vs. nonsyndromic phenotype number (Fig 7b)

import pandas as pd
#import networkx as nx
import numpy as np

import csv

#import sys
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

#sys.path.append("/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_8/")
from scipy.stats import hypergeom

import matplotlib
matplotlib.use("TKAgg")
print(matplotlib.get_backend())
from matplotlib import pyplot as plt

from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles


import ddot
from ddot import Ontology

from mlxtend.evaluate import permutation_test
import random
#from statsmodels.sandbox.stats.multicomp import multipletests

def load_nonbrain_pred_genes():
	pred_file='nonbrain_pred_genes_above_4.67.csv'
	pred_df=pd.read_csv(pred_file, index_col=[0])
	pred_genes=pred_df['genes'].tolist()
	pred_genes=[x.upper() for x in pred_genes]
	return pred_genes

def find_sfari_syndromic_genes():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Disease_genes/Autism/SFARI-Gene_genes.csv')
	syn = df[df['genetic-category'].str.contains('Syndromic', regex=False, case=False, na=False)]
	genes=syn['gene-symbol'].tolist()
	df=df[df['syndromic']==1]
	genes2=df['gene-symbol'].tolist()
	genes=list(set(genes+genes2))
	return genes

def find_sfari_nonsyndromic_genes():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Disease_genes/Autism/SFARI-Gene_genes.csv')
	syndromic=find_sfari_syndromic_genes()
	all_genes=df['gene-symbol'].tolist()
	genes=list(set(all_genes)-set(syndromic))
	print (len(genes))
	return genes

def load_synapse_negatives():
	synapse_genes=pd.read_csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/synapse_negatives.csv')
	synapse_negatives=synapse_genes['genes'].tolist()
	return synapse_negatives

def load_negative_pool():
	pred=load_nonbrain_pred_genes()
	training_neg=load_synapse_negatives()
	neg=pd.read_csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/negative_pool.csv')
	neg=neg['genes'].tolist()
	neg=list(set(neg)-set(pred)-set(training_neg))
	return neg

def find_phe_no(hpo, list1, list2):
	phe_list1=[]
	for item in list1:
		if item in hpo.genes:
			#pos.append(item)
			phe_no=len(hpo.gene_2_term[item])
			phe_list1.append(phe_no)

	print (np.mean(phe_list1))
	print (len(phe_list1))

	phe_list2=[]
	for item in list2:
		if item in hpo.genes:
			phe_no=len(hpo.gene_2_term[item])
			phe_list2.append(phe_no)

	random.seed(4)
	random.shuffle(phe_list2)

	print (np.mean(phe_list2))
	#random.shuffle(phe_list1)
	phe_list2=phe_list2[:len(phe_list1)]
	print (len(phe_list2))
	return phe_list1, phe_list2

def load_pheno_hpo():
	hpo=Ontology.from_table('/Users/karenmei/Documents/Synapse_Ontology/HPO/making_HPO/HPO_parent_child.txt')
	#print (hpo)
	hpo=hpo.propagate(direction='forward', gene_term=True, term_term=False)
	hpo=hpo.focus(branches=['Phenotypic abnormality'])
	return hpo

def hpo_delete_nervous(hpo):
	nervous=hpo.focus(branches='Abnormality of the nervous system')
	#hpo=hpo.focus(branches=['Phenotypic abnormality'])
	hpo=hpo.delete(to_delete=nervous.terms)
	return hpo

def hpo_delete_nervous_head(hpo):
	nervous=hpo.focus(branches='Abnormality of the nervous system')
	head_neck=hpo.focus(branches='Abnormality of head or neck')
	#hpo=hpo.focus(branches=['Phenotypic abnormality'])
	delete_terms=list(set(head_neck.terms+nervous.terms))
	hpo=hpo.delete(to_delete=delete_terms)
	return hpo

def make_pheno_df(phe_list1, phe_list2):
	p_value = permutation_test(phe_list1, phe_list2,
	                           method='approximate',
	                           num_rounds=10000,
	                           seed=0)
	print(p_value)
	df=pd.DataFrame({'Group': 'Syndromic', 'Phe_No': phe_list1})
	df=df.set_index('Group')
	df2=pd.DataFrame({'Group': 'Non-Syndromic', 'Phe_No': phe_list2})
	df2=df2.set_index('Group')

	final=pd.concat([df, df2], axis=0)
	return final
	#print (final)
	
	#final.to_csv('pred_neg_nb_R_new.csv')

#make figure 7b:
def compare_syndromic_ex_nervous():
	hpo=load_pheno_hpo()
	hpo=hpo_delete_nervous(hpo)
	syndromic=find_sfari_syndromic_genes()
	syndromic.sort()
	nonsyndromic=find_sfari_nonsyndromic_genes()
	nonsyndromic.sort()
	phe_list1, phe_list2=find_phe_no(hpo, syndromic, nonsyndromic)
	#plot_pheno_diff(phe_list1, phe_list2)
	final=make_pheno_df(phe_list1, phe_list2)
	final.to_csv('syndromic_nonsyndromic_nb_pred_R.csv')

#make supplemental figure:

def compare_syndromic_ex_nervous_head():
	hpo=load_pheno_hpo()
	hpo=hpo_delete_nervous_head(hpo)
	syndromic=find_sfari_syndromic_genes()
	syndromic.sort()
	nonsyndromic=find_sfari_nonsyndromic_genes()
	nonsyndromic.sort()
	phe_list1, phe_list2=find_phe_no(hpo, syndromic, nonsyndromic)
	#plot_pheno_diff(phe_list1, phe_list2)
	final=make_pheno_df(phe_list1, phe_list2)
	#final.to_csv('syndromic_nonsyndromic_nohead_R.csv')

#make supplemental figure:

def compare_synapse_ex_nervous():
	hpo=load_pheno_hpo()
	hpo=hpo_delete_nervous(hpo)
	nb=load_nonbrain_pred_genes()
	print (len(nb))
	nb.sort()
	neg=load_negative_pool()
	neg.sort()
	print (neg[:5])
	print (len(neg))
	phe_list1, phe_list2=find_phe_no(hpo, nb, neg)
	#plot_pheno_diff(phe_list1, phe_list2)
	final=make_pheno_df(phe_list1, phe_list2)
	#final.to_csv('pred_neg_nb_R_new.csv')

def compare_synapse_ex_nervous_head():
	hpo=load_pheno_hpo()
	hpo=hpo_delete_nervous(hpo)
	nb=load_nonbrain_pred_genes()
	print (len(nb))
	nb.sort()
	neg=load_negative_pool()
	neg.sort()
	print (neg[:5])
	print (len(neg))
	phe_list1, phe_list2=find_phe_no(hpo, nb, neg)
	#plot_pheno_diff(phe_list1, phe_list2)
	final=make_pheno_df(phe_list1, phe_list2)
	#final.to_csv('pred_neg_nb_nohead_R.csv')

if __name__ == '__main__':
	compare_syndromic_ex_nervous()
	compare_syndromic_ex_nervous_head()
	compare_synapse_ex_nervous()
	compare_synapse_ex_nervous_head()


	