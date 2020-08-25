#Goal: find the systems that syndromic autism genes enrich for (Fig 7c)

import pandas as pd
#import networkx as nx
import numpy as np

import csv

#import sys
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

#sys.path.append("/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_8/")
from scipy.stats import hypergeom

# import matplotlib
# matplotlib.use("TKAgg")
# print(matplotlib.get_backend())
# from matplotlib import pyplot as plt

# from matplotlib_venn import venn3, venn3_circles
# from matplotlib_venn import venn2, venn2_circles


import ddot
from ddot import Ontology

from mlxtend.evaluate import permutation_test
import random
from statsmodels.sandbox.stats.multicomp import multipletests


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

def load_pheno_hpo():
	hpo=Ontology.from_table('/Users/karenmei/Documents/Synapse_Ontology/HPO/making_HPO/HPO_parent_child.txt')
	#print (hpo)
	hpo=hpo.propagate(direction='forward', gene_term=True, term_term=False)
	hpo=hpo.focus(branches=['Phenotypic abnormality'])
	return hpo

def find_hypergeometric(subont_gene_no, query_list, overlap, hpo):
	query_gene_no=list(set(query_list)&set(hpo.genes))
	print (len(query_gene_no))
	M=len(hpo.genes)
	print (M)
	#M=20000
	N=subont_gene_no
	n=len(query_gene_no)
	x=overlap
	pval = hypergeom.sf(x-1, M, n, N)

	rv = hypergeom(M, n, N)
	distr = np.arange(0, n+1)
	prob = rv.pmf(distr)
	maximum=np.max(prob)
	result = np.where(prob == maximum)
	result=result[0]
	fold=x/result
	fold=fold.tolist()
	print ('Fold Enrichment', fold)
	print ('hypergeometric p-value', pval)
	return fold, pval

def enrich_phe_no_df(genes):
	hpo=load_pheno_hpo()
	phe_types=hpo.parent_2_child['Phenotypic abnormality']

	print (phe_types)
	phe_overlap_list=[]
	subont_tot_genes=[]
	fold_list=[]
	pval_list=[]
	for item in phe_types:
		print (item)
		ont=hpo.focus(branches=[item])
		ont_genes=ont.genes
		overlap=list(set(genes)&set(ont_genes))
		#percent=len(overlap)/len(overlap_hpo)*100
		phe_overlap_list.append(len(overlap))
		subont_tot_genes.append(len(ont.genes))
		fold, pval=find_hypergeometric(len(ont.genes), genes, len(overlap), hpo)
		fold_list.append(fold[0])
		pval_list.append(pval)

	df=pd.DataFrame({'Phenotypes': phe_types, 'Pheno_Overlap': phe_overlap_list, 'Phe_All': subont_tot_genes, 'Fold_Enrichment': fold_list, 'PVal_List': pval_list})
	#print (df)
	df=df.set_index('Phenotypes')
	print (df)
	#df.to_csv('%s_pheno_no_df'%name)
	return df

def corr_p(df):
	plist=df['PVal_List'].tolist()
	corr=multipletests(plist, alpha=0.05, method='fdr_bh')
	df['Corr_P']=corr[1]
	df['Corr_True_False']=corr[0]
	print (df)
	return df

if __name__ == '__main__':
	syndromic=find_sfari_syndromic_genes()
	df=enrich_phe_no_df(syndromic)
	corr_df=corr_p(df)
	corr_df.to_csv('syndromic_pheno_no_df')



	
