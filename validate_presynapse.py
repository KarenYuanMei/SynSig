import ddot
from ddot import Ontology

import pandas as pd
#import networkx as nx
import numpy as np

def load_syngo_genes():
	syngo=Ontology.from_table('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Metrics/SynGO_BP.txt')
	syngo_bp_genes=syngo.genes
	syngo=Ontology.from_table('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Metrics/SynGO_CC.txt')
	syngo_cc_genes=syngo.genes
	#syngo_genes=list(set(syngo_bp_genes+syngo_cc_genes))
	#training=load_training_genes()
	#syngo_genes=list(set(syngo_genes)-set(training))
	return syngo_bp_genes, syngo_cc_genes

def load_syngo_presynapse():
	syngo=Ontology.from_table('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Metrics/SynGO_CC.txt')
	syngo=syngo.propagate(direction='forward', gene_term=True, term_term=False)
	pre=syngo.focus(branches=['presynapse'])
	pre_genes=pre.genes
	return pre_genes


def load_pred_genes():
	pred_file='pred_genes_above_4.7.csv'
	pred_df=pd.read_csv(pred_file, index_col=[0])
	pred_genes=pred_df['genes'].tolist()
	pred_genes=[x.upper() for x in pred_genes]
	return pred_genes

def load_adult_ctx():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Validation_proteomics/Weijun_proteomics/weijun_ctx_uniprot.csv', sep='\t')
	#print (df)
	genes=df['To'].tolist()
	#training=load_training_genes()
	genes=[x.upper() for x in genes]
	#genes=list(set(genes)-set(training))
	return genes

def load_adult_str():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Validation_proteomics/Weijun_proteomics/weijun_str_uniprot.csv', sep='\t')
	#print (df)
	genes=df['To'].tolist()
	#training=load_training_genes()
	genes=[x.upper() for x in genes]
	#genes=list(set(genes)-set(training))
	return genes

def load_fetal_brain():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Validation_proteomics/Coba_human_fetal_2020/coba_fetal_brain.csv')
	print (df)
	genes=df['Norm_Symbol'].tolist()
	#training=load_training_genes()
	genes=[x.upper() for x in genes]
	#genes=list(set(genes)-set(training))
	return genes

def load_ngn2():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Validation_proteomics/Coba_NGN2_2020/Coba_NGN2.csv')
	genes=df['Norm_Symbol'].tolist()
	#training=load_training_genes()
	genes=[x.upper() for x in genes]
	#genes=list(set(genes)-set(training))
	return genes

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

bp, cc=load_syngo_genes()
print (len(bp), len(cc))

pre_genes=load_syngo_presynapse()
print (len(pre_genes))

ctx=load_adult_ctx()
overlap=list(set(pre_genes)&set(ctx))
print (len(overlap))

ctx=load_adult_str()
overlap=list(set(pre_genes)&set(ctx))
print (len(overlap))

ctx=load_adult_ctx()
overlap=list(set(pre_genes)&set(ctx))
print (len(overlap))

ctx=load_adult_str()
overlap=list(set(pre_genes)&set(ctx))
print (len(overlap))