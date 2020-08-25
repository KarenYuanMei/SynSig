#Goal: Fig6c

import csv
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("TKAgg")
print(matplotlib.get_backend())
from matplotlib import pyplot as plt
from scipy import stats
from sklearn import preprocessing
import random


from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles
#import venn
import ddot
from ddot import Ontology
from scipy.stats import hypergeom


def load_training_genes():
	filename='/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/synapse_positives.csv'
	df=pd.read_csv(filename, index_col=[0])
	pos=df['genes'].tolist()

	filename='/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/synapse_negatives.csv'
	df=pd.read_csv(filename, index_col=[0])
	neg=df['genes'].tolist()

	training_genes=list(set(pos+neg))
	return training_genes

def load_nb_genes():
	nb=pd.read_csv('nonbrain_pred_genes_above_4.67.csv', index_col=[0])
	print (nb)
	pred_genes=nb['genes'].tolist()
	genes=[x.upper() for x in pred_genes]

	training=load_training_genes()
	pred_genes=list(set(genes)-set(training))
	return pred_genes


def load_adult_ctx():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Validation_proteomics/Weijun_proteomics/weijun_ctx_uniprot.csv', sep='\t')
	#print (df)
	genes=df['To'].tolist()
	training=load_training_genes()
	genes=[x.upper() for x in genes]
	genes=list(set(genes)-set(training))
	return genes

def load_adult_str():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Validation_proteomics/Weijun_proteomics/weijun_str_uniprot.csv', sep='\t')
	#print (df)
	genes=df['To'].tolist()
	training=load_training_genes()
	genes=[x.upper() for x in genes]
	genes=list(set(genes)-set(training))
	return genes

def load_fetal_brain():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Validation_proteomics/Coba_human_fetal_2020/coba_fetal_brain.csv')
	#print (df)
	genes=df['Norm_Symbol'].tolist()
	training=load_training_genes()
	genes=[x.upper() for x in genes]
	genes=list(set(genes)-set(training))
	return genes

def load_ngn2():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Validation_proteomics/Coba_NGN2_2020/Coba_NGN2.csv')
	genes=df['Norm_Symbol'].tolist()
	training=load_training_genes()
	genes=[x.upper() for x in genes]
	genes=list(set(genes)-set(training))
	return genes

def load_syngo_genes():
	syngo=Ontology.from_table('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Metrics/SynGO_BP.txt')
	syngo_bp_genes=syngo.genes
	syngo=Ontology.from_table('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Metrics/SynGO_CC.txt')
	syngo_cc_genes=syngo.genes
	syngo_genes=list(set(syngo_bp_genes+syngo_cc_genes))

	training=load_training_genes()
	syngo_genes=list(set(syngo_genes)-set(training))
	return syngo_genes


def find_synsysnet():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Synapse_Genes/SynSysNet_genes.csv')
	#print (df)
	genes=df['gene_name'].tolist()
	print (len(genes))

	training=load_training_genes()
	genes=list(set(genes)-set(training))
	return genes

def find_synDB():
	df=pd.read_csv('/Users/karenmei/Documents/BrainHierarchyDataSource/SynDB_Master.csv')
	print (df)
	genes=df['Symbol'].tolist()

	training=load_training_genes()
	genes=list(set(genes)-set(training))
	return genes

def find_GO_synapse():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Synapse_Scripts/GO_Synapse.csv')
	print (df)
	genes=df['genes'].tolist()

	training=load_training_genes()
	genes=list(set(genes)-set(training))
	return genes


def plot_overlap_nb():
	syngo=load_syngo_genes()
	synsysnet=find_synsysnet()
	synDB=find_synDB()
	go_synapse=find_GO_synapse()

	db=list(set(syngo+synsysnet+synDB+go_synapse))


	nb=load_nb_genes()
	#pred=load_pred_genes()

	print (len(nb))

	#overlap=list(set(nb)&set(pred))
	#print (len(overlap))

	adult=list(set(load_adult_str())&set(load_adult_str()))
	fetal=list(set(load_fetal_brain())&set(load_ngn2()))

	#adult=list(set(load_adult_str()+load_adult_str()))
	#fetal=list(set(load_fetal_brain()+load_ngn2()))

	val=adult+fetal

	v=venn3([set(nb),set(val), set(db)], set_labels=('Non-Brain Predicted',  'Proteomics Screen', 'Synapse Databases'), set_colors=('dodgerblue', 'red', 'lightgray'),alpha=0.7)
	venn3_circles([set(nb),set(val), set(db)], linestyle='solid', linewidth=0.5, color='k');
	for text in v.set_labels:
	    text.set_fontweight('bold')
	for text in v.set_labels:
	    text.set_fontsize(30)
	for text in v.subset_labels:
	    text.set_fontsize(30)
	plt.show()
	plt.close()

def find_hypergeometric(genes, pred_no_training):
	overlap=list(set(genes)&set(pred_no_training))
	M=10682
	#M=20000
	N=len(genes)
	n=len(pred_no_training)
	x=len(overlap)
	pval = hypergeom.sf(x-1, M, n, N)

	rv = hypergeom(M, n, N)
	distr = np.arange(0, n+1)
	#print (N, n, x)
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

if __name__ == '__main__':
	plot_overlap_nb()

	syngo=load_syngo_genes()
	nb=load_nb_genes()
	find_hypergeometric(syngo, nb)

	stria=load_adult_str()
	ctx=load_adult_ctx()
	adult=list(set(stria+ctx))
	find_hypergeometric(adult, nb)

	fetal_brain=load_fetal_brain()
	ngn2=load_ngn2()
	fetal=list(set(fetal_brain+ngn2))
	find_hypergeometric(fetal, nb)




