#Goal: #1: analyze the new experimental results and compare the enrichment of SynGO and Synch
		#2: find the previously unvalidated fetal specific genes that are now validated by experiments
		#3: find the previously unvalidated adult specific genes that are now validated by experiments
		#4: draw the intersection between fetal and adult and predicted genes for figure
import pandas as pd
import numpy as np
import ddot
from ddot import Ontology
import csv
import sys
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

from scipy.stats import hypergeom

import matplotlib
from matplotlib import pyplot as plt

from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles

def load_pred_genes():
	pred_file='pred_genes_above_4.7.csv'
	pred_df=pd.read_csv(pred_file, index_col=[0])
	pred_genes=pred_df['genes'].tolist()
	pred_genes=[x.upper() for x in pred_genes]
	return pred_genes

def load_training_genes():
	filename='synapse_positives.csv'
	df=pd.read_csv(filename, index_col=[0])
	pos=df['genes'].tolist()

	filename='synapse_negatives.csv'
	df=pd.read_csv(filename, index_col=[0])
	neg=df['genes'].tolist()

	genes=list(set(pos+neg))
	return genes

def find_pred_no_training():
	pred_genes=load_pred_genes()
	training_genes=load_training_genes()
	pred_no_training=list(set(pred_genes)-set(training_genes))
	return pred_no_training

def load_fetal_brain():
	df=pd.read_csv('../experimental_validation/coba_fetal_brain.csv')
	#print (df)
	genes=df['Norm_Symbol'].tolist()
	training=load_training_genes()
	genes=[x.upper() for x in genes]
	genes=list(set(genes)-set(training))
	return genes

def load_ngn2():
	df=pd.read_csv('../experimental_validation/Coba_NGN2.csv')
	genes=df['Norm_Symbol'].tolist()
	training=load_training_genes()
	genes=[x.upper() for x in genes]
	genes=list(set(genes)-set(training))
	return genes

def find_hypergeometric(genes, pred_no_training):

	overlap=list(set(genes)&set(pred_no_training))
	M=10683
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
	print (len(genes))
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

def load_syngo_genes():
	syngo=Ontology.from_table('../prev_databases/SynGO_BP.txt')
	syngo_bp_genes=syngo.genes
	syngo=Ontology.from_table('../prev_databases/SynGO_CC.txt')
	syngo_cc_genes=syngo.genes
	syngo_genes=list(set(syngo_bp_genes+syngo_cc_genes))
	training=load_training_genes()
	syngo_genes=list(set(syngo_genes)-set(training))
	return syngo_genes

def find_GO_synapse():
	df=pd.read_csv('../prev_databases/GO_Synapse.csv')
	#print (df)
	genes=df['genes'].tolist()
	training=load_training_genes()
	genes=list(set(genes)-set(training))
	return genes

def load_fetal_data():
	fetal_brain=load_fetal_brain()
	ngn2=load_ngn2()
	fetal_overlap=list(set(fetal_brain)&set(ngn2))
	return fetal_brain, ngn2, fetal_overlap

def load_adult_data():
	adult_ctx=load_adult_ctx()
	adult_str=load_adult_str()
	adult_overlap=list(set(adult_ctx)&set(adult_str))
	return adult_ctx, adult_str, adult_overlap

def load_prev_databases():
	syngo=load_syngo_genes()
	synsysnet=find_synsysnet()
	synDB=find_synDB()
	go_synapse=find_GO_synapse()
	return syngo, synsysnet, synDB, go_synapse


#fetal datasets enrich for synsig
def find_fetal_enrichment():
	fetal_brain, ngn2, overlap=load_fetal_data()
	syngo, synsysnet, synDB, go_synapse=load_prev_databases()
	
	pred=load_pred_genes()
	print ('pred', len(pred))

	#enrich for synsig
	print ('fetal_brain, SynSig')
	find_hypergeometric(fetal_brain, pred)
	print ('ngn2, SynSig')
	find_hypergeometric(ngn2, pred)
	print ('overlap, SynSig')
	find_hypergeometric(overlap, pred)

	#fetal datasets enrich for previous databases:
	print ('fetal_brain, go_synapse')
	find_hypergeometric(fetal_brain, go_synapse)
	print ('ngn2, go_synapse')
	find_hypergeometric(ngn2, go_synapse)
	print ('overlap, go_synapse')
	find_hypergeometric(overlap, go_synapse)
	print ('fetal_brain, syngo')
	find_hypergeometric(fetal_brain, syngo)
	print ('ngn2, syngo')
	find_hypergeometric(ngn2, syngo)
	print ('overlap, syngo')
	find_hypergeometric(overlap, syngo)


def find_overlap_with_synsig():

	#load all datasets
	fetal_brain, ngn2, fetal_overlap=load_fetal_data()
	adult_ctx, adult_str, adult_overlap=load_adult_data()
	syngo, synsysnet, synDB, go_synapse=load_prev_databases()
	pred=load_pred_genes()

	adult_all=list(set(adult_ctx+adult_str))
	fetal_all=list(set(fetal_brain+ngn2))

	db=list(set(syngo+synsysnet+synDB+go_synapse))
	db_pred=list(set(db)&set(pred))
	#overlap_pred=list(set(overlap)&set(pred))
	adult_all_pred=list(set(adult_all)&set(pred))
	fetal_all_pred=list(set(fetal_all)&set(pred))
	return fetal_all_pred, adult_all_pred, db_pred

def find_fetal_synsig_only():
	fetal_all_pred, adult_all_pred, db_pred=find_overlap_with_synsig()
	fetal_only=list(set(fetal_all_pred)-set(adult_all_pred)-set(db_pred))
	print (len(fetal_only))
	df=pd.DataFrame({'genes': fetal_only})
	df.to_csv('fetal_only_val.csv')
	return fetal_only

def plot_fetal_adult_db_within_synsig():
	fetal_all_pred, adult_all_pred, db_pred=find_overlap_with_synsig()
	#plot the overlap of all fetal, all_adult and all database synapse genes:

	v=venn3([set(fetal_all_pred), set(adult_all_pred), set(db_pred)], set_labels=('Fetal Synapse', 'Adult Synapse', 'Synapse Databases'), set_colors=('red', 'gray', 'lightblue'), alpha=0.7)
	c=venn3_circles([set(fetal_all_pred), set(adult_all_pred), set(db_pred)], linestyle='solid', linewidth=0.5, color="black")
	for text in v.set_labels:
	    text.set_fontweight('bold')
	for text in v.set_labels:
	    text.set_fontsize(25)
	for text in v.subset_labels:
	    text.set_fontsize(25)

	plt.show()
	plt.close()


#find high-confidence fetal synsig genes (supported by both datasets and not in either of the other two datasets) and adult synsig genes (same logic)

def find_hc_fetal_adult_val():
	fetal_brain, ngn2, fetal_overlap=load_fetal_data()
	adult_ctx, adult_str, adult_overlap=load_adult_data()
	pred=load_pred_genes()

	fetal_overlap=list(set(fetal_brain)&set(ngn2)&set(pred))
	adult_overlap=list(set(adult_ctx)&set(adult_str)&set(pred))

	fetal_all_overlap=list(set(fetal_brain+ngn2)&set(pred))
	adult_all_overlap=list(set(adult_ctx+adult_str)&set(pred))

	syngo, synsysnet, synDB, go_synapse=load_prev_databases()
	db=list(set(syngo+synsysnet+synDB+go_synapse))

	fetal_specific_val=list(set(fetal_overlap)-set(adult_all_overlap)-set(db))
	print (len(fetal_specific_val), fetal_specific_val)

	adult_specific_val=list(set(adult_overlap)-set(fetal_all_overlap)-set(db))
	print (len(adult_specific_val), adult_specific_val)
	return fetal_specific_val, adult_specific_val

if __name__ == '__main__':

	#find how fetal synapse genes are enriched by different databses:
	find_fetal_enrichment()

	#find the synsig genes validated by only fetal experimental datasets:
	find_fetal_synsig_only()

	#plot the overlap between fetal, adult, database within synsig genes:
	plot_fetal_adult_db_within_synsig()

	#find high-confidence synsig genes specific to each period:
	find_hc_fetal_adult_val()