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

import sys
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

#sys.path.append("/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_8/")
from scipy.stats import hypergeom

import matplotlib
#matplotlib.use("TKAgg")
#print(matplotlib.get_backend())
from matplotlib import pyplot as plt

from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles
#import venn

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
	df=pd.read_csv('/../experimental_validation/coba_fetal_brain.csv')
	#print (df)
	genes=df['Norm_Symbol'].tolist()
	training=load_training_genes()
	genes=[x.upper() for x in genes]
	genes=list(set(genes)-set(training))
	return genes

def load_ngn2():
	df=pd.read_csv('/experimental_validation/Coba_NGN2.csv')
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
	df=pd.read_csv('../experimental_validation/SynSysNet_genes.csv')
	#print (df)
	genes=df['gene_name'].tolist()
	print (len(genes))
	training=load_training_genes()
	genes=list(set(genes)-set(training))
	return genes

def find_synDB():
	df=pd.read_csv('../experimental_validation/SynDB_Master.csv')
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



fetal_brain=load_fetal_brain()
#print (len(fetal_brain))
ngn2=load_ngn2()
#print (len(ngn2))
overlap=list(set(fetal_brain)&set(ngn2))
fetal_all=list(set(fetal_brain+ngn2))

pred=load_pred_genes()
print ('pred', len(pred))

print ('fetal_brain, SynSig')
find_hypergeometric(fetal_brain, pred)

print ('ngn2, SynSig')
find_hypergeometric(ngn2, pred)

print ('overlap, SynSig')
find_hypergeometric(overlap, pred)


adult_ctx=load_adult_ctx()
adult_str=load_adult_str()
adult=list(set(adult_ctx)&set(adult_str))
adult_all=list(set(adult_ctx+adult_str))

syngo=load_syngo_genes()
synsysnet=find_synsysnet()
synDB=find_synDB()
go_synapse=find_GO_synapse()

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


db=list(set(syngo+synsysnet+synDB+go_synapse))
db=list(set(db)&set(pred))
overlap_pred=list(set(overlap)&set(pred))
adult_pred=list(set(adult_all)&set(pred))
fetal_all_pred=list(set(fetal_all)&set(pred))

fetal_only=list(set(fetal_all)-set(adult_all)-set(db))
fetal_only_val=list(set(fetal_only)&set(pred))

#print (fetal_only_val)
df=pd.DataFrame({'genes': fetal_only_val})
df.to_csv('fetal_only_val.csv')

#plot the overlap of all fetal, all_adult and all database synapse genes:

v=venn3([set(fetal_all_pred), set(adult_pred), set(db)], set_labels=('Fetal Synapse', 'Adult Synapse', 'Synapse Databases'), set_colors=('red', 'gray', 'lightblue'), alpha=0.7)
c=venn3_circles([set(fetal_all_pred), set(adult_pred), set(db)], linestyle='solid', linewidth=0.5, color="black")
for text in v.set_labels:
    text.set_fontweight('bold')
for text in v.set_labels:
    text.set_fontsize(25)
for text in v.subset_labels:
    text.set_fontsize(25)

plt.show()
plt.close()


fetal_overlap=list(set(fetal_brain)&set(ngn2)&set(pred))
adult_overlap=list(set(adult_ctx)&set(adult_str)&set(pred))

fetal_all_overlap=list(set(fetal_brain+ngn2)&set(pred))
adult_all_overlap=list(set(adult_ctx+adult_str)&set(pred))

fetal_specific_val=list(set(fetal_overlap)-set(adult_all_overlap)-set(db))
print (len(fetal_specific_val), fetal_specific_val)

adult_specific_val=list(set(adult_overlap)-set(fetal_all_overlap)-set(db))
print (len(adult_specific_val), adult_specific_val)
