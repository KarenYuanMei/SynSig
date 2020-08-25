#Goal: Generate Fig 3d

import pandas as pd
#import networkx as nx
import numpy as np
#import os
import ddot
from ddot import Ontology
from collections import defaultdict
import collections
import csv

import sys
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

sys.path.append("/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_8/")
from scipy.stats import hypergeom

import matplotlib
matplotlib.use("TKAgg")
print(matplotlib.get_backend())
from matplotlib import pyplot as plt

from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
from matplotlib_venn import venn2, venn2_circles
import venn

def load_training_genes():
	filename='synapse_positives.csv'
	df=pd.read_csv(filename, index_col=[0])
	genes=df['genes'].tolist()
	return genes

def load_adult_ctx():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Validation_proteomics/Weijun_proteomics/weijun_ctx_uniprot.csv', sep='\t')
	print (df)
	genes=df['To'].tolist()
	training=load_training_genes()
	genes=[x.upper() for x in genes]
	genes=list(set(genes)-set(training))
	return genes

def load_adult_str():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Validation_proteomics/Weijun_proteomics/weijun_str_uniprot.csv', sep='\t')
	print (df)
	genes=df['To'].tolist()
	training=load_training_genes()
	genes=[x.upper() for x in genes]
	genes=list(set(genes)-set(training))
	return genes

def load_pred_genes():
	pred_file='pred_genes_above_4.7.csv'
	pred_df=pd.read_csv(pred_file, index_col=[0])
	pred_genes=pred_df['genes'].tolist()
	pred_genes=[x.upper() for x in pred_genes]
	return pred_genes

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

def load_syngo_genes():
	syngo=Ontology.from_table('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Metrics/SynGO_BP.txt')
	syngo_bp_genes=syngo.genes
	syngo=Ontology.from_table('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Metrics/SynGO_CC.txt')
	syngo_cc_genes=syngo.genes
	syngo_genes=list(set(syngo_bp_genes+syngo_cc_genes))
	training=load_training_genes()
	syngo_genes=list(set(syngo_genes)-set(training))
	return syngo_genes

def find_GO_synapse():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Synapse_Scripts/GO_Synapse.csv')
	print (df)
	genes=df['genes'].tolist()
	training=load_training_genes()
	genes=list(set(genes)-set(training))
	return genes

df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/brain_RNA_big_gene_pool_pipeline/exp_val.csv', index_col=[0])

print (df)

lit = df.drop(df[df.columns[-4:]], axis=1)

lit=lit.set_index('Genes')

lit['sum']=lit.sum(axis=1)
val=lit[lit['sum']>=1]

print (val)

val_genes=list(val.index)

ctx=load_adult_ctx()
stria=load_adult_str()
adult=list(set(ctx+stria))

pred=load_pred_genes()

v=venn3([set(adult), set(val_genes), set(pred)], set_labels=('Adult Brain \n Experimental Validation', 'Previous Screens Validation', 'Predicted'), set_colors=('coral', 'skyblue', 'lightgreen'), alpha = 0.7)
for text in v.set_labels:
    text.set_fontweight('bold')
for text in v.set_labels:
    text.set_fontsize(30)
for text in v.subset_labels:
    text.set_fontsize(30)

plt.show()
plt.close()

new=list(set(adult)-set(val_genes))
new_pred_val=list(set(new)&set(pred))

print (new_pred_val)

unval=list(set(pred)-set(adult))
val_by_lit=list(set(unval)&set(val_genes))
print (len(val_by_lit))


syngo=load_syngo_genes()
synDB=find_synDB()
synsysnet=find_synsysnet()
go_syn=find_GO_synapse()

db=list(set(syngo+synDB+synsysnet+go_syn))

v=venn3([set(adult), set(db), set(pred)], set_labels=('Adult Brain \n Synapse Validation', 'Synapse Databases', 'SynSig'), set_colors=('coral', 'skyblue', 'lightgreen'), alpha = 0.7)
for text in v.set_labels:
    text.set_fontweight('bold')
for text in v.set_labels:
    text.set_fontsize(25)
for text in v.subset_labels:
    text.set_fontsize(25)

plt.show()
plt.close()

v=venn3_unweighted([set(adult), set(db), set(pred)], set_labels=('Adult Brain \n Synapse Validation', 'Synapse Databases', 'SynSig'), set_colors=('gray', 'lightgray', 'red'), alpha = 0.7)
for text in v.set_labels:
    text.set_fontweight('bold')
for text in v.set_labels:
    text.set_fontsize(25)
for text in v.subset_labels:
    text.set_fontsize(25)

plt.show()
plt.close()

v=venn3([set(adult), set(db), set(pred)], set_labels=('Adult Brain \n Synapse Validation', 'Synapse Databases', 'SynSig'), set_colors=('gray', 'lightgray', 'red'), alpha = 0.7)
for text in v.set_labels:
    text.set_fontweight('bold')
for text in v.set_labels:
    text.set_fontsize(25)
for text in v.subset_labels:
    text.set_fontsize(25)

plt.show()
plt.close()


overlap=list(set(adult)&set(pred))
new=list(set(overlap)-set(db))
print (new)
new_adult=pd.DataFrame({'genes': new})
new_adult.to_csv('new_adult_pred.csv')

#find the new validated proteins outside of dabatases:

not_db=list(set(adult)-set(db))
val_not_db=list(set(not_db)&set(pred))
#print (val_not_db)

old=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/validate_new_exp/adult/new_adult_pred.csv')
old=old['genes'].tolist()

overlap=list(set(val_not_db)&set(old))

print (len(overlap))
def find_sfari_genes():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Disease_genes/Autism/SFARI-Gene_genes.csv')
	#print (df)
	autism=df['gene-symbol'].tolist()
	return autism

autism=find_sfari_genes()
autism_overlap=list(set(val_not_db)&set(autism))

print (autism_overlap)
