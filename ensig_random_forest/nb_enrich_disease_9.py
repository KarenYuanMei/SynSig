#Goal: determine how well does non-brain predicted synapse genes enrich for disease genes

import numpy as np
import pandas as pd
import csv

import matplotlib
matplotlib.use("TKAgg")
from matplotlib import pyplot as plt


plt.rcParams.update({'font.size': 22})
plt.rcParams["font.family"] = "Arial"

from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles

def load_training_genes():
	filename='/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/synapse_positives.csv'
	df=pd.read_csv(filename, index_col=[0])
	pos=df['genes'].tolist()

	filename='/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/synapse_negatives.csv'
	df=pd.read_csv(filename, index_col=[0])
	neg=df['genes'].tolist()

	genes=list(set(pos+neg))
	return genes

def load_brain_pred_genes():
	pred_file='/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/pred_genes_above_4.7.csv'
	pred_df=pd.read_csv(pred_file, index_col=[0])
	pred_genes=pred_df['genes'].tolist()
	pred_genes=[x.upper() for x in pred_genes]
	print ('brain', len(pred_genes))
	training=load_training_genes()

	pred_genes=list(set(pred_genes)-set(training))
	print ('brain', len(pred_genes))
	return pred_genes

def load_nonbrain_pred_genes():
	pred_file='nonbrain_pred_genes_above_4.67.csv'
	pred_df=pd.read_csv(pred_file, index_col=[0])
	pred_genes=pred_df['genes'].tolist()
	pred_genes=[x.upper() for x in pred_genes] 
	print (len(pred_genes))

	training=load_training_genes()

	pred_genes=list(set(pred_genes)-set(training))
	print (len(pred_genes))
	return pred_genes

def load_decipher_genes():
	dd=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Synapse_diseases/decipher.csv')
	dd=dd[dd['DDD category']=="confirmed"]
	genes=dd['gene symbol'].tolist()
	genes=list(set(genes))
	print ('decipher', len(genes))

	training=load_training_genes()

	genes=list(set(genes)-set(training))
	return genes

def find_decipher_syndromic():
	dd=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Synapse_diseases/decipher.csv')
	dd=dd[dd['DDD category']=="confirmed"]
	#print (dd)
	new=dd[dd['disease name'].str.contains("SYNDROME")]
	print (new)
	genes=new['gene symbol'].tolist()
	return genes

# def find_sfari_syndromic_genes():
# 	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Disease_genes/Autism/SFARI-Gene_genes.csv')
# 	#df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Disease_genes/Autism/SFARI_genes_051120.csv')
# 	print (df)
# 	autism=df['gene-symbol'].tolist()
# 	df=df[df['syndromic']>0]
# 	syndromic=df['gene-symbol'].tolist()

# 	training=load_training_genes()

# 	syndromic=list(set(syndromic)-set(training))
# 	return syndromic

# def find_sfari_nonsyndromic_genes():
# 	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Disease_genes/Autism/SFARI-Gene_genes.csv')
# 	#df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Disease_genes/Autism/SFARI_genes_051120.csv')
# 	print (df)
# 	autism=df['gene-symbol'].tolist()
# 	df=df[df['syndromic']==0]
# 	nonsyndromic=df['gene-symbol'].tolist()

# 	training=load_training_genes()

# 	nonsyndromic=list(set(nonsyndromic)-set(training))
# 	return nonsyndromic


def find_sfari_syndromic_genes():
	#df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Disease_genes/Autism/SFARI_genes_051120.csv')
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Disease_genes/Autism/SFARI-Gene_genes.csv')
	
	syn = df[df['genetic-category'].str.contains('Syndromic', regex=False, case=False, na=False)]
	print (syn)
	genes=syn['gene-symbol'].tolist()


	df=df[df['syndromic']==1]
	genes2=df['gene-symbol'].tolist()

	print (len(genes), len(genes2), len(list(set(genes)&set(genes2))))

	genes=list(set(genes+genes2))

	training=load_training_genes()

	genes=list(set(genes)-set(training))
	
	return genes

def find_sfari_nonsyndromic_genes():
	#df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Disease_genes/Autism/SFARI_genes_051120.csv')
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Disease_genes/Autism/SFARI-Gene_genes.csv')
	
	df = df[~df['genetic-category'].str.contains('Syndromic', regex=False, case=False, na=False)]
	print (df)
	genes=df['gene-symbol'].tolist()

	df=df[df['syndromic']==0]
	genes2=df['gene-symbol'].tolist()

	print (len(genes), len(genes2), len(list(set(genes)&set(genes2))))

	genes=list(set(genes+genes2))

	training=load_training_genes()

	genes=list(set(genes)-set(training))
	
	return genes


def find_nonbrain_disease_enrichment():

	nb_pred=load_nonbrain_pred_genes()
	print (len(nb_pred))

	dd_syn=load_decipher_genes()

	#how well does nonbrain predicted genes enrich for general developmental disease genes:
	venn2([set(nb_pred), set(dd_syn)], set_labels = ('Non-Brain', "Decipher Genes"))
	plt.show()
	plt.close()

	#how well does nonbrain predicted genes enrichment for autism nonsyndromic:
	nonsyndromic=find_sfari_nonsyndromic_genes()
	print ('sfari', len(nonsyndromic))

	#how well does nonbrain predicted genes enrichment for autism syndromic:
	aut_syndromic=find_sfari_syndromic_genes()
	print ('autism syndromic', len(aut_syndromic))

	venn2([set(nb_pred), set(aut_syndromic)], set_labels = ('NonBrain', "aut_syndromic"))
	plt.show()
	plt.close()

	overlap=list(set(nb_pred)&set(aut_syndromic))
	print ('overlap with nb and syndromic', overlap)

	venn2([set(nb_pred), set(nonsyndromic)], set_labels=('Non-Brain Predicted Synapse', 'Nonsyndromic'))
	plt.show()
	plt.close()


if __name__ == '__main__':
	load_brain_pred_genes()
	find_nonbrain_disease_enrichment()
