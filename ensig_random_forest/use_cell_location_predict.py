#Goal: to produce Fig6D

import csv
import math
import numpy as np

from scipy.stats.stats import pearsonr  

from mlxtend.evaluate import permutation_test

import pandas as pd

from collections import Counter
from collections import defaultdict

import matplotlib
matplotlib.use("TKAgg")
from matplotlib import pyplot as plt

import seaborn as sns; sns.set()
from scipy import stats

from numpy.random import seed 
from numpy.random import randn 
from scipy.stats import mannwhitneyu 

import ddot
from ddot import Ontology

import random
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix

plt.style.use('seaborn-deep')



def load_cell_locations_df():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/cell_location/Thul_cell_location.csv')
	df=df[df['Reliability']!='Uncertain']

	#print (df)
	df=df.set_index('Gene')
	return df

def load_syngo_genes():
	syngo=Ontology.from_table('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Metrics/SynGO_BP.txt')
	syngo_bp_genes=syngo.genes
	syngo=Ontology.from_table('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Metrics/SynGO_CC.txt')
	syngo_cc_genes=syngo.genes
	syngo_genes=list(set(syngo_bp_genes+syngo_cc_genes))
	return syngo_genes


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
	#print (df)
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


def load_synapse_positives():
	
	adult_ctx=load_adult_ctx()
	adult_str=load_adult_str()
	fetal_brain=load_fetal_brain()
	ngn2=load_ngn2()

	overlap=list(set(adult_ctx)&set(adult_str)&set(fetal_brain)&set(ngn2))
	synapse_genes=overlap

	#synapse_genes=load_syngo_genes()

	return synapse_genes

#load the negative synapse genes (training)
def load_synapse_negatives():
	synapse_genes=load_synapse_positives()
	negatives=pd.read_csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/negative_pool.csv')
	negatives=negatives['genes'].tolist()
	synapse_negatives=list(set(negatives)-set(synapse_genes))
	#random.shuffle(negatives)
	#synapse_negatives=negatives[:len(synapse_genes)]
	#print (len(synapse_negatives))

	return synapse_negatives

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def divide_list(gene_list, name):
	df=load_cell_locations_df()
	
	overlap=list(set(list(df.index))&set(gene_list))
	print ('overlap', len(overlap))
	overlap.sort()

	random.seed(4)
	random.shuffle(overlap)

	pool=overlap[:450]
	pool_df=pd.DataFrame({'genes': pool})
	pool_df.to_csv('regression_%s_pool.csv'%name)

	pool=pd.read_csv('regression_%s_pool.csv'%name)
	pool=pool['genes'].tolist()
	sublists=chunks(pool, 90)
	#print ('sublists', list(sublists))
	subs=[]
	for item in sublists:
		#print (len(item))
		sub=item
		subs.append(item)
	#training_pos=positives[:301]
	#held_out_pos=positives[301:]
	return subs

def find_model_input(df, pos, neg):
	cols=['Nucleus', 'Nucleoplasm', 'Nuclear bodies', 'Nuclear speckles', 'Nuclear membrane', 'Nucleoli', 'Nucleoli (Fibrillar center)', 'Cytosol', 'Cytoplasmic bodies', 'Lipid droplets', 'Mitochondria', 'Microtubules', 'Microtubule organizing center', 'Centrosome','Cytokinetic bridge', 'Midbody', 'Intermediate filaments', 'Actin filaments', 'Focal Adhesions', 'Endoplasmic reticulum', 'Golgi apparatus', 'Vesicles', 'Plasma membrane', 'Cell Junctions']

	df=df[cols]

	pos_df=df.loc[pos]
	pos_df['group']=1
	#print (pos_df)

	neg_df=df.loc[neg]
	neg_df['group']=0
	#print (neg_df)

	training=pd.concat([pos_df, neg_df])
	y=training['group'].tolist()
	training=training[cols]
	x=training.values
	return x, y


def find_training_heldout(pos_sublists, neg_sublists, i):

	#print (i)
	#print (len(pos_sublists))
	held_pos=pos_sublists[i]
	#print (len(held_pos))
	training_pos = [x for j,x in enumerate(pos_sublists) if j!=i] 
	training_pos = [item for sl in training_pos for item in sl]
	#training_pos=list(set(positives)-set(held_pos))
	#print (len(training_pos))

	held_neg=neg_sublists[i]
	training_neg = [x for j,x in enumerate(neg_sublists) if j!=i] 
	training_neg = [item for sl in training_neg for item in sl]

	training_x,training_y=find_model_input(df, training_pos, training_neg)
	held_x, held_y=find_model_input(df, held_pos, held_neg)
	return training_x, training_y, held_x, held_y

def plot_confusion_matrix(cm):
	fig, ax = plt.subplots(figsize=(8, 8))
	ax.imshow(cm)
	ax.grid(False)
	ax.xaxis.set(ticks=(0, 1), ticklabels=('Predicted 0s', 'Predicted 1s'))
	ax.yaxis.set(ticks=(0, 1), ticklabels=('Actual 0s', 'Actual 1s'))
	ax.set_ylim(1.5, -0.5)
	for i in range(2):
	    for j in range(2):
	        ax.text(j, i, cm[i, j], ha='center', va='center', color='red')
	plt.show()
	plt.close()

def find_model_predictions(pos_sublists, neg_sublists, locations):
	accuracies=[]
	frames=[]

	for i in range(5):
		training_x, training_y, held_x, held_y=find_training_heldout(pos_sublists, neg_sublists, i)

		model = LogisticRegression(solver='liblinear', random_state=0, penalty='l1').fit(training_x, training_y)
		model.predict_proba(held_x)
		accuracy=model.score(held_x, held_y)
		accuracies.append(accuracy)

		#locations=['Nucleus', 'Nucleoplasm', 'Nuclear bodies', 'Nuclear speckles', 'Nuclear membrane', 'Nucleoli', 'Nucleoli (Fibrillar center)', 'Cytosol', 'Cytoplasmic bodies', 'Lipid droplets', 'Mitochondria', 'Microtubules', 'Microtubule organizing center', 'Centrosome','Cytokinetic bridge', 'Midbody', 'Intermediate filaments', 'Actin filaments', 'Focal Adhesions', 'Endoplasmic reticulum', 'Golgi apparatus', 'Vesicles', 'Plasma membrane', 'Cell Junctions']
		thetaLasso=model.coef_
		coef=thetaLasso[0].tolist()
		#print (thetaLasso)
		reg_df=pd.DataFrame({'Locations': locations, 'Lasso': coef})
		reg_df=reg_df.set_index('Locations')
		frames.append(reg_df)
		cm = confusion_matrix(held_y, model.predict(held_x))
		#plot_confusion_matrix(cm)

	df=pd.concat(frames, axis=1)
	final=df
	mean=df.mean(axis=1)
	stderr=df.sem(axis=1)
	final['mean']=mean
	final['stderr']=stderr
	return accuracies, final


def find_solid_locations():
	df=load_cell_locations_df()
	df= df.iloc[:,:-5]
	df=df.drop(['ENSG', 'Uniprot'], axis=1)
	df=df.sum(axis=0)
	#df=df[df.columns[-5:]]
	print (df[df>30])
	locations=list(df[df>30].index)
	print (locations)
	return locations

if __name__ == '__main__':

	locations=find_solid_locations()

	df=load_cell_locations_df()
	
	positives=load_synapse_positives()
	
	pos_sublists=divide_list(positives, 'pos')

	negatives=load_synapse_negatives()
	neg_sublists=divide_list(negatives, 'neg')

	accuracies, final=find_model_predictions(pos_sublists, neg_sublists, locations)
		
	print (accuracies)
	print (np.mean(accuracies))

	#print (final)
	final=final.sort_values(by=['mean'])
	print (final)
	final.to_csv('regression_coef_cell_location.csv')
		