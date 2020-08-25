import pandas as pd
#import networkx as nx
import numpy as np
#import os
import ddot
from ddot import Ontology
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

from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles
import venn

from sklearn.metrics import auc

from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
plt.rcParams.update({'font.size': 22})
plt.rcParams["font.family"] = "Arial"

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

def load_fetal_brain():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Validation_proteomics/Coba_human_fetal_2020/coba_fetal_brain.csv')
	print (df)
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
	
def load_training_genes():
	filename='/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/brain_RNA_big_gene_pool_pipeline/synapse_positives.csv'
	df=pd.read_csv(filename, index_col=[0])
	genes=df['genes'].tolist()
	return genes


def load_syngo_genes():
	syngo=Ontology.from_table('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Metrics/SynGO_BP.txt')
	syngo_bp_genes=syngo.genes
	syngo=Ontology.from_table('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Metrics/SynGO_CC.txt')
	syngo_cc_genes=syngo.genes
	syngo_genes=list(set(syngo_bp_genes+syngo_cc_genes))
	return syngo_genes


def find_synsysnet():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Synapse_Genes/SynSysNet_genes.csv')
	#print (df)
	genes=df['gene_name'].tolist()
	print (len(genes))
	return genes

def find_synDB():
	df=pd.read_csv('/Users/karenmei/Documents/BrainHierarchyDataSource/SynDB_Master.csv')
	print (df)
	genes=df['Symbol'].tolist()
	return genes

def find_GO_synapse():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Synapse_Scripts/GO_Synapse.csv')
	print (df)
	genes=df['genes'].tolist()
	return genes


#Load the files with the avg predicted scores for each gene:
def load_predicted_df():
	#df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/brain_RNA_big_gene_pool_pipeline/brain_RNA_big_pool_novel_synapse_genes_avg_scores.csv', index_col=[0])
	#print ('pred', df)
	pred_file='/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/nonbrain_big_genes_pool_pipeline/nonbrain_RNA_big_pool_novel_synapse_genes_avg_scores.csv'
	df=pd.read_csv(pred_file)
	return df

def find_true_y():
	ctx=load_adult_ctx()
	stria=load_adult_str()
	fetal=load_fetal_brain()
	ngn2=load_ngn2()

	adult=list(set(ctx)&set(stria))
	fetal=list(set(fetal)&set(ngn2))


	syngo=load_syngo_genes()
	synsysnet=find_synsysnet()
	synDB=find_synDB()
	go_synapse=find_GO_synapse()
	db=list(set(syngo+synsysnet+synDB+go_synapse))
	db=list(set(syngo)&set(synsysnet)&set(synDB)&set(go_synapse))

	full=list(set(adult)&set(fetal)&set(syngo)&set(synsysnet)&set(synDB)&set(go_synapse))
	full=list(set(ctx)&set(stria)&set(fetal)&set(ngn2)&set(synsysnet))
	full=list(set(adult)&set(fetal)&set(db))
	print (len(full))
	

	df=load_predicted_df()
	avg_scores=df['avg_scores'].tolist()
	pred_genes=df['genes'].tolist()

	y_list=[]
	for item in pred_genes:
		if item in full:
			group=1
		else:
			group=0
		y_list.append(group)

	final=pd.DataFrame({'genes': pred_genes, 'avg_scores': avg_scores , 'union': y_list})
	return final


def plot_ROC(df):
	probs=df['avg_scores'].tolist()
	y=df['union'].tolist()

	fpr, tpr, thresholds = roc_curve(y, probs)

	true_false=list(zip(tpr, fpr))

	ROC=list(zip(thresholds, tpr, fpr))

	ROC_df=pd.DataFrame({'Threshold': thresholds, "True_Positives": tpr, "False_Positives": fpr})

	#ROC_df.to_csv('ROC_on_weijun.csv')

	#print (ROC)
	auc = roc_auc_score(y, probs)
	print('AUC: %.3f' % auc)

	#import matplotlib.pyplot as plt
	#plt.title('ROC Curve for Novel Synapse Predictions Using Experiments')
	plt.plot(fpr, tpr, 'g', label = 'AUC = %0.2f' %auc, linewidth=5)
	plt.legend(loc = 'lower right')
	plt.plot([0, 1], [0, 1],'--', color='grey')
	plt.plot([0, 0.0902], [0.4015, 0.4015],'--', color='darkgray')
	plt.plot([0.0902, 0.0902], [0, 0.4015],'--', color='darkgray')
	#print (fpr[517], tpr[517])
	#rnd_idx=517
	#plt.annotate('SynSig\nThreshold', color='black', xy=(fpr[rnd_idx], tpr[rnd_idx]), xytext=(fpr[rnd_idx]+0.05, tpr[rnd_idx]),
             #arrowprops=dict(facecolor='darkgray', lw=2, arrowstyle='->'),)
	plt.xlim([0, 1])
	plt.ylim([0, 1])
	plt.ylabel('True Positive Rate', fontweight='bold')
	plt.xlabel('False Positive Rate', fontweight='bold')
	plt.show()


def plot_adult_ROC():
	
	final=find_true_y()

	#print (union_df)
	plot_ROC(final)


if __name__ == '__main__':
	plot_adult_ROC()
	
