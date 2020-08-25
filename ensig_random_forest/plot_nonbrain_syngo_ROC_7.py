#Goal: 1) Draw ROC curve for predicted synapse gene scores against published experiments of synapse genes
#.     2) Find the synapse genes above a threshold

import numpy as np
import pandas as pd
import csv

import ddot
from ddot import Ontology

import matplotlib
matplotlib.use("TKAgg")
from matplotlib import pyplot as plt

from sklearn.metrics import auc

from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
plt.rcParams.update({'font.size': 22})
plt.rcParams["font.family"] = "Arial"

#Load the files with the avg predicted scores for each gene:
def load_predicted_df():
	df=pd.read_csv('nonbrain_RNA_big_pool_novel_synapse_genes_avg_scores.csv', index_col=[0])
	#print ('pred', df)
	return df

def load_training_genes():
	filename='/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/synapse_positives.csv'
	df=pd.read_csv(filename, index_col=[0])
	pos=df['genes'].tolist()

	filename='/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/synapse_negatives.csv'
	df=pd.read_csv(filename, index_col=[0])
	neg=df['genes'].tolist()

	genes=list(set(pos+neg))
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

def find_synsysnet():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Synapse_Genes/SynSysNet_genes.csv')
	#print (df)
	genes=df['gene_name'].tolist()
	#print (len(genes))
	training=load_training_genes()
	genes=list(set(genes)-set(training))
	return genes

def find_synDB():
	df=pd.read_csv('/Users/karenmei/Documents/BrainHierarchyDataSource/SynDB_Master.csv')
	#print (df)
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
	#print (df)
	genes=df['genes'].tolist()
	training=load_training_genes()
	genes=list(set(genes)-set(training))
	return genes

#find for each gene in the predicted datafile, how many databases also count it as a synapse gene
def find_union_overlap_count_pred(df, union):
	all_genes=df['genes'].tolist()
	freq=[]
	for item in all_genes:
		if item in union:
			count=1
		else:
			count=0
		freq.append(count)

	df['union']=freq

	#df.to_csv('all_genes_predicted_true_count_exp.csv')

	return df

def plot_ROC(df):
	probs=df['avg_scores'].tolist()
	y=df['union'].tolist()

	fpr, tpr, thresholds = roc_curve(y, probs)

	true_false=list(zip(tpr, fpr))

	ROC=list(zip(thresholds, tpr, fpr))

	ROC_df=pd.DataFrame({'Threshold': thresholds, "True_Positives": tpr, "False_Positives": fpr})

	ROC_df.to_csv('nonbrain_ROC_syngo.csv')

	#print (ROC)
	auc = roc_auc_score(y, probs)
	print('AUC: %.3f' % auc)

	#import matplotlib.pyplot as plt
	#plt.title('ROC Curve for Novel Synapse Predictions Using Experiments')
	plt.plot(fpr, tpr, 'g', label = 'AUC = %0.2f' %auc, linewidth=5)
	plt.legend(loc = 'lower right')
	plt.plot([0, 1], [0, 1],'--', color='grey')
	#plt.plot([0, 0.075], [0.395, 0.395],'--', color='blue')
	#plt.plot([0.075, 0.075], [0, 0.395],'--', color='blue')
	#print (fpr[603])
	#rnd_idx=603
	#plt.annotate('Synch\nThreshold', color='blue', xy=(fpr[rnd_idx], tpr[rnd_idx]), xytext=(fpr[rnd_idx]+0.05, tpr[rnd_idx]),
    #         arrowprops=dict(facecolor='blue', lw=2, arrowstyle='->'),)
	plt.xlim([0, 1])
	plt.ylim([0, 1])
	plt.ylabel('True Positive Rate', fontweight='bold')
	plt.xlabel('False Positive Rate', fontweight='bold')
	plt.show()


def find_synapse_above_threshold(df, threshold):
	#print ('df', df)
	df=df[df['avg_scores']>threshold]
	#print ('thresholded', df)

	novel_genes=df['genes'].tolist()
	df=pd.DataFrame({'genes': novel_genes})

	df.to_csv('nonbrain_pred_genes_above_%s.csv'%threshold)
	return novel_genes



#function that's the compiled version of the above functions
def find_frequency_val_df():
	df=load_predicted_df()
	pred_genes=df['genes'].tolist()
	new=pd.DataFrame({'Genes': pred_genes})
	new=find_pred_val_df(pred_genes, new)
	new['sum']=new.sum(axis=1)
	new.to_csv('nonbrain_pred_validation_matrix.csv')
	return new

#new=find_frequency_val_df()

def plot_lit_ROC():
	pred_val=pd.read_csv('nonbrain_pred_validation_matrix.csv', index_col=[0])
	#print (pred_val)

	fetal=load_fetal_brain()
	ngn2=load_ngn2()
	stria=load_adult_str()
	ctx=load_adult_ctx()

	#union=list(set(fetal)&set(ngn2)&set(stria)&set(ctx))

	syngo=load_syngo_genes()
	synsysnet=find_synsysnet()
	synDB=find_synDB()
	go_synapse=find_GO_synapse()

	#union=list(set(syngo)&set(synsysnet)&set(synDB)&set(go_synapse))

	#union=list(set(syngo)&set(fetal+ngn2)&set(stria+ctx))

	#union=list(set(fetal+ngn2)&set(stria+ctx))

	#union=list(set(syngo+fetal+ngn2+stria+ctx))

	#union=list(set(syngo)&set(fetal)&set(ngn2)&set(stria)&set(ctx))
	union=list(set(syngo))
	print ('union', len(union))

	pred=load_predicted_df()
	#print (pred)
	union_df=find_union_overlap_count_pred(pred, union)

	#print (union_df)
	plot_ROC(union_df)

def find_network_above_threshold(threshold):
	pred=load_predicted_df()

	novel_genes=find_synapse_above_threshold(pred, threshold)
	print ('novel', len(novel_genes))

	#pred_val_final=pred_val.set_index('Genes')
	#print (novel_genes)

	#thres_sum=pred_val_final.loc[novel_genes]

	#print (thres_sum[thres_sum['sum']>0])
	return novel_genes

if __name__ == '__main__':
	#find_frequency_val_df()
	plot_lit_ROC()
	find_network_above_threshold(4.67)

