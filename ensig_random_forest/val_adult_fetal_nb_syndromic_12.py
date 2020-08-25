#Goal: find the which non-brain predicted genes that are syndromic autism are also found in the experiments; Fig7d

import pandas as pd
#import networkx as nx
import numpy as np

import csv
import ddot
from ddot import Ontology

import matplotlib
matplotlib.use("TKAgg")
print(matplotlib.get_backend())
from matplotlib import pyplot as plt

from matplotlib_venn import venn3, venn3_unweighted
from matplotlib_venn import venn2, venn2_circles


#overlap with nb and syndromic 


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

def find_fetal():
	fetal_brain=load_fetal_brain()
	ngn2=load_ngn2()
	#fetal=list(set(ngn2)&set(fetal_brain))
	fetal=list(set(ngn2+fetal_brain))
	return fetal

def find_adult():
	ctx=load_adult_ctx()
	stria=load_adult_str()
	#adult=list(set(ctx)&set(stria))
	adult=list(set(ctx+stria))
	return adult

def find_synsysnet():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Synapse_Genes/SynSysNet_genes.csv')
	#print (df)
	genes=df['gene_name'].tolist()
	#print (len(genes))
	return genes

def find_synDB():
	df=pd.read_csv('/Users/karenmei/Documents/BrainHierarchyDataSource/SynDB_Master.csv')
	#print (df)
	genes=df['Symbol'].tolist()
	return genes

def find_GO_synapse():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Synapse_Scripts/GO_Synapse.csv')
	#print (df)
	genes=df['genes'].tolist()
	return genes

def load_syngo_genes():
	syngo=Ontology.from_table('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Metrics/SynGO_BP.txt')
	syngo_bp_genes=syngo.genes
	syngo=Ontology.from_table('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Metrics/SynGO_CC.txt')
	syngo_cc_genes=syngo.genes
	syngo_genes=list(set(syngo_bp_genes+syngo_cc_genes))
	return syngo_genes

def load_nonbrain_pred_genes():
	pred_file='nonbrain_pred_genes_above_4.67.csv'
	pred_df=pd.read_csv(pred_file, index_col=[0])
	pred_genes=pred_df['genes'].tolist()
	pred_genes=[x.upper() for x in pred_genes]
	return pred_genes

def find_sfari_syndromic_genes():
	#df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Disease_genes/Autism/SFARI_genes_051120.csv')
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Disease_genes/Autism/SFARI-Gene_genes.csv')
	
	syn = df[df['genetic-category'].str.contains('Syndromic', regex=False, case=False, na=False)]
	#print (syn)
	genes=syn['gene-symbol'].tolist()

	df=df[df['syndromic']==1]
	genes2=df['gene-symbol'].tolist()

	#print (len(genes), len(genes2), len(list(set(genes)&set(genes2))))

	genes=list(set(genes+genes2))
	return genes

def find_db_genes():
	synsysnet=find_synsysnet()
	synDB=find_synDB()
	GO_synapse=find_GO_synapse()
	syngo=load_syngo_genes()

	db=list(set(synsysnet+synDB+GO_synapse+syngo))
	return db


db=find_db_genes()
nb=load_nonbrain_pred_genes()
syndromic=find_sfari_syndromic_genes()
print ('syndromic', len(syndromic))

overlap=list(set(nb)&set(syndromic))

overlap_no_db=list(set(overlap)-set(db))

#overlap=['DDX3X', 'CTNNB1', 'WAC', 'ASH1L', 'ASPM', 'ARID1B', 'PCCA', 'KIF14', 'CNKSR2', 'CHD1', 'CHD3', 'KDM5B', 'TSC2', 'KCNQ2', 'MED13L', 'RAC1', 'KMT2A', 'CHD8', 'AFF4', 'CHD2', 'ANKRD11', 'SMARCA4', 'ZBTB20', 'TRIP12', 'MYT1L', 'GRID2', 'MTOR', 'SYNE1', 'EP300', 'DYRK1A', 'DPYD', 'SCN2A', 'ARNT2', 'DYNC1H1', 'TCF4', 'RALA', 'SCN8A', 'SLC6A1', 'MED13', 'ATRX', 'KMT2C', 'EHMT1', 'HNRNPU', 'PRKDC', 'SETD2', 'GRIN2B', 'PHIP', 'RAI1', 'NTRK3', 'TCF20']

fetal=find_fetal()
adult=find_adult()

fetal_overlap=list(set(overlap_no_db)&set(fetal))
#fetal_overlap=list(set(fetal_overlap)-set(db))
print (len(fetal_overlap))
print (fetal_overlap)

adult_overlap=list(set(overlap_no_db)&set(adult))
#adult_overlap=list(set(adult_overlap)-set(db))
print (len(adult_overlap))
print (adult_overlap)

fetal_df=pd.DataFrame({'Stage': 'Fetal Synapse: Syndromic Autism Genes','Genes':list(set(fetal_overlap))})
adult_df=pd.DataFrame({'Stage': 'Adult Synapse: Syndromic Autism Genes','Genes':list(set(adult_overlap))})

final=pd.concat([fetal_df, adult_df], axis=0)
print (final)
final.to_csv('nb_val_new_syndromic.csv')

#print (df)

pred_val=list(set(fetal+adult)&set(nb))

v=venn3_unweighted([set(pred_val),set(syndromic), set(db)], set_labels=('Proteomics Validated \n ENSig',  'Syndromic Autism', 'Synapse Databases'), set_colors=('skyblue', 'coral', 'gray'),alpha=0.7)
#venn3_circles([set(pred_val),set(syndromic), set(db)], linestyle='solid', linewidth=0.5, color='k');
for text in v.set_labels:
	#print (text)
	text.set_fontweight('bold')
for text in v.set_labels:
    text.set_fontsize(30)
for text in v.subset_labels:
	print (text)
	text.set_fontsize(30)

target=v.subset_labels[2]
target.set_fontweight('bold')
target.set_fontsize(35)
v.get_patch_by_id('110').set_color('red')
plt.show()
plt.close()
