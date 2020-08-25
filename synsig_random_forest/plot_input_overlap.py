#Goal: to plot the overlap of the positive examples; fig 2a

import pandas as pd
import networkx as nx
import numpy as np
import os


#import ddot
#from ddot import Ontology

#os.environ['KMP_DUPLICATE_LIB_OK']='True'

import matplotlib
matplotlib.use("TKAgg")
print(matplotlib.get_backend())
from matplotlib import pyplot as plt

from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles


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

def load_gene_pool_index():
	df=pd.read_csv('big_pool_genes_index.csv')
	genes=df['genes'].tolist()
	return genes


gene_pool=load_gene_pool_index()
synsysnet=find_synsysnet()
synDB=find_synDB()
GO_synapse=find_GO_synapse()

synsysnet_pool=list(set(gene_pool)&set(synsysnet))
synDB_pool=list(set(gene_pool)&set(synDB))
GO_synapse_pool=list(set(gene_pool)&set(GO_synapse))

v=venn3([set(synsysnet_pool), set(synDB_pool), set(GO_synapse_pool)], set_labels = ('SynSysNet', "SynaptomeDB", 'GO Synapse'))

for text in v.set_labels:
    text.set_fontweight('bold')
for text in v.set_labels:
    text.set_fontsize(25)
for text in v.subset_labels:
    text.set_fontsize(25)

v.get_patch_by_id('111').set_color('orange')
v.subset_labels[-1].set_fontweight('bold')
v.subset_labels[-1].set_fontsize(30)
plt.show()
plt.close()




# #SynSysNet dataset for all synapse genes:------------------------------------------------------------------------------
# SynSysNet=load_data('/Users/karenmei/Documents/Synapse_Ontology/Synapse_Genes/SynSysNet_genes.csv')

# SynSysNet_genes=[]
# for i in range(1, len(SynSysNet)):
# 	entry=SynSysNet[i][0]
# 	SynSysNet_genes.append(entry)
# #print (SynSysNet_genes[:5])

# #-----------------------------------------------------------------------------------------------
# SynDB=load_data('/Users/karenmei/Documents/BrainHierarchyDataSource/SynDB_Master.csv')

# SynDB_genes=[]
# for i in range(1, len(SynDB)):
# 	entry=SynDB[i][1]
# 	SynDB_genes.append(entry)

# Go_synapse_genes=load_data('/Users/karenmei/Documents/Synapse_Ontology/Synapse_Scripts/GO_Synapse.csv')