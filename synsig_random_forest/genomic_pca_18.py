#Goal: to make the dataframe of genomic and protein structural features for input for PCA (which is analyzed in R: see plot_R_figures/genomic_PCA.R)


import pandas as pd
import numpy as np
#import matplotlib
#matplotlib.use("TKAgg")
#from matplotlib import pyplot as plt

import csv
#import random

#from sklearn.preprocessing import StandardScaler
#from sklearn.decomposition import PCA
import seaborn as sns

def load_synapse_positives():
	#synapse_genes=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_8/input_genes/synapse_positives.csv')
	synapse_genes=pd.read_csv('synapse_positives.csv')
	synapse_genes=synapse_genes['genes'].tolist()
	return synapse_genes

def load_synapse_negatives():
	#synapse_negatives=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_8/input_genes/synapse_negatives.csv')
	synapse_genes=pd.read_csv('synapse_negatives.csv')
	synapse_negatives=synapse_genes['genes'].tolist()
	return synapse_negatives

def load_expanded_positives():
	#novel_synapse_genes=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_8/synapse_genes_above_5.3.csv', usecols=['genes'])
	novel_synapse_genes=pd.read_csv('pred_genes_above_4.7.csv', usecols=['genes'])
	novel_synapse_genes=novel_synapse_genes['genes'].tolist()
	return novel_synapse_genes

def load_expanded_negatives():
	#negative_genes=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Analyze_Synapse_Features/expanded_negative_list.csv')
	negative_genes=pd.read_csv('synsig_negative_list.csv')
	negatives=negative_genes['Genes'].tolist()
	return negatives

def load_data(data_name):
	df=pd.read_csv('../features/normalized_%s.csv'%data_name)
	df=df.set_index('Genes')
	return df

def load_genomic_features():
	#data_name_list=['gc_content','Phosphosite_hu_no', 'pFam_domain_number', 'protein_mass', 'trans_len', 'gene_length', 'exon_no', 'cds_length', 'Ensembl_isoform_no', 'Ensembl_aa_length']
	data_name_list=['Phosphosite_hu_no', 'pFam_domain_number', 'protein_mass', 'trans_len', 'gene_length', 'exon_no', 'cds_length', 'Ensembl_isoform_no', 'Ensembl_aa_length']

	frames=[]
	for item in data_name_list:
		df=load_data(item)
		df=df.set_index('Norm_Symbol')
		frames.append(df)

	new=pd.concat(frames, axis=1, join='inner')
	return new

def find_gene_df(df, gene_list, name):
	new=df.loc[gene_list]
	new['target']=name
	new=new.dropna()
	return new

def make_exp_df():
	new=load_genomic_features()
	new_idx=list(new.index)

	print (new)
	expanded_positives=load_expanded_positives()
	exp_positives=find_gene_df(new, expanded_positives, 'Expanded Positives')
	print (exp_positives)
	print (len(list(set(list(exp_positives.index)))))

	expanded_negatives=load_expanded_negatives()
	exp_negatives=find_gene_df(new, expanded_negatives, 'Expanded Negatives')
	print (exp_negatives)
	print (len(list(set(list(exp_negatives.index)))))

	print (set(list(exp_positives.index))&set(list(exp_negatives.index)))

	final=pd.concat([exp_positives, exp_negatives], axis=0)
	print ('final', final)
	print ('index', list(final.index)[:5], len(list(final.index)), len(list(set(final.index))))
	final.columns=['Phos Site No', 'Domain No', 'Protein Mass', 'Transcript Length', 'Gene Length', 'Exon No', 'CDS Length', 'Isoform No', 'Amino Acid Length', 'target']
	final.to_csv('expanded_pca_R.csv')
	return final

if __name__ == '__main__':
	make_exp_df()
