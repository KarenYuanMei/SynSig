#Goal: compare synsig genes vs. non-synsig genes in cell lines

import csv

import numpy as np

from mlxtend.evaluate import permutation_test
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy import stats

import pandas as pd

import matplotlib
matplotlib.use("TKAgg")
from matplotlib import pyplot as plt

import random
plt.style.use('seaborn-deep')


#find all of the different cell types
def find_cell_type():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/cell_line_expression/Thul_2017/cell_line_desc.csv')
	null_columns=df.columns[df.isnull().any()]
	types_df=df[df["Description of cell line"].isnull()][null_columns]
	types=types_df['Name'].tolist()
	types = [x for x in types if str(x) != 'nan']
	#print (types)
	types=list(set(types))
	types.sort()
	print ('number of types', len(types))
	return types

#find_cell_type()

def find_lines_per_type():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/cell_line_expression/Thul_2017/cell_line_desc.csv')
	#df=df.set_index('Name')

	idx_list=[]
	types=find_cell_type()
	for item in types:
		print (item)
		idx=df.index[df['Name']==item].tolist()
		idx_list.append(idx)
	print (idx_list)
	return idx_list
#find_lines_per_type()

def find_nonbrain_cells():
	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/cell_line_expression/Thul_2017/Thul_2017.csv')
	cell_lines=list(df.columns)

	#remove "ENSG", 'Gene'
	lines=cell_lines[2:]
	line_names=[]
	for item in lines:
		new=item[:item.index(' (')]
		line_names.append(new)
	desc=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/cell_line_expression/Thul_2017/cell_line_desc.csv')
	desc=desc.set_index('Name')
	brain_lines=['SH-SY5Y', 'U-251 MG', 'U-138 MG', 'AF22', 'U-87 MG']
	nonbrain_lines=list(set(line_names)-set(brain_lines))
	new_df=desc.loc[nonbrain_lines]
	result=new_df[new_df['Origin']!= 'Brain']
	nonbrain_cells=list(result.index)
	df=pd.DataFrame({'Cell_Lines': nonbrain_cells})
	df.to_csv('NonBrain_Cells.csv')
	return nonbrain_cells

#load the positive synapse genes (training)

def load_synapse_positives():
	synapse_genes=pd.read_csv('synapse_positives.csv')
	synapse_genes=synapse_genes['genes'].tolist()
	return synapse_genes

#load the negative synapse genes (training)
def load_synapse_negatives():
	synapse_negatives=pd.read_csv('synapse_negatives.csv')
	synapse_negatives=synapse_negatives['genes'].tolist()
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

def find_cell_pos_neg(positives, negatives, data):
	
	positives=list(set(positives)&set(list(data.index)))
	print (len(positives))

	negatives=list(set(negatives)&set(list(data.index)))
	print (len(negatives))
	return positives, negatives

def normalize_cell_df():
	data=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/cell_line_expression/Thul_2017/Thul_2017.csv')

	cols=list(data.columns)

	new_headers=[]
	for item in cols:
		if "(TPM)" in item:
			new=item[:item.index(' (TPM)')]
		else:
			new=item
		new_headers.append(new)

	data.columns = new_headers
	data=data.set_index('Gene')
	#print (data)
	return data

#plot the distributions of the synapse vs. non-synapse proteins
def plot_distributions(positives, negatives, tissue, name):

	bins=np.histogram(np.hstack((positives,negatives)), bins=40)[1] #get the bin edges
	
	plt.hist(positives, bins, alpha=0.5, edgecolor='black', linewidth=0.6)
	plt.hist(negatives, bins, alpha=0.5, edgecolor='black', linewidth=0.6)
	#plt.grid(False)

	plt.xlabel('%s mRNA Expression'%tissue, fontweight='bold')
	plt.ylabel('Frequency', fontweight = 'bold')
	#plt.xscale('log')
	plt.yscale('log')
	plt.title('Expression for Synapse Genes. vs. Non-Synapse Genes', fontweight = 'bold')
	#plt.title('%s Expression for Synapse Genes.'%tissue, fontweight = 'bold')
	plt.legend(labels=['Mean Synapse Genes is %s'%round(np.mean(positives), 2), 'Mean Non-Synapse Genes is %s'%round(np.mean(negatives), 2)])
	
	#plt.show()
	#plt.close()
	plt.savefig('%s.pdf'%str(name), bbox_inches='tight')
	plt.close()

def correct_fdr(df, col_name):
	col_list=np.array(df[col_name].tolist())
	corr_list=multipletests(col_list, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
	return corr_list

def analyze_cell_lines(data, positives, negatives, nonbrain_cells, group):
	fc_list=[]

	sem_list=[]

	p_list=[]
	for item in nonbrain_cells:
		nonbrain_df=data[item]
		#print (nonbrain_df)
		pos=nonbrain_df.loc[positives]
		print (pos)
		pos_list=pos.tolist()
		print (pos_list[:5])
		print ('pos', len(pos_list), np.mean(pos_list))

		neg=nonbrain_df.loc[negatives]
		neg_list=neg.tolist()
		#neg_list=[x for x in neg_list if x != 'nan']
		print ('neg', len(neg_list), np.mean(neg_list))
		#print (neg_list)

		p_value=permutation_test(pos_list, neg_list, method='approximate', num_rounds=10000, seed=0)
		print ('permutation test p', p_value)

		name=item[:5]

		fc=np.mean(pos_list)/np.mean(neg_list)
		fc_list.append(fc)

		norm_pos=pos_list/np.mean(neg_list)
		sem=stats.sem(norm_pos)
		sem_list.append(sem)

		p_list.append(p_value)

		print (len(nonbrain_cells), len(fc_list), len(sem_list), len(p_list))

	df=pd.DataFrame({'Cell_Line': nonbrain_cells, 'Fold_Change': fc_list, 'SEM': sem_list, 'Significance': p_list})
	print (df)
	return df

def make_analysis_df(positives, negatives, group):
	data=normalize_cell_df()
	pos, neg=find_cell_pos_neg(positives, negatives, data)
	nonbrain_cells=find_nonbrain_cells()
	df=analyze_cell_lines(data, pos, neg, nonbrain_cells, group)
	corr_df=correct_fdr(df, 'Significance')
	print (corr_df)
	df['Corr_Significance']=corr_df[1]
	print (df)
	df.to_csv('%s_Cell_Lines_Exp.csv'%group)
	return df

		#plot_distributions(pos_list, neg_list, item, name)
#data=normalize_cell_df()
if __name__ == '__main__':

	positives=load_synapse_positives()
	negatives=load_synapse_negatives()

	df1=make_analysis_df(positives, negatives, 'Training')

	positives=load_expanded_positives()
	print (positives[:5])
	negatives=load_expanded_negatives()
	df2=make_analysis_df(positives, negatives, 'SynSig')
	
	df1['Group']='Training'

	df2['Group']='SynSig'
	final=pd.concat([df1, df2])
	final.to_csv('final_cells_fc.csv')

	print (final)



