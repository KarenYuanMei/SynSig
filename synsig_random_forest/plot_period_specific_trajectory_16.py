#Goal: 
#overlap: the intersection of syndromic autism, fetal experiments, and nonbrain predicted synapse
#Plot the trajectory of overlap genes
import csv
import pandas as pd
import numpy as np

import matplotlib
#matplotlib.use("TKAgg")
#print(matplotlib.get_backend())
from matplotlib import pyplot as plt
from scipy import stats
from sklearn import preprocessing
import random


def add_gene_names(df):
	geneid=df['Geneid'].tolist()
	symbols=[]
	for item in geneid:
		new=item[item.index('|')+1:]
		symbols.append(new)
	df.insert(loc=0, column='Genes', value=symbols)
	df=df.set_index('Genes')
	return df

def get_ctx_df(df):
	df=add_gene_names(df)
	cols=list(df.columns)
	ctx=['Geneid']
	for item in cols:
		if "DFC" in item:
			ctx.append(item)
	df=df[ctx]
	df=df.drop(['Geneid'], axis=1)
	return df

# def make_github_source():
# 	df=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/Data_Source/PsychEncode_Ensembl_NoZero.csv')
# 	print (df)
# 	df=add_gene_names(df)
# 	print (df)
# 	df=get_ctx_df(df)
# 	print (df)

# 	fetal=pd.read_csv('fetal_specific.csv')
# 	fetal=fetal['genes'].tolist()
# 	adult=pd.read_csv('adult_specific.csv')
# 	adult=adult['genes'].tolist()

# 	specific_genes=list(set(fetal+adult))

# 	df=df.loc[specific_genes]

# 	df.rename(columns=lambda x: x[:x.index('.DFC')], inplace=True)
# 	df.to_csv('PsychEncode_Ensembl_DFC.csv')
# 	return df

def key_dict():
	key_file=pd.read_csv('../other_resources/PsychEncode_key.csv', skiprows=3)
	key_file=key_file.drop(key_file.index[[0]])

	dfc_key=key_file[key_file['Regioncode']=='DFC']

	subjects=key_file['Braincode'].tolist()
	age=key_file['Days'].tolist()
	zipped_key=list(zip(subjects, age))
	key_dict=dict(zipped_key)
	return key_dict

def translate_df():

	df=pd.read_csv('../other_resources/PsychEncode_Ensembl_DFC.csv')
	df=df.set_index('Genes')
	print (df)

	age_dict=key_dict()

	cols=list(df.columns)

	trans_cols=[]
	for item in cols:
		trans=age_dict[item]
		trans_cols.append(trans)

	df.columns=trans_cols
	df.sort_index(axis=1, level=None, ascending=True, inplace=True, kind='quicksort', na_position='last', sort_remaining=True, by=None)
	#df.to_csv('DFC_expression.csv')
	return df

def find_overlap_df(gene_list):
	df=translate_df()
	overlap=list(set(list(df.index))&set(gene_list))
	fetal_df=df.loc[overlap]
	return fetal_df

def plot_age_mean(df, name):
	
	df.index.names = ['Age']
	#new=new.reset_index()
	cols=list(df.columns)

	mod_df=df

	mod_df['mean']=mod_df.mean(axis=1)
	mod_df['stderr']=df.sem(axis=1)
	#print (df)
	df.to_csv('%s_trajectory.csv'%name)

	df=df.reset_index()
	f = plt.figure()
	plt.plot(df['Age'], df['mean'])
	#sns.lineplot(x='Age', y='mean', data=df, err_style='band')
	plt.ylim([0, 70])
	plt.fill_between(df['Age'], df['mean']-df['stderr'], df['mean']+df['stderr'], alpha=0.5)
	plt.xscale('log')
	plt.show()
	f.savefig("%s_trajectory.pdf"%name, bbox_inches='tight')


def find_fetal_trajectory():
	df=translate_df()
	overlap=pd.read_csv('../other_resources/fetal_specific.csv')
	overlap=overlap['genes'].tolist()
	#overlap=['SUPT16H', 'XRCC5', 'SSRP1', 'XRCC6', 'TMOD3', 'FAT3', 'KHDRBS1', 'PCDH10', 'RAI14', 'ILF2', 'SLC25A1', 'PRMT1', 'OCIAD1', 'CIRBP', 'UBQLN2', 'G3BP1', 'UBQLN1', 'APEX1', 'SLC25A6', 'TUBG1', 'SNRNP200', 'RUNDC3A', 'SNRPN', 'SRSF1', 'TMPO', 'DHX15', 'PRPF19', 'PRPF8', 'AGO1']
	synd_df=find_overlap_df(overlap)
	#print (synd_df)
	new=synd_df.transpose()
	new=new.groupby(new.index).mean()
	new = new.groupby(by=new.columns, axis=1).mean()
	plot_age_mean(new, "fetal_specific")

def find_adult_trajectory():
	
	df=translate_df()
	overlap=pd.read_csv('../other_resources/adult_specific.csv')
	overlap=overlap['genes'].tolist()

	#overlap=['ANLN', 'OTUD7A', 'REEP2', 'MATK', 'TCEAL5', 'MYL12B', 'PREPL', 'CPNE9', 'JPH3', 'SLITRK5']
	print ('adult', len(overlap))
	#df=pd.DataFrame({'genes': overlap})
	#df.to_csv('adult_specific.csv')
	synd_df=find_overlap_df(overlap)
	#print (synd_df)
	new=synd_df.transpose()
	new=new.groupby(new.index).mean()
	new = new.groupby(by=new.columns, axis=1).mean()
	plot_age_mean(new, "adult_specific")

if __name__ == '__main__':
	find_fetal_trajectory()
	find_adult_trajectory()
