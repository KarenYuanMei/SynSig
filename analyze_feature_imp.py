import numpy as np
import pandas as pd
import csv
from scipy.stats.stats import pearsonr
from collections import defaultdict
from itertools import combinations, combinations_with_replacement
from itertools import product
from scipy import spatial

import networkx as nx
import os
import ddot
from ddot import Ontology

import pickle

import matplotlib
matplotlib.use("TKAgg")
from matplotlib import pyplot as plt

#import seaborn as sns; sns.set()
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import explained_variance_score, mean_absolute_error, r2_score
from scipy.stats.stats import pearsonr, spearmanr
import pylab
from sklearn.datasets import make_regression


frames=[]
for i in range(5):
	df=pd.read_csv('Feature_Importance_%s.csv'%i, index_col=[0])
	df=df.set_index('Features')
	frames.append(df)

#print (frames)
new=pd.concat(frames, axis=1)
#print (new)

new['mean']=new.mean(axis=1)
#print (new)

new=new.sort_values(by='mean')
#print (new)

new.loc['sum']=new.sum()
#print(new)

def find_subdf(cat_list):
	exp=[]
	features=list(new.index)
	for item in cat_list:
		for feature in features:
			if item in feature:
				exp.append(feature)

	subdf=new.loc[exp]

	#print (subdf)

	subdf.loc['sum']=subdf.sum()
	print (subdf)
	return subdf

exp_list=['exp', 'RNA']
df=find_subdf(exp_list)

gene_list=['chr', 'cds_length', 'trans', 'isoform_no', 'exon','gc_content', 'gene_length']
df=find_subdf(gene_list)

protein_list=['protein', 'domain', 'aa_length']
df=find_subdf(protein_list)

dynamics=['mentha', 'phos', 'Phos']
df=find_subdf(dynamics)

brain=['RNA']
df=find_subdf(brain)

nb=['isoform_exp']
df=find_subdf(nb)

feature_imp=[0.132, 0.486, 0.144, 0.066, 0.12]
feature_cat=['Brain Expression', 'Non-Brain Expression', 'Genomic Sequence', 'Protein Architecture', 'Interaction Dynamics']

df=pd.DataFrame({'Features': feature_cat, 'Importance': feature_imp})
print (df)

ax = df.plot.bar(x='Features', y='Importance', rot=0, color='purple')
ax.grid(b=False)
plt.show()
plt.close()