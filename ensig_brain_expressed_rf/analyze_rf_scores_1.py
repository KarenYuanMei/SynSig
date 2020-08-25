import numpy as np
from igraph import *
import pandas as pd
import sys
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

#-----load your ontology onto HiView-------------------------------------------------------------------------
#code for uploading to HiView taken from DDOT package: https://github.com/michaelkyu/ddot/blob/master/examples/Tutorial.ipynb


#make sure that you insert the path to your ddot folder
sys.path.append("/Users/karenmei/Documents/ddot/")

import ddot
from ddot import Ontology

import networkx as nx

import matplotlib
matplotlib.use("TKAgg")
from matplotlib import pyplot as plt

import seaborn as sns; sns.set()
import pandas as pd
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import explained_variance_score, mean_absolute_error, r2_score
from scipy.stats.stats import pearsonr, spearmanr
import pylab
from sklearn.datasets import make_regression

#construct the random forest so that when doing the 5X cross validation, the model is not seeing 20% of the genes, not just rows--------------------
import random
import pickle
#from define_gene_objects_rf_5 import PairOfGenes
import seaborn as sns
import matplotlib.patches as mpatches


def plot_all_predicted():
	frames=[]
	for i in range(5):

		data=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/nonbrain_features_brain_genes_pipeline/random_forest/ypredict_ytest_%s.csv'%i, index_col=[0])

		data=data.drop(['Gene1', 'Gene2'], axis=1)

		frames.append(data)

	df=pd.concat(frames)

	print (df)

	g=sns.jointplot(x=df.ytest, y=df.ypredict, kind='kde')
	
	y_test=df['ytest'].tolist()
	yfit=df['ypredict'].tolist()
	pearson_corr=pearsonr(y_test, yfit)[0]
	pearson_corr=np.round(pearson_corr,2)
	p_value=pearsonr(y_test, yfit)[1]

	r = pearson_corr
	p = p_value
	# if you choose to write your own legend, then you should adjust the properties then
	phantom, = g.ax_joint.plot([], [], linestyle="", alpha=0)
	g.ax_joint.legend([phantom],['r={:f}, p={:f}'.format(r,p)])
	#plt.title('RF Performance: Predicted vs. True',fontweight='bold')

	plt.show()
	plt.close()

#plot_all_predicted()

frames=[]
for i in range(5):
	data=pd.read_csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/nonbrain_features_brain_genes_pipeline/random_forest/ypredict_ytest_%s.csv'%i, index_col=[0])
	frames.append(data)

df=pd.concat(frames)

#print (df)

positive_filename='/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/brain_features_brain_genes_pipeline/synapse_positives.csv'
negative_filename='/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/brain_features_brain_genes_pipeline/synapse_negatives.csv'


def load_gene_lists(filename):
	df=pd.read_csv(filename)
	gene_list=df['genes'].tolist()
	return gene_list

positives=load_gene_lists(positive_filename)
#print (positives[:10])
negatives=load_gene_lists(negative_filename)

pos_pos_df=df[df["Gene2"].isin(positives)]

pos_pos_df['group']='positive'

print (pos_pos_df)



def plot_sample(df1, color1, df2, color2):
	# Sample 1000 random lines
	df1_sample=df1.sample(1000)
	df2_sample=df2.sample(1000)

	#df1_sample=df1
	#df2_sample=df2
	 
	# Make the plot with this subset
	plt.plot( 'ytest', 'ypredict', data=df1_sample, linestyle='', marker='o', color=color1)

	plt.plot( 'ytest', 'ypredict', data=df2_sample, linestyle='', marker='o', color=color2)
	plt.grid(False)

	red_patch = mpatches.Patch(color='red', label='Negative Genes')
	blue_patch=mpatches.Patch(color='blue', label='Synapse Genes')
	plt.legend(handles=[red_patch, blue_patch])
	 
	# titles
	plt.xlabel('True Semantic Similarity Scores')
	plt.ylabel('Predicted Semantic Similarity Scores')
	plt.title('Correlations of Predicted vs. True Scores', fontweight='bold')
	plt.show()
	 


neg_neg_df=df[df["Gene2"].isin(negatives)]


neg_neg_df['group']='negative'

print (neg_neg_df)

final=pd.concat([pos_pos_df, neg_neg_df], axis=0)

#print (final)
plot_sample(neg_neg_df, 'r', pos_pos_df, 'b')



def plot_df(df, color):
	g=sns.jointplot(x=df.ytest, y=df.ypredict, kind='kde', color=color, xlim=(-1,14), ylim=(-1, 11), joint_kws=dict(shade_lowest=False))
	
		
	y_test=df['ytest'].tolist()
	yfit=df['ypredict'].tolist()
	pearson_corr=pearsonr(y_test, yfit)[0]
	pearson_corr=np.round(pearson_corr,2)
	p_value=pearsonr(y_test, yfit)[1]

	r = pearson_corr
	p = p_value
	# if you choose to write your own legend, then you should adjust the properties then
	phantom, = g.ax_joint.plot([], [], linestyle="", alpha=0)
	g.ax_joint.legend([phantom],['r={:f}, p={:f}'.format(r,p)])
	#plt.title('RF Performance: Predicted vs. True',fontweight='bold')
	#plt.xlim(0, 14)
	#plt.ylim(0,14)
	plt.grid(False)


	#plt.show()
	#plt.close()

plot_df(pos_pos_df, 'b')

plot_df(neg_neg_df, 'r')
plt.show()
plt.close()
plot_df(final, 'grey')
plt.show()
plt.close()




# ax=sns.kdeplot(neg_neg_df.ytest, neg_neg_df.ypredict, cmap="Reds", shade=True, shade_lowest=False, alpha=0.6)
# ax = sns.kdeplot(pos_pos_df.ytest, pos_pos_df.ypredict, cmap="Greens", shade=True, shade_lowest=False, alpha=0.6)
# plt.title('RF Performance: Predicted vs. True',fontweight='bold')
# phantom, = ax.plot([], [], linestyle="", alpha=0)
# ax.legend([phantom],['synapse-synapse r=0.23, synapse-negative p=0.22'])
# plt.show()

#g = sns.FacetGrid(final, col="group", hue="group")
#g = (g.map(plt.scatter, "ytest", "ypredict", edgecolor="w"))

#plt.show()
 
