#Goal: Predicted vs. Reference semantic similarity scores correlation: for all of the validation sets



import sys
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

#-----load your ontology onto HiView-------------------------------------------------------------------------
#code for uploading to HiView taken from DDOT package: https://github.com/michaelkyu/ddot/blob/master/examples/Tutorial.ipynb

#make sure that you insert the path to your ddot folder
sys.path.append("/Users/karenmei/Documents/ddot/")

import matplotlib
matplotlib.use("TKAgg")
from matplotlib import pyplot as plt
matplotlib.rcParams.update({'font.size': 22})


import seaborn as sns;sns.set(style="white",font_scale=2)

sns.set_style("ticks")

from scipy.stats.stats import pearsonr, spearmanr
#from sklearn.datasets import make_regression

#construct the random forest so that when doing the 5X cross validation, the model is not seeing 20% of the genes, not just rows--------------------
import random


import numpy as np, pandas as pd; np.random.seed(0)
#sns.set(font_scale=2)
#import seaborn as sns; sns.set(style="white", color_codes=True)

def sum_df():
	frames=[]
	for i in range(5):

		data=pd.read_csv('ypredict_ytest_%s.csv'%i, index_col=[0])

		data=data.drop(['Gene1', 'Gene2'], axis=1)

		frames.append(data)

	df=pd.concat(frames)

	print (df)
	return df

def plot_all_predicted():
	df=sum_df()
	df.rename(columns={'ytest': 'Reference Scores', 'ypredict': 'Predicted Scores'}, inplace=True)
	
	#g=sns.jointplot(x=df['Predicted Scores'], y=df['Reference Scores'], kind='kde')
	g=sns.kdeplot(df['Predicted Scores'], df['Reference Scores'], cmap='Purples', shade=True, cbar=True)
	
	y_test=df['Reference Scores'].tolist()
	yfit=df['Predicted Scores'].tolist()
	pearson_corr=pearsonr(y_test, yfit)[0]
	pearson_corr=np.round(pearson_corr,3)
	p_value=pearsonr(y_test, yfit)[1]

	r = pearson_corr
	p = p_value

	print (r, p)

	print (spearmanr(y_test, yfit)[0])
	# if you choose to write your own legend, then you should adjust the properties then
	#phantom, = g.ax_joint.plot([], [], linestyle="", alpha=0)
	#g.ax_joint.legend([phantom],['r={:f}, p={:f}'.format(r,p)])
	#g.ax_joint.legend([phantom],['r={:f}, p<0.00001'.format(r,p)])
	#plt.title('RF Performance: Predicted vs. True',fontweight='bold')

	#plt.xlim(2, 11)
	plt.ylim(0, 14)

	plt.show()
	plt.close()

if __name__ == '__main__':
	plot_all_predicted()







