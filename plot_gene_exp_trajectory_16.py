import csv
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("TKAgg")
print(matplotlib.get_backend())
from matplotlib import pyplot as plt
from scipy import stats

def add_gene_names(df):
	geneid=df['Geneid'].tolist()
	symbols=[]
	for item in geneid:
		new=item[item.index('|')+1:]
		symbols.append(new)
	#print (symbols[:5])
	df.insert(loc=0, column='Genes', value=symbols)
	df=df.set_index('Genes')
	#print (df)
	return df

def get_ctx_df(df):
	df=add_gene_names(df)
	#print (df)

	cols=list(df.columns)
	ctx=['Geneid']
	for item in cols:
		if "DFC" in item:
			ctx.append(item)
	df=df[ctx]
	df=df.drop(['Geneid'], axis=1)
	#print (df)
	return df
#df=df.set_index('Geneid')

def key_dict():
	key_file=pd.read_csv('/Users/karenmei/Documents/BrainHierarchyDataSource/PsychEncode_key.csv', skiprows=3)
	#key_file.columns = key_file.iloc[2]
	key_file=key_file.drop(key_file.index[[0]])

	dfc_key=key_file[key_file['Regioncode']=='DFC']

	subjects=key_file['Braincode'].tolist()
	age=key_file['Days'].tolist()
	zipped_key=list(zip(subjects, age))
	key_dict=dict(zipped_key)

	#print (key_dict['HSB112'])
	return key_dict

def load_fetal():
	fetal=pd.read_csv('fetal_new_val.csv')
	fetal=fetal['genes'].tolist()
	return fetal

def load_adult():
	adult=pd.read_csv('adult_new_val.csv')
	adult=adult['genes'].tolist()
	return adult

def translate_df():

	df=pd.read_csv('/Users/karenmei/Documents/BrainHierarchyDataSource/PsychEncode_Ensembl_NoZero.csv')

	df=get_ctx_df(df)

	df.rename(columns=lambda x: x[:x.index('.DFC')], inplace=True)

	age_dict=key_dict()

	cols=list(df.columns)

	trans_cols=[]
	for item in cols:
		trans=age_dict[item]
		trans_cols.append(trans)

	df.columns=trans_cols
	df.sort_index(axis=1, level=None, ascending=True, inplace=True, kind='quicksort', na_position='last', sort_remaining=True, by=None)
	df.to_csv('DFC_expression.csv')
	return df

def find_overlap_df(gene_list):
	df=translate_df()
	#print (df)

	overlap=list(set(list(df.index))&set(gene_list))
	# print (len(overlap))

	fetal_df=df.loc[overlap]
	return fetal_df

def plot_age_mean(df, name):
	
	df.index.names = ['Age']
	#new=new.reset_index()
	cols=list(df.columns)

	mod_df=df

	mod_df['mean']=mod_df.mean(axis=1)
	mod_df['stderr']=df.sem(axis=1)
	print (df)
	new.to_csv('%s_trajectory.csv'%name)


	df=df.reset_index()
	plt.plot(df['Age'], df['mean'])
	#sns.lineplot(x='Age', y='mean', data=df, err_style='band')
	#plt.ylim([0, 25])
	plt.fill_between(df['Age'], df['mean']-df['stderr'], df['mean']+df['stderr'], alpha=0.5)
	plt.xscale('log')
	plt.show()

#import seaborn as sns; sns.set()
#import matplotlib.pyplot as plt
#fmri = sns.load_dataset("fmri")
#print (fmri)
#ax = sns.lineplot(x="timepoint", y="signal", data=fmri)
#plt.show()
from sklearn import preprocessing

df=translate_df()

fetal=load_fetal()
#fetal=['PARP1', 'SYT11', 'RAI14', 'PTPRO', 'SMG1', 'LMNA', 'REPS1', 'FN1', 'PRPF19', 'CIRBP', 'DCX', 'NUP93', 'HNRNPU', 'RBM14','TSC1', 'CALD1', 'CNOT1', 'DLG5', 'ILF3', 'HCFC1']
fetal_df=find_overlap_df(fetal)

print (fetal_df)
new=fetal_df.transpose()
new=new.groupby(new.index).mean()
#new=new.reset_index()
#new=new.groupby(np.arange(len(new))//3).mean()
new = new.groupby(by=new.columns, axis=1).mean()
#new=new.set_index('index')

plot_age_mean(new, "fetal")

#x = fetal_df.values #returns a numpy array
#min_max_scaler = preprocessing.MinMaxScaler()
#x_scaled = min_max_scaler.fit_transform(x)
#df = pd.DataFrame(x_scaled)
#print (df)

adult=['BRAF', 'CPNE9', 'MPRIP', 'WNK1', 'ABCE1', 'PRKAG2', 'SULT4A1', 'TIAM1', 'APOO', 'RAP1B', 'SYBU', 'TECR', 'UBB','APBA1', 'KIF1B', 'ANLN', 'PPP2CB', 'PTGES2', 'SRRM2', 'CLSTN3', 'MATK']
adult_df=find_overlap_df(adult)
new=adult_df.transpose()

new=new.groupby(new.index).mean()
new = new.groupby(by=new.columns, axis=1).mean()

#plot_age_mean(new)
#plt.show()
adult=load_adult()
adult_df=find_overlap_df(adult)
new=adult_df.transpose()
new=new.groupby(new.index).mean()
new = new.groupby(by=new.columns, axis=1).mean()
#adult_0=new.loc[7206.0].tolist()
#plt.hist(fetal_0, normed=True, alpha=0.5)
#plt.hist(adult_0, normed=True, alpha=0.5)
#plt.show()

#print (new)
plot_age_mean(new, "adult")

#t2, p2 = stats.ttest_ind(fetal_0,adult_0)
#print("t = " + str(t2))
#print("p = " + str(p2))