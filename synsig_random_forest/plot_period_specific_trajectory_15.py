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

def key_dict():
	key_file=pd.read_csv('/../other_resources/PsychEncode_key.csv', skiprows=3)
	key_file=key_file.drop(key_file.index[[0]])

	dfc_key=key_file[key_file['Regioncode']=='DFC']

	subjects=key_file['Braincode'].tolist()
	age=key_file['Days'].tolist()
	zipped_key=list(zip(subjects, age))
	key_dict=dict(zipped_key)
	return key_dict

def find_regions():
	key_file=pd.read_csv('/../other_resources/PsychEncode_key.csv', skiprows=3)
	key_file=key_file.drop(key_file.index[[0]])

	regions=key_file['Regioncode'].tolist()
	regions=list(set(regions))
	return regions

#print (find_regions())

def translate_df():

	df=pd.read_csv('/../other_resources/PsychEncode_Ensembl_NoZero.csv')

	df=get_ctx_df(df)
	#print (df)

	df.rename(columns=lambda x: x[:x.index('.DFC')], inplace=True)

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
	#print (df)
	df.to_csv('%s_trajectory.csv'%name)

	df=df.reset_index()
	plt.plot(df['Age'], df['mean'])
	#sns.lineplot(x='Age', y='mean', data=df, err_style='band')
	plt.ylim([0, 70])
	plt.fill_between(df['Age'], df['mean']-df['stderr'], df['mean']+df['stderr'], alpha=0.5)
	plt.xscale('log')
	plt.show()



def find_fetal_trajectory():
	#find_fetal_diff()
	#find_adult_diff()
	df=translate_df()
	#syndromic=find_sfari_syndromic_genes()
	#overlap=['HSPA5', 'KIF3A', 'NECAP1', 'FARP1', 'STAU2', 'HSPB1', 'TUBG1', 'PFKM', 'ILF3', 'KIDINS220', 'TTC7B', 'OCIAD1', 'ANXA2', 'FN1', 'RELA', 'PI4KA', 'ATP2B1', 'YES1', 'STRAP', 'NDRG1', 'EFTUD2', 'LMNA', 'CAD', 'GDI1', 'LIMA1', 'DLG5', 'RAP1GDS1', 'UBC', 'ACTL6B', 'PRPF8', 'PFKP', 'KRT17', 'ESYT1', 'PTPRO', 'TRIM28', 'AGO1', 'PFKL', 'FTH1', 'CNOT1', 'ALDH7A1', 'KRT5', 'SEC16A', 'SLC25A6', 'RUNDC3A', 'UBQLN2', 'PRPF19', 'TUBB', 'TMPO', 'SPTBN2', 'RBMX', 'RAI14', 'SBF2', 'SNRNP200', 'SNRPN', 'APC2', 'PPP1R9A', 'PLEKHA5', 'HCFC1', 'DDAH1', 'AGAP1', 'DNAJC13', 'DYRK1A', 'CIRBP', 'PSMD7', 'THRAP3', 'CLTC', 'DSP']
	#overlap=['KRT17', 'TMPO', 'UBQLN2', 'HCFC1', 'OCIAD1', 'CNOT1', 'RAI14', 'SNRNP200', 'AGO1', 'SEC16A', 'EFTUD2', 'PRPF19', 'THRAP3', 'CIRBP', 'RUNDC3A', 'SLC25A6', 'ILF3', 'TUBG1', 'DYRK1A', 'ESYT1', 'TRIM28', 'PRPF8', 'PSMD7', 'FN1', 'SNRPN']
	#overlap=['TACC2', 'ESYT1', 'PSMD7', 'SLC25A6', 'UBE2D3', 'TPM4', 'CMAS', 'DYRK1A', 'DCX', 'TRIM28', 'PSMD5', 'CIRBP', 'RUVBL1', 'DDX39B', 'PPME1', 'AGO1', 'KRT17', 'SNRPN', 'LARP1', 'RBM14', 'KRT10', 'FARSA', 'HCFC1', 'ILF3', 'PRPSAP2', 'RAI14', 'CNOT1', 'SSB', 'SEC16A', 'STMN2', 'TUBG1', 'FN1', 'EIF4A1', 'DHX9', 'PRPF19', 'RNH1', 'OCIAD1', 'HNRNPC', 'EFTUD2', 'BAG6', 'RUNDC3A', 'EMD', 'SNRNP200', 'TMPO', 'THRAP3', 'UBQLN2', 'PRPF8', 'RALGAPA1', 'REPS1', 'TARDBP', 'PLD3', 'CLINT1', 'NUP93', 'NEDD4L', 'PSMD4', 'SMG1']
	overlap=['SUPT16H', 'XRCC5', 'SSRP1', 'XRCC6', 'TMOD3', 'FAT3', 'KHDRBS1', 'PCDH10', 'RAI14', 'ILF2', 'SLC25A1', 'PRMT1', 'OCIAD1', 'CIRBP', 'UBQLN2', 'G3BP1', 'UBQLN1', 'APEX1', 'SLC25A6', 'TUBG1', 'SNRNP200', 'RUNDC3A', 'SNRPN', 'SRSF1', 'TMPO', 'DHX15', 'PRPF19', 'PRPF8', 'AGO1']

	df=pd.DataFrame({'genes': overlap})
	df.to_csv('fetal_specific.csv')
	synd_df=find_overlap_df(overlap)
	#print (synd_df)
	new=synd_df.transpose()
	new=new.groupby(new.index).mean()
	new = new.groupby(by=new.columns, axis=1).mean()
	plot_age_mean(new, "fetal_specific")

def find_adult_trajectory():
	#find_fetal_diff()
	#find_adult_diff()
	df=translate_df()
	#syndromic=find_sfari_syndromic_genes()
	#overlap=['JPH3', 'CAMSAP2', 'MID2', 'SHC3', 'RALBP1', 'ATP2C1', 'ATP6AP2', 'CACNA1S', 'SRRM1', 'NECAB1', 'TCOF1', 'FKBP3', 'SERPINH1', 'BAG4', 'OTUD7A', 'UBE2V1', 'PPIP5K1', 'AGO2', 'CLIP3', 'FKBP8', 'CDK14', 'MTMR7', 'ELMO1', 'NDRG4', 'SRRM2', 'USP11', 'HDAC6', 'UHRF1BP1L', 'SULT4A1', 'KCNB2', 'TUBG2', 'COBL', 'SEH1L', 'SPHKAP', 'ZC3H15', 'CPSF7', 'VPS13A', 'NDUFAF4', 'PIK3R1', 'CDS2', 'PTPRR', 'PRKAG2', 'SLC7A14', 'PPP1R16B', 'CLVS2', 'MPP1', 'HABP4', 'HYOU1', 'GPR162', 'TTBK1', 'PDK3', 'KIFC2', 'ANLN', 'FERMT2', 'FCHO1', 'PREPL', 'CELF4', 'YTHDC2', 'REEP2', 'JPH1', 'PTGES2', 'RASGRF1', 'JPH4', 'RBFOX3', 'UNC5A', 'KCNQ2', 'CPNE9', 'SGTB', 'SLC25A23', 'CDK5R2', 'MAP2K4', 'KSR2', 'DTNB', 'SLC9A1', 'TBC1D5', 'UBL4A', 'ABCE1', 'BRAF', 'PRRC2A', 'JAKMIP1', 'TECR', 'MYH14', 'KRT222', 'DNAJC7', 'FRMD4A', 'RGS20', 'FAM126B', 'KHDRBS2', 'KCNIP4', 'PDXK', 'EHBP1', 'ATXN2', 'HSD17B10', 'RELL2', 'UBE2M', 'AIFM1', 'KCNH7', 'STMN3', 'CBL', 'PHACTR4', 'MAST1', 'NMT1', 'CCM2', 'CCKBR', 'RANBP2', 'RANGAP1', 'MAP7D2', 'ACAA2', 'MATK', 'KIF3C', 'MCM4', 'PYCR2', 'ISCA1', 'ARHGEF12', 'NAV3', 'PEAK1', 'DDRGK1', 'NDRG3']
	#with high-confidence 10 adult:
	#overlap=['OTUD7A', 'ANLN', 'SLC9A1', 'PRKAG2', 'JPH3', 'CPNE9', 'MATK', 'PTGES2', 'REEP2', 'PREPL']
	#overlap=['SLC9A1', 'AIFM1', 'BRAF', 'MAST1', 'MYH14', 'CAMSAP2', 'OTUD7A', 'HSD17B10', 'CPNE9', 'ABCE1', 'FAM126B', 'PREPL', 'TECR', 'JPH3', 'AGO2', 'SULT4A1', 'PTGES2', 'ANLN', 'PRRC2A', 'MAP7D2', 'REEP2', 'PRKAG2', 'MATK', 'FKBP8', 'PDXK', 'ATXN2']
	#overlap=['PREPL', 'SULT4A1', 'CAMSAP2', 'REEP2', 'FAM131B', 'TECR', 'TCEAL5', 'SLITRK5', 'ANLN', 'JPH3', 'MAST1', 'MYL12B', 'MAP7D2', 'AGO2', 'HSD17B10', 'MATK', 'CPNE9', 'ABCE1', 'AIFM1', 'PRRC2A', 'FAM126B', 'OTUD7A']
	overlap=['ANLN', 'OTUD7A', 'REEP2', 'MATK', 'TCEAL5', 'MYL12B', 'PREPL', 'CPNE9', 'JPH3', 'SLITRK5']
	print ('adult', len(overlap))
	df=pd.DataFrame({'genes': overlap})
	df.to_csv('adult_specific.csv')
	synd_df=find_overlap_df(overlap)
	#print (synd_df)
	new=synd_df.transpose()
	new=new.groupby(new.index).mean()
	new = new.groupby(by=new.columns, axis=1).mean()
	plot_age_mean(new, "adult_specific")

if __name__ == '__main__':
	find_fetal_trajectory()
	find_adult_trajectory()
