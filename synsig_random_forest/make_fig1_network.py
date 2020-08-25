import pandas as pd


df=pd.read_csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/pred_genes_above_4.7.csv')

print (df)

genes=df['genes'].tolist()

new=pd.DataFrame({'Root': 'Genes', 'Genes': genes})

print (new)

new.to_csv('Fig1_network.csv')