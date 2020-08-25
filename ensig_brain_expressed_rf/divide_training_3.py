import csv

import pandas as pd



def load_genes(filename):
	df=pd.read_csv(filename)
	genes=df['genes'].tolist()
	return genes

#this function divides the genes into 5 sections to enable 5-fold cross-validation
def divide_positives_negatives(positives, negatives):
	chunk_size=int(len(positives)/5)
	pos_chunks = [positives[i * chunk_size:(i + 1) * chunk_size] for i in range((len(positives) + chunk_size - 1) // chunk_size )]
	neg_chunks = [negatives[i * chunk_size:(i + 1) * chunk_size] for i in range((len(negatives) + chunk_size - 1) // chunk_size )]
	return pos_chunks, neg_chunks

def define_training_test(pos_chunks, neg_chunks, chunk_no):
	#random.shuffle(positives)
	#positives_chunks = [positives[i * 40:(i + 1) * 40] for i in range((len(positives) + 40 - 1) // 40 )]
	#print (positives_chunks)

	test_pos=pos_chunks[chunk_no]
	training_pos=list(set(positives)-set(test_pos))

	#random.shuffle(negatives)
	#neg_chunks = [negatives[i * 40:(i + 1) * 40] for i in range((len(negatives) + 40 - 1) // 40 )]

	test_neg=neg_chunks[chunk_no]
	training_neg=list(set(negatives)-set(test_neg))

	training=training_pos+training_neg
	training_df=pd.DataFrame(training, columns=['genes'])
	training_df.to_csv('training_genes_%s.csv'%chunk_no)
	test=test_pos+test_neg
	test_df=pd.DataFrame(test, columns=['genes'])
	test_df.to_csv('test_genes_%s.csv'%chunk_no)
	#print (len(training))
	#print (len(test))
	#print ('training', training)
	#print ('test', test)
	print ('overlap', len(set(training)&set(test)))
	return training, test



#unit test:------------------------------------------------------------------------------------------------
def load_test():
	test_0=pd.read_csv('test_genes_0.csv')
	test_0=test_0['genes'].tolist()

	test_1=pd.read_csv('test_genes_1.csv')
	test_1=test_1['genes'].tolist()

	test_2=pd.read_csv('test_genes_2.csv')
	test_2=test_2['genes'].tolist()

	test_3=pd.read_csv('test_genes_3.csv')
	test_3=test_3['genes'].tolist()

	test_4=pd.read_csv('test_genes_4.csv')
	test_4=test_4['genes'].tolist()
	overlap=set(test_0)&set(test_1)&set(test_2)&set(test_3)&set(test_4)
	return overlap


#main function--------------------------------------------------------------
if __name__ == "__main__":
	
	#training, test=divide_training_test(positives, negatives, 0)

	positive_file='nb_brain_synapse_positives.csv'
	negative_file='nb_brain_synapse_negatives.csv'

	positives=load_genes(positive_file)
	negatives=load_genes(negative_file)
	print (len(positives))
	print (len(negatives))

	pos_chunks, neg_chunks=divide_positives_negatives(positives, negatives)
	print (pos_chunks)
	print (neg_chunks)
	n=[0,1,2,3,4]
	for item in n:
		training, test=define_training_test(pos_chunks, neg_chunks, item)
	#hold_out(synapse_overlap, positives, negatives)
	print (load_test())