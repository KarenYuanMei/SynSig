#goal: define the training and test genes from the input examples

import csv

import pandas as pd
import random


#load genes from a csv file
def load_genes(filename):
	df=pd.read_csv(filename)
	genes=df['genes'].tolist()
	genes.sort()
	return genes

#this function divides the genes into 5 sections to enable 5-fold cross-validation
def divide_positives_negatives(positives, negatives):
	chunk_size=int(len(positives)/5)
	pos_chunks = [positives[i * chunk_size:(i + 1) * chunk_size] for i in range((len(positives) + chunk_size - 1) // chunk_size )]
	neg_chunks = [negatives[i * chunk_size:(i + 1) * chunk_size] for i in range((len(negatives) + chunk_size - 1) // chunk_size )]
	#print (pos_chunks)
	#print (neg_chunks)
	return pos_chunks, neg_chunks

def define_training_test(pos_chunks, neg_chunks, chunk_no):
	
	test_pos=pos_chunks[chunk_no]
	training_pos=list(set(positives)-set(test_pos))

	test_neg=neg_chunks[chunk_no]
	training_neg=list(set(negatives)-set(test_neg))

	training=training_pos+training_neg
	training.sort()
	training_df=pd.DataFrame(training, columns=['genes'])
	training_df.to_csv('training_genes_%s.csv'%chunk_no)
	test=test_pos+test_neg
	test.sort()
	test_df=pd.DataFrame(test, columns=['genes'])
	test_df.to_csv('test_genes_%s.csv'%chunk_no)
	
	print ('overlap', len(set(training)&set(test)))
	return training, test



#unit test:------------------------------------------------------------------------------------------------
def load_test():

	tests=[]

	for i in range(5):
		test=pd.read_csv('test_genes_%s.csv'%i)
		test=test['genes'].tolist()
		print (len(test))
		tests.append(test)

	overlap=set(tests[0]).intersection(*tests)
	print (overlap)

	#print (tests)
	flat_tests=[item for sublist in tests for item in sublist]
	print (len(set(flat_tests)))

	return overlap


#main function--------------------------------------------------------------
if __name__ == "__main__":
	
	#load the synapse positives and negatives files

	positive_file='synapse_positives.csv'
	negative_file='synapse_negatives.csv'

	positives=load_genes(positive_file)
	print (positives[:5])
	negatives=load_genes(negative_file)
	print (negatives[:5])
	
	#randomly shuffle the positive and negative gene lists
	random.seed(4)
	random.shuffle(positives)
	random.shuffle(negatives)

	#divide the positive and negative gene lists into training and test:

	pos_chunks, neg_chunks=divide_positives_negatives(positives, negatives)
	#print (pos_chunks)
	#print (neg_chunks)
	n=[0,1,2,3,4]
	for item in n:
		training, test=define_training_test(pos_chunks, neg_chunks, item)
		print (training[:5])
		print (test[:5])
	#run unit test
	print (load_test())