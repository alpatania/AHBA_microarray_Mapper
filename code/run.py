#__author__ = 'Alice Patania'
import cPickle as pk
import pandas as pd
import numpy as np
import os
from sklearn.decomposition import PCA
from scipy.spatial import distance

idx_='dopamine' # tag that is gonna be used to save the output
path_to_gene_list='../data/'+idx_+'_list.txt' #substitute with yours
verbose = False

left_store_gen=pk.load(open('ABA/microarray_input_mapper.pk','r'))

f=open(path_to_gene_list, 'r')
gene_list=numpy.loadtxt(f,dtype=str)# reads the txt and creates and array with the list of genes
f.close()
store_gen=left_store_gen[gene_list].copy()
del path_to_gene_list, left_store_gen, f

df=store_gen.copy()
pca = PCA(n_components=2).fit(df)
pca_score = pca.explained_variance_ratio_
V = pca.components_
if verbose: print 'ratio:',pca_score, '--- variance:',pca.explained_variance_
fil_0=df.dot((V[0,:]).transpose())
fil_1=df.dot((V[1,:]).transpose())
del V, pca_score, pca, df

dis_list = distance.pdist(store_gen, metric = 'correlation')
