#__author__ = 'Alice Patania'
import cPickle as pk
import pandas as pd
import numpy as np
import os
import MapperTools
from sklearn.decomposition import PCA
from scipy.spatial import distance
import argparse, sys


if __name__ == '__main__':
    parser=argparse.ArgumentParser()

    parser.add_argument('--name', help='name of the gene list to run the data on. Options: dopamine, richiardi, full')

    args=parser.parse_args()

    idx_= args.name

    # idx_='dopamine' # tag that is gonna be used to save the output

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
    filter_dict={'comp_0':fil_0, 'comp_1':fil_1}
    del V, pca_score, pca, df

    dis_list = distance.pdist(store_gen, metric = 'correlation')

    with open('./csv_outputs/Parameters_'+idx_+'.txt', "r") as  text_file:
        lines = text_file.read().lstrip('[(').rstrip(')]').split('), (')
    params = [tuple(map(int,s.split(','))) for s in lines]
    del lines

    store_gen=pd.DataFrame(distance.squareform(dis_list),index= index_str_order, columns= index_str_order)
    import os
    if not os.path.isdir('../ouput'):
        od.makedir('../ouput')
    for nbins,overlap in params:
        print nbins, overlap
        bins_0 = mt.percentile_bins(fil_0, q = nbins/100, overlap = overlap/100)
        bins_1 = mt.percentile_bins(fil_1, q = nbins/100, overlap = overlap/100)
        bins_dict={'comp_0':bins_0, 'comp_1':bins_1}
        store_gen=pd.DataFrame(distance.squareform(dis_list),index= index_str_order, columns= index_str_order)
        adja,node_info = mapper_2D_density(store_gen,filter_dict,bins_dict,method="DBSCAN", metric = 'precomputed')
        pk.dump(node_info,open('../ouput'+idx_+'{}_{}_node_info.pk'.format(nbins,overlap),'w'))
        pk.dump(adja,open('../ouput'+idx_+'{}_{}_adja.pk'.format(nbins,overlap),'w'))
