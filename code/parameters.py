
import MapperTools as mt
import networkx as nx
import pandas as pd
import numpy as np
from networkx.drawing.nx_agraph import graphviz_layout
from itertools import combinations, product
from sklearn.decomposition import PCA
from scipy.spatial import distance
import argparse, sys
import pickle as pk
from numpy import triu_indices_from

if __name__ == '__main__':
    parser=argparse.ArgumentParser()

    parser.add_argument('--name', help='name of the gene list to run the data on. Options: dopamine, richiardi, full')

    args=parser.parse_args()

    idx_= args.name

    # idx_='dopamine' # tag that is gonna be used to save the output
    print('Starting...')
    #---------------------Getting the data
    #---------------------
    #----PERSONALIZATION-NOTE 1: to run this script on your data substitute the next block with
    #----                        store_gen = pd.DataFrame() with (n,m) with n the number of data points, m the dimension of the dataset point
    #----                        fil_0, fil_1 = pd Series() of length n, with fil_0.index == fil_1.index == store_gen.index
    #---------------------
    if idx_ != 'full':
        path_to_gene_list='../data/'+idx_+'_list.txt' #substitute with yours
    verbose = False

    left_store_gen=pk.load(open('../data/microarray_input_mapper.pk','r'))

    f=open(path_to_gene_list, 'r')
    gene_list=np.loadtxt(f,dtype=str)# reads the txt and creates and array with the list of genes
    f.close()
    store_gen=left_store_gen[gene_list].copy()
    del path_to_gene_list, left_store_gen, f
    print('Filters...')
    df=store_gen.copy()
    pca = PCA(n_components=2).fit(df)
    pca_score = pca.explained_variance_ratio_
    V = pca.components_
    if verbose: print 'ratio:',pca_score, '--- variance:',pca.explained_variance_
    fil_0=df.dot((V[0,:]).transpose())
    fil_1=df.dot((V[1,:]).transpose())
    filter_dict={'comp_0':fil_0, 'comp_1':fil_1}
    del V, pca_score, pca, df
    print('Distance...')
    dis_list = distance.pdist(store_gen, metric = 'correlation')
    print('Computing average bin size...')
    #-------------------------Computing average bin size for percentile bins with density from 1% to 50% of the total dataset size
    bin_TEST = {}
    for nbins in  range(1,50):
        #print('\t',nbins)
        bins_0 = mt.percentile_bins(fil_0, q = nbins/100., overlap = 0.)
        bins_1 = mt.percentile_bins(fil_1, q = nbins/100., overlap = 0.)
        bins_dict={'pca_comp_0':bins_0, 'pca_comp_1':bins_1}
        #print('\tsumming them')
        entropy = 0
        bin_size = 0
        for (a,b),(i,j) in product(bins_dict['pca_comp_0'],bins_dict['pca_comp_1']):
            idx_0 = fil_0.index[((fil_0 <= b) & (fil_0 >= a)).nonzero()]
            idx_1 = fil_1.index[((fil_1 <= j) & (fil_1 >= i)).nonzero()]
            current_bin = idx_0.intersection(idx_1)
            if current_bin.tolist():
                bin_size += current_bin.size
        bin_TEST[nbins] = bin_size/(len(bins_dict['pca_comp_0'])*len(bins_dict['pca_comp_1']))
    #-------------------List of bin densities that on average have enough points to run the algorithm
    print('listing winning bin size...')
    q = np.array(bin_TEST.keys())
    PASS_nbins = []
    bin_size= np.array(bin_TEST.values())
    for s in np.unique(bin_size[bin_size>=5]): #window size at least 5 points
        PASS_nbins.append(q[bin_size==s].max())
        PASS_nbins.append(q[bin_size==s].min())
    PASS_nbins = np.array(sorted(list(set(PASS_nbins))))
    #-------------------Initialize dictionaries for the collection of statistics
    print('Running the mapper...')
    mapper_medsize = {}
    mapper_avsize = {}
    mapper_avdeg = {}
    mapper_cc = {}
    mapper_nonisol_cc = {}
    mapper_density_grid_f = {}
    mapper_density_diag_f = {}
    mapper_density_grid_c = {}
    mapper_density_diag_c = {}
    mapper_std = {}
    mapper_entropy = {}
    mapper_entropy_1 = {}
    mapper_entropy_Ed = {}
    #--------------------Running the mapper for all parameters and recording the statistics
    filter_dict={'pca_comp_0':fil_0, 'pca_comp_1':fil_1}
    struct=pk.load(open('../data/sampleID_to_structure.pk','r'))
    index_str_order = struct.loc[store_gen.index].sort_values('structure_3_id').index.get_values()
    store_gen=pd.DataFrame(distance.squareform(dis_list),index= index_str_order, columns= index_str_order)
    cluster="DBSCAN"
    for nbins,overlap in  product(PASS_nbins,range(5,90,5)):
        print nbins, overlap
        bins_0 = mt.percentile_bins(fil_0, q = nbins/100., overlap = overlap/100)
        bins_1 = mt.percentile_bins(fil_1, q = nbins/100., overlap = overlap/100)
        bins_dict={'pca_comp_0':bins_0, 'pca_comp_1':bins_1}
        adja,node_info = mt.mapper_2D_density(store_gen,filter_dict,bins_dict,method=cluster,metric='precomputed')
        #---
        gg=nx.Graph()
        gg.add_nodes_from(node_info.keys())
        gg.add_edges_from(adja.keys())
        #
        from scipy.stats import entropy
        mapper_avsize[(nbins,overlap)] = len(gg.nodes())/nx.number_connected_components(gg)
        mapper_medsize[(nbins,overlap)] = np.median(map(len,nx.connected_components(gg)))
        mapper_avdeg[(nbins,overlap)] = np.average(gg.degree().values())
        mapper_nonisol_cc[(nbins,overlap)] = np.sum(np.array(map(len,nx.connected_components(gg)))>1)
        mapper_cc[(nbins,overlap)] = nx.number_connected_components(gg)
        #
        mapper_density_grid_f[(nbins,overlap)] = gg.number_of_edges()/sum(map(lambda x:2*np.floor(np.sqrt(len(x)))**2-2*np.floor(np.sqrt(len(x))),nx.connected_components(gg)))
        mapper_density_diag_f[(nbins,overlap)] = gg.number_of_edges()/sum(map(lambda x:4*np.floor(np.sqrt(len(x)))**2-6*np.floor(np.sqrt(len(x)))+2,nx.connected_components(gg)))
        mapper_density_grid_c[(nbins,overlap)] = gg.number_of_edges()/sum(map(lambda x:2*np.ceil(np.sqrt(len(x)))**2-2*np.ceil(np.sqrt(len(x))),nx.connected_components(gg)))
        mapper_density_diag_c[(nbins,overlap)] = gg.number_of_edges()/sum(map(lambda x:4*np.ceil(np.sqrt(len(x)))**2-6*np.ceil(np.sqrt(len(x)))+2,nx.connected_components(gg)))
        #
        S_cc = np.array(map(len,nx.connected_components(gg)))
        mapper_std[(nbins,overlap)] = np.std(S_cc[S_cc>1])
        mapper_entropy[(nbins,overlap)] = entropy(S_cc)
        mapper_entropy_1[(nbins,overlap)] = entropy(S_cc[S_cc>1])
        S_cc_Ed = np.array(map(lambda x: nx.number_of_edges(nx.subgraph(gg,x)),nx.connected_components(gg)))
        mapper_entropy_Ed[(nbins,overlap)] = entropy(S_cc_Ed)

        del gg, adja, node_info
    print('saving...')
    #--------------------Saving the statistics in a file
    EUCLIDEAN = pd.concat([pd.Series(mapper_avdeg, name='av. deg.'),
                           pd.Series(mapper_cc, name='numb. conn. comp.'),
                           pd.Series(mapper_nonisol_cc, name='numb. conn. comp. NOT points'),
                           pd.Series(mapper_avsize, name='av. size conn. comp.'),
                           pd.Series(mapper_std, name='stand. dev. of size conn. comp. NOT points'),
                           pd.Series(mapper_entropy_Ed, name='entropy numb. edges in conn. comp.'),
                           pd.Series(mapper_entropy_1, name='entropy numb. nodes in conn. comp.'),
                           pd.Series(mapper_density_grid_f, name='edge density grid floor(sqrt(n))'),
                           pd.Series(mapper_density_diag_f, name='edge density grid+diag floor(sqrt(n))'),
                           pd.Series(mapper_density_grid_c, name='edge density grid ceiling(sqrt(n))'),
                           pd.Series(mapper_density_diag_c, name='edge density grid+diag ceiling(sqrt(n))')],axis=1, ignore_index=False)
    import os
    if not os.path.isdir('../output'):
        os.makedirs('../output')
    EUCLIDEAN.to_csv('../output/Mapper_parameter_study_'+idx_+'.csv')
