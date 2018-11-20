import mapper,cmappertools
import networkx as nx
from itertools import combinations
from networkx.drawing.nx_agraph import graphviz_layout
from __future__ import division
from itertools import product
import MapperTools as mt
bin_TEST = {}
for nbins in  range(1,50):
    bins_0 = mt.percentile_bins(fil_0, q = nbins/100, overlap = 0.)
    bins_1 = mt.percentile_bins(fil_1, q = nbins/100, overlap = 0.)
    bins_dict={'pca_comp_0':bins_0, 'pca_comp_1':bins_1}
    #
    entropy = 0
    bin_size = 0
    for (a,b),(i,j) in product(bins_dict['pca_comp_0'],bins_dict['pca_comp_1']):
        idx_0 = fil_0.index[((fil_0 <= b) & (fil_0 >= a)).nonzero()]
        idx_1 = fil_1.index[((fil_1 <= j) & (fil_1 >= i)).nonzero()]
        current_bin = idx_0.intersection(idx_1)
        if current_bin.tolist():
            bin_size += current_bin.size
    bin_TEST[nbins] = bin_size/(len(bins_dict['pca_comp_0'])*len(bins_dict['pca_comp_1']))

q = np.array(bin_TEST.keys())
PASS_nbins = []
bin_size= np.array(bin_TEST.values())
for s in np.unique(bin_size[bin_size>=5]): #window size at least 5 points
    PASS_nbins.append(q[bin_size==s].max())
    PASS_nbins.append(q[bin_size==s].min())
PASS_nbins = np.array(sorted(list(set(PASS_nbins))))

from itertools import product
from scipy.spatial import distance
from scipy.stats import entropy
#
from MapperTools import cluster_DBSCAN
def th_avdeg(L,node_info):
    tot_th_deg = 0
    if sqrt(len(L))% 1 !=0:
        raise ValueError('number of level set indices is not correct')
    from itertools import product
    def neigh_numbers(x):
        nn = 0
        for i,j in product([x[0]-1,x[0],x[0]+1],[x[1]-1,x[1],x[1]+1]):
            if (i,j)==x:
                continue
            if (i<0) or (j<0) or (i>=sqrt(len(L))) or (j>=sqrt(len(L))):
                continue
            nn += len(L[(i,j)])
        return nn
    for x,ns in L.iteritems():
        if (x[0] == 0) or (x[0] == sqrt(len(L))-1) or (x[1] == 0) or (x[1] == sqrt(len(L))-1):
            if (x[0] == x[1]) or (x == (0,sqrt(len(L))-1)) or (x == (sqrt(len(L))-1,0)):
                tot_th_deg += neigh_numbers(x)*len(ns)
            else:
                tot_th_deg += neigh_numbers(x)*len(ns)
        else:
            tot_th_deg += neigh_numbers(x)*len(ns)
    return tot_th_deg/len(node_info)
#
mapper_medsize = {}#
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
from itertools import product
#
filter_dict={'pca_comp_0':fil_0, 'pca_comp_1':fil_1}
store_gen=pd.DataFrame(distance.squareform(dis_list),index= index_str_order, columns= index_str_order)
cluster="DBSCAN"
for nbins,overlap in  product(PASS_nbins,range(5,90,5)):
    print nbins, overlap
    bins_0 = mt.percentile_bins(fil_0, q = nbins/100, overlap = overlap/100)
    bins_1 = mt.percentile_bins(fil_1, q = nbins/100, overlap = overlap/100)
    bins_dict={'pca_comp_0':bins_0, 'pca_comp_1':bins_1}
    adja,node_info = mapper_2D_density(store_gen,filter_dict,bins_dict,method=cluster,metric='precomputed')
    #---
    gg=nx.Graph()
    gg.add_nodes_from(node_info.keys())
    gg.add_edges_from(adja.keys())
    #adja,node_info,level = dist_TEST[(nbins,overlap)]['euclidean']
    #
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
EUCLIDEAN.to_csv('Mapper_parameter_study_dopamine.csv')
vfunc = np.vectorize(lambda x: 1./x)
