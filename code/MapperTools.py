from collections import defaultdict
from itertools import combinations_with_replacement,product

def mapper_2D_density(data_dict,filter_dict,bins_dict,method= "kmeans", verb= False, n_edges= False, bins= False,**kwargs):
    """
    Computing the mapper summary in 2 dimensions.

    Parameters
    ----------

    data_dict: pandas.DataFrame
        Input data to perform the clustering on.
    filter_dict: dict(pandas.Series)
        pandas.Series() with indices = data_dict.indices, containg the filter values.
    bins_dict: dict(list(tuples))
        The output of the choosen mapper_tools.bins function.
    method: (Optional)
        Clustering method to be used. Can be {"kmeans","affinity","HDBSCAN"}.

    Returns
    -------
    adja : defaultdict(dict)
        the keys are tuples (i,j) where i,j is and existing edge in the graph. the values are the wights of the edges computed as the numer of elements in common between clusters i and j.
    node_info : defaultdict(dict)
        the keys are the nodes of the graph. the values are the set of data points in each node.
    """
    def creat_adja_notefficient(node_info, level_idx, **kwargs):
        adja = dict()
        from itertools import product
        levels=level_idx
        #
        TOTAL_MAPPER_EDGES=0
        check=[]
        for i,j in product(range(len(bins_0)),range(len(bins_1))):
            if j==0:
                if i==len(bins_0)-1:
                    poss_edges=[(i,j+1)]
                else:
                    poss_edges=zip([i+1,i+1,i],[j,j+1,j+1])
            elif j==len(bins_1)-1:
                if i==len(bins_0)-1:
                    continue
                else:
                    poss_edges=zip([i+1,i+1],[j-1,j])#i+1 only
            else:
                if i==len(bins_0)-1:
                    poss_edges=[(i,j+1)]
                else:
                    poss_edges=zip([i+1,i+1,i+1,i],[j-1,j,j+1,j+1])
            this_lvl_node = level_idx[(i,j)]
            for (i_compare,j_compare) in poss_edges:
                TOTAL_MAPPER_EDGES += len(levels[(i,j)])*len(levels[(i_compare,j_compare)])
                compare_lvl_node = level_idx[(i_compare,j_compare)]
                for i1 in compare_lvl_node:
                    for i2 in this_lvl_node:
                        a = int(i1)
                        b = int(i2)
                        if len(node_info[a] & node_info[b]) > 0:
                            adja[(a, b)] = len(node_info[a] & node_info[b])
                            adja[(b, a)] = len(node_info[a] & node_info[b])
        if n_edges:
            return adja.copy(),TOTAL_MAPPER_EDGES*2
        else:
            return adja.copy()
    #---------------------------
    if type(filter_dict)==pd.core.series.Series:
        return mapper_single_density(data_dict,filter_dict,bins_dict,points=points,method = method,verb=verb)
    elif len(filter_dict.keys())==1:
        return mapper_single_density(data_dict,filter_dict[filter_dict.keys()[0]],bins_dict[bins_dict.keys()[0]],points=points,method = method,verb=verb)
    #---------------------------
    global present_samples,undefined_samples
    present_samples = []
    undefined_samples = []
    #---------------------------
    node_info = defaultdict(dict) #--dict{node <- list of data samples in the node}
    level_idx = defaultdict(dict) #--dict{bin id <- list of graph nodes created in this bin}
    #---------------------------
    global filter_fcn_0, filter_fcn_1
    if len(bins_dict[bins_dict.keys()[0]])<len(bins_dict[bins_dict.keys()[1]]):
        #---------------------------the first of the nesting for loop is going to be the longest
        filter_fcn_1 = filter_dict[filter_dict.keys()[0]].copy()
        filter_fcn_0 = filter_dict[filter_dict.keys()[1]].copy()
        bins_1 = bins_dict[bins_dict.keys()[0]]
        bins_0 = bins_dict[bins_dict.keys()[1]]
    else:
        filter_fcn_0 = filter_dict[filter_dict.keys()[0]].copy()
        filter_fcn_1 = filter_dict[filter_dict.keys()[1]].copy()
        bins_0 = bins_dict[bins_dict.keys()[0]]
        bins_1 = bins_dict[bins_dict.keys()[1]]
    for level_i in range(len(bins_0)):
        #---------------------------get the list of indices in window level_i for the x axis
        idx_0_i = filter_fcn_0[(filter_fcn_0 <= bins_0[level_i][1]) & (filter_fcn_0>= bins_0[level_i][0])].index
        if len(idx_0_i)==0:#---------if bin is empty don't work on this axis anymore
            if verb: print "\x1b[31m no points in the %d(%f) row\x1b[0m"%(level_i,bins_0[level_i][0])
            for level_j in range(len(bins_1)):
                level_idx[(level_i,level_j)]=[]
            continue
        #-------------------------------------------CONSTRUCTING THE RECTANGULAR BINS
        for level_j in range(len(bins_1)):
            #-----------------------------------get the list of indices in window level_j for the y axis
            idx_1_j = filter_fcn_1[(filter_fcn_1 <= bins_1[level_j][1]) & (filter_fcn_1 >= bins_1[level_j][0])].index
            if len(idx_1_j)==0:#check_i==1
                if verb: print "\x1b[31m no points in the (%d,%d) quadrant\x1b[0m"%(level_i,level_j)
                level_idx[(level_i,level_j)]=[]
                continue
            #----------------------------------------------COMPUTING CLUSTERS & CONSTRUCTING THE NODE INFO DICTIONARY
            result = defaultdict(dict) #--where the cluster is going to be saved
            level = (level_i,level_j) #---current bin id
            idx = idx_0_i & idx_1_j #-----samples in the bin

            result = inside_mapper_pb(idx, data_dict, level_idx, node_info, method, level, (bins_0[level_i][0],bins_1[level_j][0]),verb=verb,**kwargs)
            node_info,level_idx[level] = result
            del result, idx, level
    #---------------------------------------------END FOR
    #---------------------------adding samples not in any cluster to node_info
    n_nodes = len(node_info)
    for j_,n in enumerate(list(set(undefined_samples).difference(set(present_samples)))):
        j_ += 1
        curr_color = j_ + n_nodes
        if curr_color in node_info:
            raise ValueError('Please raise this issue in the git repo. It should not happen.')
        node_info[curr_color] = set([n])
    #---------------------------
    print('done.')
    if bins & n_edges:
        adja, edges = creat_adja_notefficient(node_info, level_idx)
        return adja,node_info,level_idx, edges
    elif bins & (not n_edges):
        adja = creat_adja_notefficient(node_info, level_idx)
        return adja,node_info,level_idx
    elif n_edges & (not bins):
        adja, edges = creat_adja_notefficient(node_info, level_idx)
        return adja,node_info,edges
    else:
        adja = creat_adja_notefficient(node_info, level_idx)
        return adja,node_info

def inside_mapper_pb(idx, data_dict, level_idx, node_info, method, level, b, verb=False, forced_cutoff=True, **kwargs):
    #----------------- ARE THERE ENOUGH POINTS TO APPLY CLUSTERING?
    if len(idx)<5:
        if verb: print "\x1b[31m less 5 points \x1b[0m"
        if verb: print ("\tEstimated number of clusters: %d' "% len(idx))
        for n in set(idx):
            if not n in present_samples:
                undefined_samples.append(n)
        return node_info,level_idx[level]
    #-----------------COMPUTE THE CLUSTERING
    now=data_dict.loc[idx] #-- select the data present in the bin
    if method== "kmeans":
        labels, n_clusters_, L =cluster_kmeans(now,range(1,min([len(idx),8])),**kwargs)
    elif method== "affinity":
        if kwargs.get("affinity")=='precomputed':
            now=now[idx]
        X=now.as_matrix()
        labels, n_clusters_=cluster_affinity(X,**kwargs)
    elif method == "HDBSCAN":
        if kwargs.get("metric")=='precomputed':
            now=now[idx]
        X=now.as_matrix()
        labels, n_clusters_=cluster_HDBSCAN(X,**kwargs)
    elif method == "DBSCAN":
        if forced_cutoff:
            if kwargs.get("metric")=='precomputed':
                X=now[idx].as_matrix()
                distance=X
            else:
                X=now.as_matrix()
                distance=pairwise_distances(X,X,**kwargs)
            measure=distance[triu_indices_from(distance,k=1)]
            cut_off=np.average(measure)
            labels, n_clusters_=cluster_DBSCAN(X,eps=cut_off,**kwargs)
        else:
            if kwargs.get("metric")=='precomputed':
                now=now[idx]
            X=now.as_matrix()
            labels, n_clusters_=cluster_DBSCAN(X,**kwargs)
    else:
        raise AttributeError('The chosen clustering method is not valid')
    #-------------------- IF THERE ARE NO CLUSTERS
    if verb: print labels
    if verb: print n_clusters_
    level_idx[level] = []
    if n_clusters_==0:
        if verb: print "\x1b[31m no cluster \x1b[0m"
        undefined_samples.extend(idx)
        return node_info, level_idx[level]
    #-------------------- ADDING NODES TO THE GRAPH
    if verb: print ("\tEstimated number of clusters: %d' "% n_clusters_)
    num_graph_nodes = len(node_info) #last node id added to the graph
    #--------------------
    assert type(labels) == numpy.ndarray, "the clustering algorithm does not return a numpy array. type {}".format(type(labels))
    for k in set(labels):
        if k==-1: # sklearn signs as -1 the outliers that could not be part of a cluster
            for n in set(now.index[labels == k]):
                if not n in present_samples:
                    undefined_samples.append(n)
        else:
            curr_node_id = k + num_graph_nodes
            level_idx[level].append(curr_node_id)
            node_info[curr_node_id] = set(now.index[labels == k]);
            present_samples.extend(node_info[curr_node_id])
    if verb: print 'nodes in the graph before:',num_graph_nodes,'and now:', len(node_info)
    return node_info,level_idx[level]

def create_node_matrix(node_info):
    data={}
    for n,d in node_info.iteritems():
        for a,b in combinations_with_replacement(d,2):
            data.setdefault(a,dict()).setdefault(b,0)
            data[a][b] +=1
            if a!=b:
                data.setdefault(b,dict()).setdefault(a,0)
                data[b][a] +=1
    return data

def cluster_HDBSCAN(X, **kwargs):
    '''
    Computes the hdbscan clustering method for data matrix X using the hdbscan module.

    Parameters
    ----------
    X : numpy.appary
        Multidimension numpy.array of data points

    Returns
    -------
    labels, n_clusters_
    '''
    import hdbscan
    #cluster_selection_method='leaf'
    clusterer = hdbscan.HDBSCAN(min_cluster_size=2,min_samples=1, prediction_data=True, **kwargs).fit(X)
    #print clusterer.prediction_data_
    #soft_clusters = hdbscan.all_points_membership_vectors(clusterer)
    labels_ = clusterer.labels_#[np.argmax(x)+1 for x in soft_clusters]
    #translate = dict(zip(set(labels),range(len(set(labels)))))
    n_clusters_ = len(set(labels_))
    #labels_ = map(lambda x: translate[x], labels)
    return np.array(labels_), n_clusters_

def mapper_single_density(data_dict,filter_fcn,bins,points="row",method = "kmeans",verb=False,**kwargs):
    '''
    Computes the mapper graph for data in data_dict, using filters in filter_dict, and bins in bins_dict.

    Parameters
    ----------
    data_dict: pandas.DataFrame
        Input data to perform the clustering on.
    filter_fcn: pandas.Series
         pandas.Series() with indices = data_dict.indices, containg the filter values.
    bins: dict(list(tuples))
        The output of the choosen mapper_tools.bins function.
    points(OPTIONAL): str
        Which axis has to be treated as patients, default="row".
        If "column" is given then the algorithm will be computed for the transposed DataFrame.
    method : str
        Clustering method to be used. Can be {"kmeans","affinity","HDBSCAN"}.
        Refer to sklearn documentation for more information on the differences between the methods.

    Returns
    -------
    adja : defaultdict(dict)
        the keys are tuples (i,j) where i,j is and existing edge in the graph. the values are the wights of the edges computed as the numer of elements in common between clusters i and j.
    node_info : defaultdict(dict)
        the keys are the nodes of the graph. the values are

    '''
    #---------------------------
    present_samples = []
    undefined_samples = []
    #---------------------------
    node_info = defaultdict(dict)
    adja = defaultdict(dict)
    level_idx = defaultdict(dict)
    if points=="column":
        eg=eg.T
    #---------------------------
    for level_i in range(0,len(bins)): #filter_dict.keys()[0]
        #---------------------------constructing bin on the x axis with range i
        idx = filter_fcn[(filter_fcn <= bins[level_i][1]) & (filter_fcn>= bins[level_i][0])].dropna().index
        if len(idx)==0:#---------if bin is empty don't work on this axis anymore
            if verb: print "\x1b[31m no points in the %d(%f) bin\x1b[0m"%(level_i,bins[level_i][0])
            if verb: print "nothing to do here, moving onto the next level"
            level_idx[level_i] = []
            continue
        #------------#------------#------------#------------#------------#------------#------------#------------
        num_points = len(idx)
        if verb: print 'number of points in level', level_i,':', num_points
        #------------#------------#------------#------------#------------#------------#------------#------------
        dam=0.8
        now=data_dict.loc[idx]
        if num_points<5:
            labels, n_clusters_ = range(num_points), num_points
        elif method == "kmeans":
            X=now.as_matrix()
            if verb: print "using kmeans"
            labels, n_clusters_, rand__=cluster_kmeans(now,range(1,5),**kwargs)
        elif method == "affinity":
            try:
                now=now[idx]
            except:
                print 'sanity check: i could not make matrix square'
            X=now.as_matrix()
            labels, n_clusters_=cluster_affinity(X,**kwargs)
        elif method == "HDBSCAN":
            X=now.as_matrix()
            if verb: print "using HDBSCAN"
            labels, n_clusters_=cluster_HDBSCAN(X,**kwargs)
        elif method == "DBSCAN":
            X=now.as_matrix()
            if verb: print "using DBSCAN"
            if kwargs.get("metric")=='precomputed':
                X=now[idx].as_matrix()
            labels, n_clusters_=cluster_DBSCAN(X,**kwargs)
        elif method == "agglomerative":
            X=now.as_matrix()
            if verb: print "using Agglomerative clustering with cosine similarity"
            labels, n_clusters_=cluster_agglomerative(X,**kwargs)
        elif method == "spectral":
            X=now.as_matrix()
            if verb: print "using spectral clustering with cosine similarity"
            labels, n_clusters_=cluster_spectral(X,**kwargs)
        else:
            raise ValueError('The chosen clustering method is not valid.')
        #if n_clusters_==0:
        #    if len(labels)<1:
        #        if verb: print 'no clusters'
        #        level_idx[level_i] = []
        #        continue

        from itertools import cycle

        if verb: print ("\tEstimated number of clusters: %d' "% n_clusters_)

        node_colors = labels#color_graph(G)
        num_colors = len(set(node_colors))

        num_graph_nodes = len(node_info)#init at level 1 then adds up the old sets
        level_idx[level_i] = []
        correction = 0
        if verb: print 'prima del for che crea node info', node_colors, num_colors
        for j in sorted(set(node_colors),reverse=True):
            if j==-1:
                for n in set(now.index[node_colors == j]):
                    if not n in present_samples:
                        undefined_samples.append(n)
            else:
                curr_color = j + num_graph_nodes
                if verb: print curr_color
                level_idx[level_i].append(curr_color)
                #node_info[curr_color]['level'] = level_i
                #node_info[curr_color]['fnval'] = bins[level_i][0]
                #print (j,len(now.index[node_colors == j]))
                node_info[curr_color] = set(now.index[np.array(node_colors) == j])
                present_samples.extend(node_info[curr_color])
        if level_i > 0:
            prev_lvl_idx = level_idx[level_i - 1]
            this_lvl_idx = level_idx[level_i]
            #print ('prev and this _lvl_idx',prev_lvl_idx,this_lvl_idx)
            for i1 in prev_lvl_idx:
                for i2 in this_lvl_idx:
                    a = int(i1)
                    b = int(i2)
                    #print('\t',a,b,node_info[a]['set'] & node_info[b]['set'])
                    if len(node_info[a] & node_info[b]) > 0:
                        #print ('\t',a,b,True)
                        adja[(a, b)] = len(node_info[a] & node_info[b])
                        adja[(b, a)] = len(node_info[a] & node_info[b])
    #---------------------------adding undefined_samples to node_info
    for j_,n in enumerate(list(set(set(undefined_samples).difference(set(present_samples))))):
        j_+=1
        curr_color = j_ + len(node_info)
        node_info[curr_color] = set([n])
    print('done')
    return adja,node_info
