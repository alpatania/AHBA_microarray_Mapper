import networkx as nx
import pandas as pd
import pickle as pk

def shortest_path(g,seed=None,max_iter=1000):
    '''
    computes the shortest path from a seed node to every other node in the graph.

    Parameters
    ----------
    g: nx.Graph
    seed: id starting node in the graph

    Returns
    -------
    short_path: dict
    '''
    C=list(nx.connected_components(g))
    max_cc=sorted(C,key=len,reverse=True)[0]
    if seed is None:
        seed=list(max_cc)[0]
    time_node_old=[seed]
    short_path={seed:[0]}
    t=1
    while t<max_iter:
        time_node_new=[]
        for n in time_node_old:
            for m in g.neighbors(n):
                if m not in short_path:
                    short_path[m]=[t]
                    time_node_new.append(m)
        time_node_old=time_node_new
        if len(short_path.keys())==len(max_cc):
            print 'network percolated'
            break
        t+=1
    return short_path

if __name__ == '__main__':
    parser=argparse.ArgumentParser()

    parser.add_argument('--name', help='name of the gene list to run the data on. Options: dopamine, richiardi, full')

    args=parser.parse_args()

    idx_= args.name

    verbose = False

    with open('./csv_outputs/Parameters_'+idx_+'.txt', "r") as  text_file:
        lines = text_file.read().lstrip('[(').rstrip(')]').split('), (')
    params = [tuple(map(int,s.split(','))) for s in lines]
    del lines

    for nbins,overlap in params:
        node_info = pk.load(open('../ouput'+idx_+'{}_{}_node_info.pk'.format(nbins,overlap),'r'))
        adja = pk.load(open('../ouput'+idx_+'{}_{}_adja.pk'.format(nbins,overlap),'r'))
        #-----------Looking for the list of nodes containg elements from the seed area
        #---------------------------------
        #---------PERSONALIZE: Modify the id of the ROI to use as seed to use this code for your own purpouses.
        #---------------------------------
        seed_searcher=struct.reset_index().set_index('structure_2_id')[['donor','sample_id']]
        A=seed_searcher.ix[9066]#ventral tegmental area
        B=seed_searcher.ix[9072]#substantia nigra
        set_A=zip(A.values[:,0],A.values[:,1])
        set_B=zip(B.values[:,0],B.values[:,1])
        seed_A=[]
        seed_B=[]
        #----------Running shortest path for every node in the list area
        for n in node_info:
            if len(node_info[n].intersection(set_A))>0:
                seed_A.append(n)
            if len(node_info[n].intersection(set_B))>0:
                seed_B.append(n)
        del A,B,set_A,set_B
        gg=nx.Graph()
        gg.add_edges_from(adja.keys())

        NODE_VISIT_TW = pd.DataFrame(np.zeros((gg.number_of_nodes(),len(seed_A))),index=gg.nodes(), columns=seed_A)
        #---------------------------------
        #-------------This is commented out, but it registers which edges are most used by the shortest paths.
        #EDGE_USED = pd.DataFrame(np.zeros((gg.number_of_nodes(),gg.number_of_nodes())),index=gg.nodes(),columns=gg.nodes())
        #---------------------------------
        verbose= True
        for seed in seed_A:
            if seed in gg.nodes():
                walk_node=shortest_path(gg,seed=seed)
                NODE_VISIT_TW[seed][walk_node.keys()] = zip(*walk_node.values())[0]
            else:
                if verbose: print seed, len(seed_A)
        pk.dump(NODE_VISIT_TW,open('NODE_VISIT_SP_vta_{}_{}.pk'.format(nbins,overlap),'w'))
        NODE_VISIT_TW = pd.DataFrame(np.zeros((gg.number_of_nodes(),len(seed_B))),index=gg.nodes(), columns=seed_B)
        #---------------------------------
        #-------------This is commented out, but it registers which edges are most used by the shortest paths.
        #EDGE_USED = pd.DataFrame(np.zeros((gg.number_of_nodes(),gg.number_of_nodes())),index=gg.nodes(),columns=gg.nodes())
        #---------------------------------

        for seed in seed_B:
            if seed in gg.nodes():
                walk_node=shortest_path(gg,seed=seed)
                NODE_VISIT_TW[seed][walk_node.keys()] = zip(*walk_node.values())[0]
            else:
                if verbose: print seed, len(seed_B)
        pk.dump(NODE_VISIT_TW,open('NODE_VISIT_SP_sn_{}_{}.pk'.format(nbins,overlap),'w'))
        #---------------------------------
        #-------------This is commented out, but it registers which edges are most used by the shortest paths.
        #pk.dump(EDGE_USED,open('EDGE_USED_SP_vta_{}_{}.pk'.format(nbins,overlap),'w'))
        #---------------------------------
        del NODE_VISIT_TW
