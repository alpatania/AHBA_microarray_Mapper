for nbins,overlap in params:#parameters:
    print nbins, overlap
    node_info = pk.load(open(idx_+'_{}_{}_node_info.pk'.format(nbins,overlap),'r')
     ## having the graph
    seed_searcher=struct.reset_index().set_index('structure_2_id')[['donor','sample_id']]
    A=seed_searcher.ix[9066]#ventral tegmental area
    B=seed_searcher.ix[9072]#substantia nigra
    set_A=zip(A.values[:,0],A.values[:,1])
    set_B=zip(B.values[:,0],B.values[:,1])
    seed_A=[]
    seed_B=[]
    for n in node_info:
        if len(node_info[n].intersection(set_A))>0:
            seed_A.append(n)
        if len(node_info[n].intersection(set_B))>0:
            seed_B.append(n)
    del A,B,set_A,set_B
    gg=nx.Graph()
    gg.add_edges_from(adja.keys())
    pos= graphviz_layout(gg)
    plt.figure(figsize=(15,5))
    plt.subplot(121)
    nx.draw_networkx_edges(gg,pos= pos)
    nx.draw_networkx_nodes(gg.subgraph(seed_A),pos= pos, node_size=10)
    plt.axis('off')
    plt.title('VTA')
    plt.subplot(122)
    nx.draw_networkx_edges(gg,pos= pos)
    nx.draw_networkx_nodes(gg.subgraph(seed_B),pos= pos, node_size=10)
    plt.axis('off')
    plt.title('Substantia nigra')
    plt.savefig('img_results/seeds_{}_{}.pdf'.format(nbins, overlap),format='pdf')
