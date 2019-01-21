import pandas as pd
import numpy as np
import argparse, sys

if __name__ == '__main__':
    parser=argparse.ArgumentParser()

    parser.add_argument('--name', help='name of the gene list to run the data on. Options: dopamine, richiardi, full')

    args=parser.parse_args()

    idx_= args.name

    verbose = False

    read_param = pd.read_csv('../git-meesss/AHBA_microarray_Mapper/output/Mapper_parameter_study_'+idx_+'.csv', index_col=[0,1])

    vfunc = np.vectorize(lambda x: 1./x)

    PASS_nbins = read_param.index.get_level_values(None).unique().tolist()

    if verbose: print(read_param.columns)

    data2 = np.array([read_param['av. size conn. comp.'][x].sort_index().values.tolist() for x in sorted(PASS_nbins, reverse= True)])
    data7 = np.array([read_param['numb. conn. comp.'][x].sort_index().values.tolist() for x in sorted(PASS_nbins, reverse= True)])
    data_nodes = np.multiply(data7,data2)
    data8 = np.array([read_param['numb. conn. comp. NOT points'][x].sort_index().values.tolist() for x in sorted(PASS_nbins, reverse= True)])
    data_isola = data7-data8
    data4 = np.multiply(vfunc(data_nodes),data_isola)

    data_lattice = np.array([read_param['edge density grid floor(sqrt(n))'][x].sort_index().values.tolist() for x in sorted(PASS_nbins, reverse= True)])

    parameters = np.where(data_lattice<.8, np.where((data4<.5),np.where(data8>1,data_lattice,0),0),0).nonzero()
    parameters = zip(map(lambda x: read_param.index.levels[0][len(read_param.index.levels[0])-x-1],parameters[0]),map(lambda x: read_param.index.levels[1][x],parameters[1]))

    np.savetxt('../data/Parameters_'+idx_+'_gene_list.txt',parameters)
