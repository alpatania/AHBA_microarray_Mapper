# TDA summary for gene-expression: Mapper Algorithm in 2D
Authors: Alice Patania, Pierluigi Selvaggi, Mattia Veronese, Ottavia Dipasquale, Paul Expert, and Giovanni Petri

## Abstract
Understanding how gene expression translates to and affects human behaviour is one of the ultimate aims of neuroscience.  In this paper, we present a pipeline based on Mapper, a topological simplification tool, to produce and analyze genes co-expression data.  We first validate the method by reproducing key results from the literature on the Allen Human Brain Atlas, and the correlations between resting-state fMRI and gene co-expression maps.  We then analyze adopamine-related gene-set and find that co-expression networks produced by Mapper returned a structure that matches the well-known anatomy of the dopaminergic pathway.  Our results suggest that topological network descriptions can be a powerful tool to explore the relationships between genetic pathways and their association with brain function and its perturbation due to illness and/or pharmacological challenge.

> **DISCLAIMER:**
> Running all the scripts in this repository is going to give the list of all the results found in the paper, but not the figures or the standard exploratory analysis ( i.e. the histograms and KS tests ). 
> I am willing to change this decision if anyone needs it, write to me or start an issue
  
## To do list:
  - [x] put up the datasets
  - [x] make a parameters selection script
  - [x] make all code into scripts that can be run from command line
  - [ ] make a script to compute the agreement matrix
  - [x] make the shortest path script
  - [x] write a tutorial on how to run all the code
  - [ ] find a way to put the dataset that are too big for git
  - [ ] add the list of sample ids used by Richiardi et al. in their paper
  - [ ] add dependencies
  
## Content:
1. data:
    - dataset normalized: Download the data used in the study [here](https://figshare.com/s/9f9806df6a5a73cc18bf).
    - the two list of genes used in the study
      - `dopamine.txt`
      - `richiardi.txt`
3. code:
    - `MapperTools.py`: All the functions needed to build the graph
    - `parameters.py`: Computes the statistics used for the choice of parameters.  
    takes as input the dataset id (`dopamine`, `richiardi`, or `full`) and saves the statistics in a csv in the folder `output`.
    - `selection.py`: Selects the optimal parameters using the output from `parameters.py`.  
    takes as input the dataset id (`dopamine`, `richiardi`, or `full`) and saves the parameters in a txt in the folder `output`.
    - `run.py`: Builds the graph for the optimal parameters found by `selection.py`.  
    takes as input the dataset id (`dopamine`, `richiardi`, or `full`) and saves the adjacency matrix and node information in 2 pickled dictionaries in the folder `output`.
    - `agreement_matrix.py`: Computes the agreement matrix for the different graph built by `run.py`.  
    takes as input the dataset id (`dopamine`, `richiardi`, or `full`) and saves the matrix a pickled pandas DataFrame in the folder `output`.
    - `shortest_path.py`: Computes the shortest path from the nodes containing samples of VGA and substantia nigra to the rest of the brain.  
    takes as input the dataset id (`dopamine`, `richiardi`, or `full`) and saves the information for each node in a pickled dictionary in the folder `output`.
    
## Citing
If you make use of this work in your research please cite the following [paper](https://www.biorxiv.org/content/10.1101/476382v1):

Patania, Alice, Pierluigi Selvaggi, Mattia Veronese, Ottavia DiPasquale, Paul Expert, and Giovanni Petri. "Topological gene-expression networks recapitulate brain anatomy and function." bioRxiv (2018): 476382.

### Bibtex

@article{patania2018topological,  
  title={Topological gene-expression networks recapitulate brain anatomy and function},  
  author={Patania, Alice and Selvaggi, Pierluigi and Veronese, Mattia and DiPasquale, Ottavia and Expert, Paul and Petri, Giovanni},  
  journal={bioRxiv},  
  pages={476382},  
  year={2018},  
  publisher={Cold Spring Harbor Laboratory}  
}

## How to reproduce the results:

```shell
$ python parameters.py name_gene_list
$ python selection.py name_gene_list
$ python run.py name_gene_list
$ python agreement matrix.py name_gene_list
```
with `name_gene_list` is one of (`dopamine`, `richiardi`, or `full`)

The file `shortest_path.py` can be run with any output from `run.py` and `selection.py`. In the paper we only looked at the outcomes from the dopamine related mappers, but it can be run on any other output.

### Dependencies
An up-to-date Python 3.5 distribution, with the standard packages provided by the anaconda distribution is required.  
In particular, the code was tested with:  
pandas (__version__) etc  
