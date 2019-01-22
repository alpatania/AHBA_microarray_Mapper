# paper_NetNeuro_code
Code for reproducing results in the manuscript: 
"Topological gene-expression networks recapitulate brain anatomy and function" by Patania A. et al. (2018)

## Content:
1. data:
    - dataset normalized *this is not actually possible on git because of the file size, i'll find another way, in the meantime write me an email*
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
    - `shortest_path.py`: Computes the shortest path from the nodes containing samples of VGA and substantia nigra to the rest of the brain.  
    takes as input the dataset id (`dopamine`, `richiardi`, or `full`) and saves the information for each node in a pickled dictionary in the folder `output`.
    - `agreement_matrix.py`: Computes the agreement matrix for the different graph built by `run.py`.  
    takes as input the dataset id (`dopamine`, `richiardi`, or `full`) and saves the matrix a pickled pandas DataFrame in the folder `output`.
    
## To do list:
  - [x] put up the datasets
  - [x] make a parameters selection script
  - [x] make the shortest path script
  - [x] make all code into scripts that can be run from command line
  - [ ] write a tutorial on how to run all the code
  - [ ] find a way to put the dataset that are too big for git
  - [ ] make a script to compute the agreement matrix
