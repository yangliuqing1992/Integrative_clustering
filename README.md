# Integrative_clustering
Research work for integrative clustering analysis

## Codes
`Codes` folder contains:
1. code for numerical example
    - `code_simulation.R` can run on a simulated data as an example.
    - Please refer `readme_threshold.txt` in `code for selecting threshold` folder first to generate thresholds.  
3. code for analyzing mouse embryo data 
    - `code_mouse_embryo.R` is for analysis at number of clusters `K=14`. It takes input `X.csv`, `A.csv`, `Xaff.csv` and output labels for each methods and print out BIC and CBIC for `K=14`.
    - All the input data sets are in folder `inputs`. `X_to_affinity.py` takes `X.csv` as input and calculate an affinity matrix `Xaff.csv`, which used to calculate the initial labels to the integrative method. 
    - `outputs` folder contains a excel file contains the labels from multiple methods together with corresponding gene names. Meawhile columns `ig`, `hc`, `pl_sc` are from integrative method, hierarchical clustering, SBM respectively.
