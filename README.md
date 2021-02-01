# GraphFP:a dynamic inference software for reconstructing the cell state-transition complex energy landscape from time series single-cell transcriptomic data

GraphFP can reconstruct cell state potential energies which measure their differentiation potencies, 
faithfully chart the probability flows between paired cell states over the dynamic processes of cell differentiation, 
and accurately quantify the stochastic dynamics of the cell type frequencies on probability simplex in continuous time.

The main part of GraphFP is performed in R. Before using GraphFP, please install deSolve,dplyr,igraph,RColorBrewer,corrplot,ggplot2 packages in advance. 

# Usage
GraphFP runs as follows: 

1. GraphFP takes processed scRNAseq data as input, which includes the measurement time of each cell and the cluster labels (cell type/state) of each cell. Clustering is not included in our software and we also assume that input cluster labels are already appropriate for our model. R script file "Preprocess.R" in folder Data shows that the information we need is extracted from the processed data by Tempora, which is reserved as input when training parameters.

2. Run the R script file, named "main.R" to estimate the parameters (linear potential energy Phi and Interaction matrix W) and further visualize the results by derived parameres. 



