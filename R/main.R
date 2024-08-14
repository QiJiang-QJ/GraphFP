setwd('/Users/qj/Desktop/Github')
load('Data/Data.RData')
source('preprocess.R')
source('modeling.R')
source('visualization.R')


#CelltypeNames <- unique(Data$typeName)
CelltypeNames <- c('1-Neurons','2-Young Neurons','3-APs/RPs','4-IPs','5-APs/RPs','6-Young Neurons','7-IPs')
Timepoints <- sort(unique(Data$timepoint))
f <- length(Timepoints)
N <- length(CelltypeNames)
Plot_DataWithLabels(Data,N)


#### Define Transition Graph (Complete Graph)
Graph <- make_full_graph(N)
V(Graph)$name <- CelltypeNames
GraphTypeOrder <- as.numeric(substr(CelltypeNames,1,1))
num_E <- length(E(Graph))


### compute D & A & B
Edgelist <- get.edgelist(Graph)
D <- rep(1,num_E)
A <- get_incidentMatrix(Edgelist,V(Graph)$name)
B <- get_incidentMatrix_D(A)


### Estimate probablities at t1...tf
ActualProbs <- Estimate_Probs(Data,V(Graph)$name,Timepoints)


### Estimate parameters ( linear potential energy Phi & Interaction matrix W)
lambda=rep(1000,(f-1))
beta=1e-03

initial_Phi <- rep(0,N)
names(initial_Phi) <- CelltypeNames
initial_W <- matrix(0,N,N)
colnames(initial_W) <- CelltypeNames
rownames(initial_W) <- CelltypeNames
List_Result <- Estimate_THETA( lambda=lambda, beta=beta, initial_Phi=initial_Phi, initial_W=initial_W, Graph=Graph ,probs=ActualProbs )
### Estimate parameters using L1 regularization
# List_Result <- Estimate_THETA_regularization( lambda=lambda, beta=beta, initial_Phi=initial_Phi, initial_W=initial_W, Graph=Graph ,probs=ActualProbs )


### Analysis result and visualization
Phi <- List_Result$parameter$Phi
W <- List_Result$parameter$W
 
options(warn=-1)
orderlabels <- rev(order(Phi))
Visualize_Phi(Phi)
Visualize_Phi_TSNE(Phi)
Visualize_W(W)

Plot_StochasticDynamic(Phi,W)
Plot_PotentialEnergy(Phi,W)
Plot_FreeEnergy(Phi,W,c(Timepoints[1],50))

Plot_EdgeFlow7(Phi,W,Timepoints,labels=CelltypeNames,main='Probability Flow')
