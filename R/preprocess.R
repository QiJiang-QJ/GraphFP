library('igraph')
library('RColorBrewer')

Plot_DataWithLabels <- function(Data){
  Colors <- brewer.pal(n = N, name = "Set3")
  Colors[2] <- heat.colors(10)[8]
  
  label <- Data$celltype
  label <- as.character(label)
  centers <- matrix(0,7,2)
  for(i in 1:N){
    Clusters <- Data[which(label==i),]
    Clusters <- Clusters[,c(3:4)]
    centers[i,] <- apply(Clusters,2,mean)
  }
  for(i in 1:length(label)){
    if(label[i]==1){label[i] <- Colors[1]}
    if(label[i]==2){label[i] <- Colors[2]}
    if(label[i]==3){label[i] <- Colors[3]}
    if(label[i]==4){label[i] <- Colors[4]}
    if(label[i]==5){label[i] <- Colors[5]}
    if(label[i]==6){label[i] <- Colors[6]}
    if(label[i]==7){label[i] <- Colors[7]}
  }
  par(mfrow=c(1,1))
  par(pin=c(1,1),mai=c(1,1,1,1))
  plot(Data[,3],Data[,4],col=label,pch=20,cex=1,xlab = 'tSNE_1',ylab='tSNE_2',cex.lab=1.8)
  for(i in 1:7){
    x=centers[i,1]
    y=centers[i,2]
    text(x,y,CelltypeNames[i],font = 2,cex=1.3)
  }
}

get_incidentMatrix <- function(EEdgelist,Graph_type_order){
  rowname_of_EEdgelist <- c()
  N_EE <- dim(EEdgelist)[1]
  for(i_EE in 1:N_EE){
    rowname_of_EEdgelist <- c(rowname_of_EEdgelist,paste(EEdgelist[i_EE,1],EEdgelist[i_EE,2],sep = '-'))
  }
  EEdgelist <- as.data.frame(EEdgelist)
  rownames(EEdgelist) <- rowname_of_EEdgelist
  Edges_id <- 1:length(rowname_of_EEdgelist)
  EEdgelist$Edge_id <- Edges_id
  dataframeE <- EEdgelist
  incidenceMtrix <- table(dataframeE$Edge_id[row(dataframeE[-3])], unlist(dataframeE[-3]))
  colnam <- colnames(incidenceMtrix)
  newnew <- incidenceMtrix
  colnames(newnew) <- Graph_type_order
  for(i in 1:length(Graph_type_order)){
    loca <- which(colnam==Graph_type_order[i])
    newnew[,i]<- incidenceMtrix[,loca]
  }
  return(newnew)
}

get_incidentMatrix_D <- function(matrix){
  n_row <- dim(matrix)[1]
  for (i in 1:n_row) {
    a <- matrix[i,]
    location <- which(a==1)[1]
    matrix[i,location] <- -1
  }
  return(matrix)
}


Estimate_Probs <- function(data,orders,times){
  actualProbs <- matrix(0,nrow = length(times),ncol=length(orders))
  for (i in 1:dim(actualProbs)[1]) {
    Prob <- table(data[which(data$timepoint==times[i]),]$typeName)
    Prob <- Prob/sum(Prob)
    OrderProb <- Prob
    k <- 1
    for (j in orders) {
      OrderProb[k] <- Prob[j]
      k <- k+1
    }
    actualProbs[i,] <- OrderProb
  }
  colnames(actualProbs) <- orders
  rownames(actualProbs) <- times
  return(actualProbs)
}
