library('corrplot')
library('ggplot2')
library('dplyr')
library('RColorBrewer')

Manual_Critical <- 1

Visualize_Phi <- function(Phi){
  Colors <- brewer.pal(n = length(Phi), name = "Set3")
  Colors[2] <- heat.colors(10)[8]
  
  SortPhi <- rep(0,N)
  j <- 1
  for (i in orderlabels) {
    SortPhi[j] <- Phi[i]
    j <- j+1
  }
  names(SortPhi) <- names(Phi)[orderlabels]
  barplot(SortPhi, legend.text = names(SortPhi),args.legend=locator(1),   main = expression(bold(Phi)),cex.main=4,horiz = FALSE,cex.axis=2,cex.names = 2.5,names.arg = GraphTypeOrder[orderlabels],col = Colors[GraphTypeOrder[orderlabels]])
}

Visualize_Phi_TSNE <- function(Phi){
  Data$Phi <- rep(0,dim(Data)[1])
  for(i in 1:dim(Data)[1]){
    Data$Phi[i] <- Phi[Data$celltype[i]]
  }
  Colors <- brewer.pal(n = N, name = "YlOrRd")
  par(mar=c(2,2,2,2))
  par(pin=c(3,3),mai=c(3,3,3,3))
  ggplot(Data,aes(x=tSNE_1,y=tSNE_2,colour=Phi))+
    geom_point(size=2,shape=16)+ggtitle(expression(bold(Phi))) +
    theme(plot.title = element_text(hjust = 0.5,size = 50,face = 'bold'))+
    scale_color_gradient(low = Colors[4],high = Colors[1])
}

Visualize_W <- function(W){
  SortW <- matrix(0,N,N)
  k <- 1
  l <- 1
  for (i in orderlabels) {
    l <- 1
    for (j in orderlabels) {
      SortW[k,l] <- W[i,j]
      l <- l+1
    }
    k <- k+1
  }
  colnames(SortW) <- GraphTypeOrder[orderlabels]
  rownames(SortW) <- GraphTypeOrder[orderlabels]
  corrplot(SortW,is.corr=FALSE,tl.srt=0,addCoef.col = "grey",main = '    W',cex.main=3,mar=c(2,2,2,2),cl.lim = c(min(SortW), max(SortW)),tl.cex = 1.8)
}

log_one <- function(a){
  label_0_a <- which(a==0)
  for(i in label_0_a){
    a[i] <- a[i]+addMin
  }
  new_loga <- log(a) %>% as.numeric()
  return(new_loga)
}

F <- function(Phi,W,p,beta){
  bigf <- (Phi + as.numeric(W %*% p) + beta*(log_one(p))) %>% as.numeric()
  return(bigf)
}

Solve_p <- function(Phi,W,init_p,init_t,t,Edgelist){
  p <- init_p
  critical <- (1/(length(p)-1))/Manual_Critical
  delta_t <- 1/Integral_step
  MM <- (t-init_t) %/% delta_t
  last_t <- (t-init_t) %% delta_t
  if(MM > 0){
    for (j in 1:MM) {
      count_time <- 0
      while((count_time)< delta_t){
        rest_time <- delta_t - count_time
        dv <- derive_dv(p,Phi,W,beta)
        new_dv <- (rest_time)*dv
        if(length(which(new_dv > critical))==0){
          count_time <- count_time + rest_time
          p <- p + (rest_time*derive_dp(p,Phi,W,beta))
        }else{
          new_step <- (rest_time/(max(dv)))*critical
          count_time <- count_time + new_step
          p <- p + (new_step*derive_dp(p,Phi,W,beta))
        }
      }
    }
  }
  count_time <- 0
  while((count_time)< last_t){
    rest_time <- last_t - count_time
    dv <- derive_dv(p,Phi,W,beta)
    new_dv <- (rest_time)*dv
    if(length(which(new_dv > critical))==0){
      count_time <- count_time + rest_time
      p <- p + (rest_time*derive_dp(p,Phi,W,beta))
    }else{
      new_step <- (rest_time/(max(dv)))*critical
      count_time <- count_time + new_step
      p <- p + (new_step*derive_dp(p,Phi,W,beta))
    }
  }
  return(p)
}

Plot_StochasticDynamic <- function(Phi,W){
  MM <- 200*(Timepoints[f] - Timepoints[1])
  Traj <- matrix(0,nrow = (MM+1),ncol=N)  ### each row is a probability distribution at time t
  ini_P <- ActualProbs[1,]
  Traj[1,] <- ini_P
  delta_tt <- (Timepoints[f]-Timepoints[1])/MM
  for (j in 1:MM) {
    init_t <- Timepoints[1] + (j-1)*delta_tt
    t <- Timepoints[1] + j*delta_tt
    next_P <- Solve_p(Phi,W,init_p = ini_P,init_t=init_t,t=t,Edgelist)
    Traj[j+1,] <- next_P
    ini_P <- next_P
  }
  time_xlab <- seq(Timepoints[1],Timepoints[f],delta_tt)+0.5
  par(mfrow=c(3,3))
  par(pin=c(1,2),mai=c(0.8,0.8,0.8,0.8))
  par(mar=c(5,5,5,5))
  Types <- names(Phi)
  for(j in 1:N){
    y_min <- min(min(Traj[,j]),min(as.numeric(ActualProbs[,j])))
    y_max <- max(max(Traj[,j]),max(as.numeric(ActualProbs[,j])))
    plot(time_xlab,y=Traj[,j],ylim=c(y_min-0.02,y_max+0.02),xlab='time',ylab='probability',main=Types[j],pch=20,cex.lab=1.5,cex.main=1.8,cex=0.2,cex.axis=1.2,type='o')
    points(Timepoints+0.5,ActualProbs[,j],pch=17,col=c('red','red','red','red'),cex=1)
  }
}

Plot_PotentialEnergy <- function(Phi,W){
  Colors <- brewer.pal(n = length(Phi), name = "Set3")
  Colors[2] <- heat.colors(10)[8]
  
  MM <- 200*(Timepoints[f] - Timepoints[1])
  Traj <- matrix(0,nrow = (MM+1),ncol=N) 
  PotentialEnergy <- matrix(0,nrow = (MM+1),ncol=N)  
  ini_P <- ActualProbs[1,]
  Traj[1,] <- ini_P
  PotentialEnergy[1,] <- F(Phi,W,ini_P,beta)
  delta_tt <- (Timepoints[f]-Timepoints[1])/MM
  for (j in 1:MM) {
    init_t <- Timepoints[1] + (j-1)*delta_tt
    t <- Timepoints[1] + j*delta_tt
    next_P <- Solve_p(Phi,W,init_p = ini_P,init_t=init_t,t=t,Edgelist)
    Traj[j+1,] <- next_P
    PotentialEnergy[j+1,] <- F(Phi,W,next_P,beta)
    ini_P <- next_P
  }
  time_xlab <- seq(Timepoints[1],Timepoints[f],delta_tt)
  par(mfrow=c(1,1))
  par(pin=c(1,1),mai=c(1,1,1,3))
  y_min <- min(PotentialEnergy)
  y_max <- max(PotentialEnergy)
  plot(time_xlab,PotentialEnergy[,1],ylim=c(y_min,y_max),xlab='Embryonic time',ylab='',main = expression(bold(Psi)),cex.main=4,pch=20,cex=0.2,type='o',col=Colors[GraphTypeOrder[1]],lwd=4,cex.lab=2,cex.axis=1.5)
  for(j in 2:N){
    points(time_xlab,PotentialEnergy[,j],pch=20,cex=0.2,col=Colors[GraphTypeOrder[j]])
    lines(time_xlab,PotentialEnergy[,j],col=Colors[GraphTypeOrder[j]],lwd=4)
  }
  xy <- par('usr')
  legend(x=xy[2]+xinch(0.4),y=xy[4]-yinch(1), legend=CelltypeNames,
         fill = Colors[GraphTypeOrder],xpd=TRUE,bty='n',cex=1.5)
}

Plot_PotentialEnergy_BetweenTwo <- function(Phi,W,TwoTypelabels,main){
  Colors <- brewer.pal(n = length(Phi), name = "Set3")
  Colors[2] <- heat.colors(10)[8]
  
  MM <- 200*(Timepoints[f] - Timepoints[1])
  Traj <- matrix(0,nrow = (MM+1),ncol=N) 
  PotentialEnergy <- matrix(0,nrow = (MM+1),ncol=N)  
  ini_P <- ActualProbs[1,]
  Traj[1,] <- ini_P
  PotentialEnergy[1,] <- F(Phi,W,ini_P,beta)
  delta_tt <- (Timepoints[f]-Timepoints[1])/MM
  for (j in 1:MM) {
    init_t <- Timepoints[1] + (j-1)*delta_tt
    t <- Timepoints[1] + j*delta_tt
    next_P <- Solve_p(Phi,W,init_p = ini_P,init_t=init_t,t=t,Edgelist)
    Traj[j+1,] <- next_P
    PotentialEnergy[j+1,] <- F(Phi,W,next_P,beta)
    ini_P <- next_P
  }
  time_xlab <- seq(Timepoints[1],Timepoints[f],delta_tt)
  par(mfrow=c(1,1))
  par(mai=c(1,1,1,1))
  i <- TwoTypelabels[1]
  j <- TwoTypelabels[2]
  Y1 <- PotentialEnergy[,i]
  Y2 <- PotentialEnergy[,j]
  y_min <- min(Y1,Y2)
  y_max <- max(Y1,Y2)
  plot(time_xlab,Y1,main=main,ylim=c(y_min-0.02,y_max+0.02),col=Colors[i],xlab='Embryonic time',ylab=expression(Psi),cex.main=4,pch=20,cex=0.2,type='o',lwd=4,cex.lab=2.4,cex.axis=1.5)
  lines(time_xlab,Y2,col=Colors[j],lwd=8)
  legend(locator(1),legend = CelltypeNames[c(i,j)],lwd=8,col=Colors[c(i,j)],cex=1.5)
}

Plot_FreeEnergy <- function(Phi,W,Times){
  freeEnergy <- function(p){
    EEE <- Phi %*% p + p %*% W %*%p + beta*p %*% log_one(p)
    return(EEE)
  }
  MM <- 200*(Times[2] - Times[1])
  FreeEnergy <- rep(0,MM+1)
  ini_P <- ActualProbs[1,]
  FreeEnergy[1] <- freeEnergy(ini_P)
  delta_tt <- (Times[2]-Times[1])/MM
  for (j in 1:MM) {
    init_t <- Times[1] + (j-1)*delta_tt
    t <- Times[1] + j*delta_tt
    next_P <- Solve_p(Phi,W,init_p = ini_P,init_t=init_t,t=t,Edgelist)
    FreeEnergy[j+1] <- freeEnergy(next_P)
    ini_P <- next_P
  }
  time_xlab <- seq(Times[1],Times[2],delta_tt)+0.5
  plot(time_xlab,FreeEnergy,xlab='Embryonic time',ylab='',main="Free Energy of System State",cex.main=2.4,pch=20,cex=0.2,type='o',lwd=4,cex.lab=2,cex.axis=1.5)
  lines(rep(17.5,length(seq(-0.22,0.1,by=0.01))),seq(-0.22,0.1,by=0.01),lwd=4,lty=2,col='red')
  text(locator(1),'time=17.5',cex=2,col='red')
}


DeriveEdgeflow <- function(Phi,W,Timepoints,P1){
    critical <- (1/(length(P1)-1))/Manual_Critical
  
    initial_t <- Timepoints[1]
    final_t <- Timepoints[length(Timepoints)]
    derive_FlowRate <- function(p){
      Fprob <- F(Phi,W,p,beta) %>% as.numeric()
      aaF <- (B %*% Fprob) %>% as.numeric()
      MatrixLabel <- list()
      gij <- c()
      for (l in 1:num_E) {
        vertex_two <- Edgelist[l,]
        if(aaF[l]>0){
          index_output <- which(CelltypeNames==vertex_two[2])
          gij[l] <- p[index_output]
        }else{
          index_output <- which(CelltypeNames==vertex_two[1])
          gij[l] <- p[index_output]
        }
      }
      flow <- (gij*aaF)
      return(flow)
    }
    p <- P1
    total_time <- final_t - initial_t
    delta_ti <- 1/Integral_step
    M <- (final_t - initial_t) %/% delta_ti
    last_ti <- (final_t-initial_t) %% delta_ti
    count_time <- 0
    Flow <- rep(0,num_E)
    FF <- matrix(0,ncol=num_E,nrow=1)
    PP <- matrix(0,ncol=N,nrow=1)
    for (j in 1:M) {
      count_time <- 0
      OldF <- Flow
      while((count_time)< delta_ti){
        rest_time <- delta_ti - count_time
        dv <- derive_dv(p,Phi,W,beta)
        new_dv <- (rest_time)*dv
        if(length(which(new_dv > critical))==0){
          count_time <- count_time + rest_time
          Flow <- Flow + (rest_time*derive_FlowRate(p))
          p <- p + (rest_time*derive_dp(p,Phi,W,beta))
        }else{
          new_step <- (rest_time/(max(dv)))*critical
          count_time <- count_time + new_step
          Flow <- Flow + (new_step*derive_FlowRate(p))
          p <- p + (new_step*derive_dp(p,Phi,W,beta))
        }
      }
      FF <- rbind(FF,Flow-OldF)
      PP <- rbind(PP,p)
    }
    if(last_ti > 0){
      count_time <- 0
      while((count_time)< last_ti){
        rest_time <- last_ti - count_time
        dv <- derive_dv(p,Phi,W,beta)
        new_dv <- (rest_time)*dv
        if(length(which(new_dv > critical))==0){
          count_time <- count_time + rest_time
          Flow <- Flow + (rest_time*derive_FlowRate(p))
          p <- p + (rest_time*derive_dp(p,Phi,W,beta))
        }else{
          new_step <- (rest_time/(max(dv)))*critical
          count_time <- count_time + new_step
          Flow <- Flow + (new_step*derive_FlowRate(p))
          p <- p + (new_step*derive_dp(p,Phi,W,beta))
        }
      }
    }
    FlowMatrix <- matrix(0,nrow=N,ncol=N)
    rownames(FlowMatrix) <- CelltypeNames
    colnames(FlowMatrix) <- CelltypeNames
    for(j in 1:num_E){
      FF <- Flow[j]
      Edgename <- Edgelist[j,]
      if(FF > 0){
        label1 <- which(CelltypeNames==Edgename[2])
        label2 <- which(CelltypeNames==Edgename[1])
        FlowMatrix[label1,label2] <- FF
      }else{
        label1 <- which(CelltypeNames==Edgename[1])
        label2 <- which(CelltypeNames==Edgename[2])
        FlowMatrix[label1,label2] <- abs(FF)
      }
    }
    RowSum <- apply(FlowMatrix,1,sum)
    diagValue <- c()
    for(kk in 1:N){
      diagValue <- c(diagValue,max(P1[kk]-RowSum[kk],0))
    }
    for(j in 1:N){
      FlowMatrix[j,j]<-diagValue[j]
    }
    return(FlowMatrix)
}

mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}


Plot_EdgeFlow7 <- function(Phi,W,Timepoints,labels,main){
  Colors <- brewer.pal(n = length(Phi), name = "Set3")
  Colors[2] <- heat.colors(10)[8]
  
  EdgeFlowS <- list()
  P1 <- ActualProbs[1,]
  EdgeFlowS[[1]] <- DeriveEdgeflow(Phi,W,Timepoints,P1)
  for(i in 2:4){
    Times <- Timepoints[(i-1):i]
    P1 <- ActualProbs[(i-1),]
    EdgeFlowS[[i]] <- DeriveEdgeflow(Phi,W,Times,P1)
  }
  netS <- list()
  par(mfrow=c(1,1))
  for(i in 2:f){
    EdgeFlow <- EdgeFlowS[[i]]
    for(j in 1:7){
      EdgeFlow[j,j] <- 0
    }
    netS[[i-1]] <- EdgeFlow
  }
  
  Net3 <- matrix(0,28,28)
  Net3[1:7,8:14] <- netS[[1]]
  Net3[8:14,15:21] <- netS[[2]]
  Net3[15:21,22:28] <- netS[[3]]
  row.names(Net3) <- rep(colnames(net1),f)
  colnames(Net3) <- row.names(Net3)
  
  vertex.weight <- c()
  for (i in 1:f) {
    weight <- ActualProbs[i,]
    weight <- weight*20
    vertex.weight <- c(vertex.weight,weight)
  }
  
  color.use3 <- rep(Colors, f)
  color.use3.frame <- rep(Colors, f)
  shape <- rep('circle',28)
  
  
  g <- graph_from_adjacency_matrix(Net3, mode = "directed", weighted = TRUE)
  edge.start <- ends(g, es=E(g), names=FALSE)
  coords <- matrix(NA, nrow(Net3), 2)
  coords[1:7,1] <- 0
  coords[8:14,1] <- 1.6
  coords[15:21,1] <- 3.2
  coords[22:28,1] <- 4.8
  
  
  coords[1:7,2] <- seq(1.2, 0, by = -1.2/(7-1))
  coords[8:14,2] <- seq(1.2, 0, by = -1.2/(7-1))
  coords[15:21,2] <- seq(1.2, 0, by = -1.2/(7-1))
  coords[22:28,2] <- seq(1.2, 0, by = -1.2/(7-1))
  
  coords_scale<-coords
  V(g)$size<-(vertex.weight)*1.5
  V(g)$color<-rep(Colors,f)
  V(g)$frame.color <- rep(Colors,f)
  V(g)$frame.font <- 5
  V(g)$label.color <- 'black'
  E(g)$width<- E(g)$weight/max(E(g)$weight)*12
  E(g)$arrow.width<-300*E(g)$width/max(E(g)$width)
  E(g)$arrow.size<-200*E(g)$width/max(E(g)$width)
  E(g)$label.color<-'black'
  E(g)$label.cex<-2
  E(g)$color<-adjustcolor(V(g)$color[edge.start[,1]],0.6)
  
  label.dist <- rep(1.6*4.5,28)
  label.locs <- c(rep(-pi, 7), rep(0, 7))
  text.pos <- cbind(seq(-1.6/1.5, 1.6/1.5,(3.2/1.5)/3), 1.5-1.5/7)
  add.vertex.shape("fcircle", clip=igraph.shape.noclip,plot=mycircle,parameters=list(vertex.frame.color=1, vertex.frame.width=1))
  par(mfrow=c(1,1))
  par(pin=c(1,1),mai=c(1,1,1,1))
  par(mar=c(4,4,4,4))
  plot(g,vertex.label=c(labels,rep(NA,14),labels),vertex.label.cex=1.2,edge.curved=0,layout=coords,margin=0.2,vertex.shape="fcircle", vertex.frame.width = rep(1,28),
       vertex.label.degree=label.locs, vertex.label.dist=label.dist, vertex.label.family="Helvetica",main = main,cex.main=10)
  text(text.pos, c("E11.5","E13.5",'E15.5','E17.5'), cex = 1.5, col = rep("#c51b7d",4))
}

