library('Seurat')
library('dplyr')
load('Data/MouseCortex.RData')
Data <- MouseCortex@meta.data
time <- substr(Data$Time_points,2,3) %>% as.numeric()
Data <- data.frame(timepoint=time,celltype=as.numeric(Data$res.0.6),tSNE_1=MouseCortex@dr$tsne@cell.embeddings[,1],tSNE_2=MouseCortex@dr$tsne@cell.embeddings[,2])
rownames(Data) <- names(MouseCortex@ident)
CelltypeNames <- c('1-Neurons','2-Young Neurons','3-APs/RPs','4-IPs','5-APs/RPs','6-Young Neurons','7-IPs')
label <- Data$celltype
for(i in 1:length(label)){
  if(label[i]==1){label[i] <- CelltypeNames[1]}
  if(label[i]==2){label[i] <- CelltypeNames[2]}
  if(label[i]==3){label[i] <- CelltypeNames[3]}
  if(label[i]==4){label[i] <- CelltypeNames[4]}
  if(label[i]==5){label[i] <- CelltypeNames[5]}
  if(label[i]==6){label[i] <- CelltypeNames[6]}
  if(label[i]==7){label[i] <- CelltypeNames[7]}
}
Data$typeName <- label
save(Data,file='Data.RData')


