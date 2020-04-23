library(dplyr)
library(data.table)
library(tibble)
'%&%' = function(a,b) paste(a,b,sep="")

pops<-c("CHN","CAU","HIS","AFA")
for ( pop in pops ){
  cat(pop,'\n')
  expression <- fread("/home/ryan/topmed/expression/PBMC/expression_" %&% pop %&% "age_sex_adj.txt",header=T)
  
  pcs <- fread("/home/ryan/topmed/expression/PBMC/" %&% pop %&% "/genotypes/QC/PCA/unmerged_pca.eigenvec",header=TRUE)
  pcmat <- as.matrix(pcs[,-1:-2])
  
  
  finaldf<-expression %>% filter(sidno %in% pcs$IID) 
  sidno<-finaldf %>% select("sidno")
  finaldf<-finaldf %>% column_to_rownames(var="sidno")
  # str(finaldf)
  # 
  str(sidno)
  lmfunc <- function(x){resid(lm(x ~ pcmat))}
  adjmat <- apply(finaldf, 2, lmfunc)
  adjmat<- cbind.data.frame(sidno,adjmat)
  str(adjmat)
  fwrite(adjmat, file = "/home/ryan/topmed/expression/PBMC/expression_" %&% pop %&% "age_sex_adj_PC10.txt",quote=F,sep="\t")
}