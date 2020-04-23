# install.packages("preprocessCore")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("preprocessCore")
library(data.table)
library(dplyr)
library(preprocessCore) #has normalize.quantiles function
library(tibble)
rankinv <- function(x){qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
'%&%' = function(a,b) paste(a,b,sep="")
gene_expression<-fread("/home/ryan/topmed/expression/TOPMed_MESA_RNAseq_Pilot_RSEMv1.3.0.rsem_genes_tpm_no_transcriptid.txt")
PBMC_list<-fread("/home/ryan/topmed/expression/observed_expression_age_PBMC_pop_code.txt")
Mono_list<-fread("/home/ryan/topmed/expression/observed_expression_age_Mono_pop_code.txt")
Tcell_list<-fread("/home/ryan/topmed/expression/observed_expression_age_Tcell_pop_code.txt")
adjmatlist <- list() #list to store adjmat's
adjdflist <- list() #list to store adjdf's
for (e in c(1,5)){
  for (p in c('CHN')){
    for (tiss in c("PBMC")){
      ids<-get(tiss %&% "_list") %>% filter(pop==p) %>% filter(Exam==e)
      if (dim(ids)[1]==0){
        next
      }
      df <- gene_expression %>% column_to_rownames("gene_id") %>% select(one_of(unlist(ids$V1))) #%>% t()
      df<-df %>%  t()  %>% as.data.frame() %>% rownames_to_column(var="TORID") %>% inner_join(ids,.,by=c("V1"="TORID"))
      #df<-cbind.data.frame(ids,df)
      rawmat<-as.matrix(df[,10:ncol(df)])
      str(rawmat)
      #logmat <- log(rawmat) #natural log transform
      #logmat[is.infinite(logmat)]<-0
      lmfunc <- function(x){resid(lm(x ~ df$age + df$gender_coide))} #get residuals of prot after adj for age & sex
      adjmat <- apply(rawmat, 2, lmfunc) #apply lmfunc to each column of logmat
      adjdf <- adjmat %>% aggregate(by=list(df$sidno),FUN=mean,drop=F)
      name <- tiss %&% p %&% e
      adjmatlist[[name]] <- adjmat
      adjdflist[[name]] <- adjdf
      print(name)
    }
  }
}
str(adjdflist)
elementMean <- function(my.list) { 
  arr <- array(unlist(my.list),c(dim(my.list[[1]])[1],dim(my.list[[1]])[2],length(my.list)))
  rowMeans( arr , na.rm=TRUE, dims = 2 )
}
# adjdflist$P
#full join to add NA's to df if missing Exam 1 or 5
# PBMCafadf <- full_join(adjdflist$PBMCAFA1, adjdflist$PBMCAFA5, by = "Group.1")
# PBMCcaudf <- full_join(adjdflist$PBMCCAU1, adjdflist$PBMCCAU5, by = "Group.1")
# PBMChisdf <- full_join(adjdflist$PBMCHIS1, adjdflist$PBMCHIS5, by = "Group.1")
PBMCchndf <- adjdflist$PBMCCHN1

# Monoafadf <- adjdflist$MonoAFA5
# Monocaudf <- adjdflist$MonoCAU5
# Monohisdf <- adjdflist$MonoHIS5
# 
# # Tcellafadf <- adjdflist$TcellAFA5
# Tcellcaudf <- adjdflist$TcellCAU5
# Tcellhisdf <- adjdflist$TcellHIS5


# for(pop in c('CAU','HIS')){
#   cat("taking mean of PBMC exams in ",pop,'\n')
#   df <- get("PBMC" %&% tolower(pop) %&% "df")
#   df1na <- select(df,ends_with(".x"))
#   df5na <- select(df,ends_with(".y"))
#   meanmat <- elementMean(list(df1na, df5na)) ##take the mean of exam 1 and exam 5
#   colnames(meanmat)
#   #invmeanmat <- round(apply(meanmat, 2, rankinv),6) ##rank-inverse normalize and round
#   finaldf <- cbind(df[,1],as.data.frame(meanmat)) ##add sidno to df
#   colnames(finaldf) <- c("sidno", unlist(gene_expression$gene_id)) ##retrieve column names
#   fwrite(finaldf, file = "/home/ryan/topmed/expression/PBMC/expression_" %&% pop %&% "age_sex_adj.txt",quote=F,sep="\t")
#   # pcs <- fread("/home/ryan/topmed/expression/PBMC/" %&% pop %&% "/genotypes/QC/PCA/unmerged_pca.eigenvec",header=TRUE)
#   # 
#   # pcmat <- as.matrix(pcs[,-1:-2])
#   # str(pcmat)
#   # finaldf<-finaldf %>% filter(sidno %in% pcs$IID) %>% column_to_rownames("sidno")
#   # str(finaldf)
#   # 
#   # lmfunc <- function(x){resid(lm(x ~ pcmat))}
#   # adjmat <- apply(finaldf, 2, lmfunc)
#   # fwrite(adjmat, file = "/home/ryan/topmed/expression/PBMC/expression_" %&% pop %&% "age_sex_adj_PC10.txt",quote=F,sep="\t")
# } 
colnames(PBMCchndf)<- c("sidno", unlist(gene_expression$gene_id))
fwrite(PBMCchndf, file = "/home/ryan/topmed/expression/PBMC/expression_CHNage_sex_adj.txt",quote=F,sep="\t")
# pcs <- fread("/home/ryan/TOPMED/expression/CHN/PCA/unmerged_pca.eigenvec",header=TRUE)
# pcmat <- as.matrix(pcs[,-1:-2])
# lmfunc <- function(x){resid(lm(x ~ pcmat))}
# adjmat <- apply(finaldf, 2, lmfunc)
# fwrite(adjmat, file = "/home/ryan/topmed/expression/PBMC/expression_CHNage_sex_adj_PC10.txt",quote=F,sep="\t")

# for( pop in c('AFA','CAU','HIS')){
#  for (tiss in c('Mono','Tcell')){
#    df <- get(tiss %&% tolower(pop) %&% "df")
#    fwrite(df, file = "/home/ryan/topmed/expression/" %&% tiss %&% "/expression_" %&% pop %&% "age_sex_adj.txt",quote=F,sep="\t")
#  }
# }
