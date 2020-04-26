##################################################################################################

#SETUP ENVIRONMENT
cat("SETUP ENVIRONMENT\n")
##################################################################################################

library(dplyr)
library(qvalue)
library(data.table)
library(ggplot2) 


##################################################################################################

#DEFINE FUNCTIONS
cat("DEFINE FUNCTIONS\n")
##################################################################################################


"%&%" = function(a,b) paste(a,b,sep="")

match_orders<-function(rowmatrix,colmatrix, idcol=2){
  sampleOrderColMat<-colmatrix[,idcol] %>% unlist() %>% unname()
  col1RowMat<-colnames(rowmatrix)[1]
  rowmatrix<-rowmatrix %>% select(one_of(c(col1RowMat,sampleOrderColMat)))
  return(rowmatrix)
}

get_gene_corr<-function(gene_name,colmatrix,rowmatrix,skipcol=1)
{
  #gene name is the name of the gene
  #colmatrix is a matrix containing genes as columns
  #rowmatrix is a matrix containing genes as rows
  #skipcol = The first n columns of the row matrix do not contain values, skip these for correlation test  

  predicted_exp<-colmatrix %>% 
    select(gene_name) %>% ##select the gene
    unlist() %>% unname()
  
  measured_exp<-rowmatrix %>% 
    filter( gene_id == gene_name) %>% 
    select((skipcol+1):ncol(rowmatrix)) %>%
    unlist() %>% 
    unname()
  
  correlation<-cor.test(measured_exp,predicted_exp, method = "spearman")
  
  expression_corr<-list()
  expression_corr[["gene_id"]]<-gene_name
  expression_corr[["estimate"]]<-correlation$estimate
  expression_corr[["p.value"]]<-correlation$p.value
  
  return(expression_corr)
}

shift_row<-function(vec){
  if (grepl("_",vec[2])){
    vec<-c(vec[-2],NA)
    return(vec)
  } else{
    return(vec)
  }
}

reformat_data<-function(df){
  names<-colnames(df)
  tmp<-apply(X=df, MARGIN=1, FUN=shift_row) %>% t() %>% as.data.frame()
  colnames(tmp)<-names
  tmp<-tmp %>% mutate_all(as.character) %>% mutate_at(names[4:length(names)],as.numeric)
  return(tmp)
}
##################################################################################################

#SET GLOBAL VARIABLES
cat("SET GLOBAL VARIABLES\n")
##################################################################################################

test<-T

pop_list<-c("AFA","CHN","HIS")

R2_filtering<-c(-1)


pi1_matrix<-matrix(NA, nrow = length(pop_list), ncol = 1)
rownames(pi1_matrix)<-pop_list
colnames(pi1_matrix)<-"CAU"

##################################################################################################

#READ IN & PROCESS DATA
cat("READ IN & PROCESS DATA\n")
##################################################################################################

for (pop in pop_list)
  {#per each set of observed expression data
  
  observed_expression<-fread("zcat /home/ryan/topmed/proteome/proteome_with_multigene_expression/Proteome_TOPMed_" %&% pop %&% "_ln_adjAgeSex_mean_rank-inverse_adj10PCs_MEQTL_reformatted_ids.txt.gz", header = T, sep = '\t',stringsAsFactors = F) %>%
  rename_at(vars(1),function(x){return("gene_id")}) %>% mutate(gene_id=gsub("\\.[0-9]+","",gene_id))
  #How well does each model replicate
    
  models<-fread("/home/ashley/TopMed/EN_Base_Models/CAU/en_baseline_CAU_chr1_model_summaries.txt",header=T,stringsAsFactors = F,fill=TRUE)
  models<-reformat_data(models)
  
  chrom<-c(2:22)
  for (c in chrom){
    
    # print(c)
    # Extra table ----
    tmp_models<-fread("/home/ashley/TopMed/EN_Base_Models/CAU/en_baseline_CAU_chr" %&% c %&% "_model_summaries.txt",header=T,stringsAsFactors = F,fill=TRUE)
    tmp_models<-reformat_data(tmp_models)
    models<-rbind.data.frame(models,tmp_models)
    # tmp_weights <- read.table('/home/ryan/topmed/multiomic_modeling/output//CAU_PBMC_chr' %&% c %&% 
    #                         '_multiomic_models.txt_chr' %&% c %&% '_weights.txt', header = T, stringsAsFactors = F)
    # weights<-rbind.data.frame(weights,tmp_weights)
    # 
  }
    models<-models %>% select(gene_id, test_R2_avg) %>% mutate(gene_id=gsub("\\.[0-9]+","",gene_id))

    predicted_expression<-fread("/home/ryan/topmed/multiomic_modeling/imputed_results/CAU_baseline_model_" %&% pop %&% "_geno_predicted_expression.txt",header=T,stringsAsFactors = F)
    
    obsGenes<-unlist(observed_expression$gene_id)
    predGenes<-colnames(predicted_expression) %>% gsub("\\.[0-9]+","",.)
    colnames(predicted_expression)<-predGenes
    
    if(test == T) {
      print(str(obsGenes))
      print(str(predGenes))
    }
    
    gene_list<-data.frame(gene_id=intersect(obsGenes,predGenes),stringsAsFactors = F)

    models<-inner_join(gene_list,models, by = "gene_id")
    cat("Existing Gene Models:\n")
    str(models)

    ##################################################################################################
    
    #ANALYZE DATA
    cat("ANALYZE DATA\n")
    ##################################################################################################
    
    for(j in R2_filtering)
      {#Calculate pi1 at different model R2 thresholds
      
      cat("testing replication rate at R2 threshold of ",j,"\n")
      filtered_gene_R2<-models %>% 
        filter( test_R2_avg > j)
      
      if(test==T) print(str(filtered_gene_R2))
      
      if(dim(filtered_gene_R2)[1] == 0)
      {#check if there are any models that meet R2 threshold
        
        cat("WARNING: ",pi1column," had no models at R2 threshold ", j, " for genes in ",obspop,"\n")
#        cat("WARNING: possible at higher R2 thresholds, but unlikely at low thresholds unless there is a gene name mismatch. Check that your genenames are the same between files")
        cat("WARNING: ASSIGNING 0 as pi1\n")
        pi1_matrix[obspop,pi1column] <- 0
        next
        
      }
      # if(test==T) quit()
      
      predictive_correlations<-sapply(X=gene_list$gene_id,FUN=get_gene_corr,colmatrix=predicted_expression,rowmatrix=observed_expression,simplify=T,USE.NAMES = T)
      predictive_correlations<-data.frame(gene_id=unlist(predictive_correlations[1,]),estimate=unlist(predictive_correlations[2,]),p.value=unlist(predictive_correlations[3,]))
      
      str(predictive_correlations)
      
      filtered_gene_R2<- filtered_gene_R2 %>% inner_join(predictive_correlations, by = "gene_id")
      
      ##################################################################################################
      
      #WRITE OUT CORRELATION DATA
      cat("WRITE OUT CORRELATION DATA\n")
      ##################################################################################################
      
      fwrite(filtered_gene_R2,"/home/ryan/topmed/multiomic_modeling/imputed_results/correlations/" %&% pop %&% "_spearman_correlation_with_baseline_model.txt",col.names = T,sep='\t')
      corr_pvals<-filtered_gene_R2$p.value
      qobjCorr <- tryCatch({qvalue(p = corr_pvals)},
                           error=function(cond){
                             cat("Error: ", pi1column, "and ",obspop, " caused error to qvalue, assigning 0 as pi1\n")
                             cond
                           })
      if(inherits(qobjCorr, "error")) {
        pi1_matrix[obspop,pi1column] <- 0
        next
      }
      pi1<- 1 - qobjCorr$pi0
      pi12<- signif(pi1,4)
      pi1_matrix[pop,1] <- pi12
    }#close
}

##################################################################################################

#WRITE OUT SUMMARY DATA
cat("WRITE OUT SUMMARY DATA")
##################################################################################################


pi1_matrix<-data.table(pi1_matrix)
row.names(pi1_matrix)<-pop_list
fwrite(pi1_matrix, "/home/ryan/topmed/multiomic_modeling/imputed_results/correlations/baseline_pi1_R2-1.txt", sep = '\t', col.names=T,row.names = T)