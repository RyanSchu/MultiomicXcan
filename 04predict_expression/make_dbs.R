library(dplyr)
library(RSQLite)
library(qvalue)
library(data.table)
"%&%" <- function(a,b) paste(a,b, sep='')

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
driver <- dbDriver('SQLite')

#gene_annot <- read.table("../../prepare_data/expression/gencode.v19.genes.patched_contigs.parsed.txt", header = T, stringsAsFactors = F)

model_summaries<-fread("/home/ryan/topmed/multiomic_modeling/output/CAU_PBMC_chr1_multiomic_models.txt_chr1_model_summaries.txt",header=T,stringsAsFactors = F,fill=TRUE)
model_summaries<-reformat_data(model_summaries)

tiss_summary <- read.table('/home/ryan/topmed/multiomic_modeling/output/CAU_PBMC_chr1_multiomic_models.txt_chr1_pop_chr_summary.txt', header = T, stringsAsFactors = F)
n_samples <- tiss_summary$n_samples

weights <- read.table('/home/ryan/topmed/multiomic_modeling/output/CAU_PBMC_chr1_multiomic_models.txt_chr1_weights.txt', header = T, stringsAsFactors = F)

chrom<-c(2:22)
for (c in chrom){

        print(c)
        # Extra table ----
        tmp_model_summaries<-fread("/home/ryan/topmed/multiomic_modeling/output/CAU_PBMC_chr" %&% c %&% 
                                 "_multiomic_models.txt_chr" %&% c %&% "_model_summaries.txt",header=T,stringsAsFactors = F,fill=TRUE)
        tmp_model_summaries<-reformat_data(tmp_model_summaries)
        model_summaries<-rbind.data.frame(model_summaries,tmp_model_summaries)
        tmp_weights <- read.table('/home/ryan/topmed/multiomic_modeling/output//CAU_PBMC_chr' %&% c %&% 
                                '_multiomic_models.txt_chr' %&% c %&% '_weights.txt', header = T, stringsAsFactors = F)
        weights<-rbind.data.frame(weights,tmp_weights)
        
}
        # Weights Table -----

str(model_summaries)
str(weights)
weights<-weights%>% mutate(rsid="chr" %&% rsid)
        pvalues<-model_summaries$zscore_pval
        qvalues<-tryCatch(qvalue(pvalues), error = function(cond) {message('Error'); message(geterrmessage()); list()})
        model_summaries <- rename(model_summaries, 
                                  gene = gene_id, 
                                  genename = gene_name,
                                  n.snps.in.window = n_snps_in_window,
                                  n.snps.in.model = n_snps_in_model,
                                  pred.perf.R2 = rho_avg_squared,
                                  pred.perf.pval = zscore_pval
                                  )# %>% mutate(pred.perf.qval = qvalues$qvalues)
        if (length(qvalues) == 0){
          model_summaries <- model_summaries %>% mutate(pred.perf.qval = 0)
        } else {
          model_summaries <- model_summaries %>% mutate(pred.perf.qval = qvalues$qvalues)
        }
        
        conn <- dbConnect(drv = driver, '/home/ryan/topmed/multiomic_modeling/dbs/CAU_multiomic_model.db')
        dbWriteTable(conn, 'extra', model_summaries, overwrite = TRUE)
        dbGetQuery(conn, "CREATE INDEX gene_model_summary ON extra (gene)")

        weights <- rename(weights,gene = gene_id, weight = beta, ref_allele = ref, eff_allele=alt)
        dbWriteTable(conn, 'weights', weights, overwrite = TRUE)
        dbGetQuery(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
        dbGetQuery(conn, "CREATE INDEX weights_gene ON weights (gene)")
        dbGetQuery(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
        
        # Sample_info Table ----
        sample_info <- data.frame(n_samples = n_samples, population ="CAU", tissue = "PBMC")
        dbWriteTable(conn, 'sample_info', sample_info, overwrite = TRUE)
        
        # Construction Table ----
        construction <- tiss_summary %>%
                          select(cv_seed)
        dbWriteTable(conn, 'construction', construction, overwrite = TRUE)
        
        ### write models with no NA
        # 
        # conn2 <- dbConnect(drv = driver, '/home/ryan/topmed/multiomic_modeling/dbs/CAU_multiomic_model.db')
        # complete_model_summaries<-model_summaries[complete.cases(model_summaries),]
        # dbWriteTable(conn2, 'extra', complete_model_summaries, overwrite = TRUE)
        # dbGetQuery(conn2, "CREATE INDEX gene_model_summary ON extra (gene)")
        # 
        # # Weights Table -----
        # weights <- read.table('/home/ryan/topmed/proteome/dapg_net_CORRECT_WINDOWS/' %&% window %&% 
        #                         '_windows/' %&% pop %&% '_PIP_' %&% pip %&% '_clus_' %&% clus %&% 
        #                         '_weights.txt', header = T, stringsAsFactors = F)
        # 
        # 
        # weights <- rename(weights, gene = gene_id, weight = beta, ref_allele = ref, eff_allele=alt)
        # dbWriteTable(conn2, 'weights', weights, overwrite = TRUE)
        # dbGetQuery(conn2, "CREATE INDEX weights_rsid ON weights (rsid)")
        # dbGetQuery(conn2, "CREATE INDEX weights_gene ON weights (gene)")
        # dbGetQuery(conn2, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
        # 
        # # Sample_info Table ----
        # sample_info <- data.frame(n_samples = n_samples, population = pop, tissue = "PBMC")
        # dbWriteTable(conn2, 'sample_info', sample_info, overwrite = TRUE)
        # 
        # # Construction Table ----
        # construction <- tiss_summary %>%
        #   select(cv_seed)
        # dbWriteTable(conn2, 'construction', construction, overwrite = TRUE)
        # 
        # ### Write only significant models
        # conn3 <- dbConnect(drv = driver, '/home/ryan/topmed/proteome/dapg_net_CORRECT_WINDOWS/dbs/' %&% pop %&% '_'%&% window %&% '_WG_PIP_' %&% pip %&% '_cluster_filt_' %&% clus %&% '_rho0.1_zpval0.05.db')
        # signif_models<-filter(complete_model_summaries, pred.perf.pval < 0.05 , rho_avg > 0.1)
        # 
        # dbWriteTable(conn3, 'extra', signif_models, overwrite = TRUE)
        # dbGetQuery(conn3, "CREATE INDEX gene_model_summary ON extra (gene)")
        # 
        # # Weights Table -----
        # weights <- read.table('/home/ryan/topmed/proteome/dapg_net_CORRECT_WINDOWS/' %&% window %&% 
        #                         '_windows/' %&% pop %&% '_PIP_' %&% pip %&% '_clus_' %&% clus %&% 
        #                         '_weights.txt', header = T, stringsAsFactors = F)
        # weights <- rename(weights, gene = gene_id, weight = beta, ref_allele = ref, eff_allele=alt)
        # weights<- weights %>% filter ( gene %in% signif_models$gene)
        # dbWriteTable(conn3, 'weights', weights, overwrite = TRUE)
        # dbGetQuery(conn3, "CREATE INDEX weights_rsid ON weights (rsid)")
        # dbGetQuery(conn3, "CREATE INDEX weights_gene ON weights (gene)")
        # dbGetQuery(conn3, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
        # 
        # 
        # # Sample_info Table ----
        # sample_info <- data.frame(n_samples = n_samples, population = pop, tissue = "PBMC")
        # dbWriteTable(conn3, 'sample_info', sample_info, overwrite = TRUE)
        # 
        # # Construction Table ----
        # construction <- tiss_summary %>%
        #   select( cv_seed)
        # dbWriteTable(conn3, 'construction', construction, overwrite = TRUE)
        # 
        # ### write table with no NAs and no empty models
        # 
        # # conn2 <- dbConnect(drv = driver, '/home/ryan/topmed/proteome/dapg_net/dbs/' %&% pop %&% '_WG_PIP_' %&% pip %&% '_' %&% clus %&% "_noNA__noEmptyModels_unfiltered.db")
        # # complete_model_summaries<-model_summaries[complete.cases(model_summaries),]
        # # dbWriteTable(conn2, 'extra', complete_model_summaries, overwrite = TRUE)
        # # dbGetQuery(conn2, "CREATE INDEX gene_model_summary ON extra (gene)")
        # # 
        # # 
        # # # Weights Table -----
        # # weights <- read.table('/home/ryan/topmed/proteome/dapg_net/PIP_' %&% pip %&% 
        # #                         '_' %&% clus %&% '/' %&% pop %&% 
        # #                         '/' %&% pop %&% '_chr1_PIP_' %&% pip %&% '_' %&% clus %&% '_chr1_weights.txt', header = T, stringsAsFactors = F)
        # # for (i in 2:22) {
        # #   weights <- rbind(weights,
        # #                    read.table('/home/ryan/topmed/proteome/dapg_net/PIP_' %&% pip %&% 
        # #                                 '_' %&% clus %&% '/' %&% pop %&% 
        # #                                 '/' %&% pop %&% '_chr' %&% as.character(i) %&% '_PIP_' %&% pip %&% '_' %&% clus %&% '_chr' %&% as.character(i) %&% '_weights.txt', header = T, stringsAsFactors = F))
        # # }
        # # weights <- rename(weights, gene = gene_id, weight = beta, ref_allele = ref, eff_allele=alt)
        # # dbWriteTable(conn2, 'weights', weights, overwrite = TRUE)
        # # dbGetQuery(conn2, "CREATE INDEX weights_rsid ON weights (rsid)")
        # # dbGetQuery(conn2, "CREATE INDEX weights_gene ON weights (gene)")
        # # dbGetQuery(conn2, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
        # # 
        # # # Sample_info Table ----
        # # sample_info <- data.frame(n_samples = n_samples, population = pop, tissue = "PBMC")
        # # dbWriteTable(conn2, 'sample_info', sample_info, overwrite = TRUE)
        # # 
        # # # Construction Table ----
        # # construction <- tiss_summary %>%
        # #   select(chrom, cv_seed) %>%
        # #   rename(chromosome = chrom)
        # # dbWriteTable(conn2, 'construction', construction, overwrite = TRUE)
