suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages((library(reshape2)))
suppressMessages(library(methods))
suppressMessages(library(doMC))
suppressMessages(library(doRNG))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))

"%&%" <- function(a,b) paste(a,b, sep = "")

get_filtered_snp_annot <- function(snp_annot_file_name) {
  snp_annot <- read.table(snp_annot_file_name, header = T, stringsAsFactors = F) %>%
    filter(!((ref == 'A' & alt == 'T') |
               (ref == 'T' & alt == 'A') |
               (ref == 'C' & alt == 'G') |
               (ref == 'G' & alt == 'C')) &
             !(is.na(snp_ID))) %>%
    distinct(varid, .keep_all = TRUE) 
  # str(snp_annot)
  snp_annot
}

get_maf_filtered_genotype <- function(genotype_file_name,  maf, samples) {
  gt_df <- read.table(genotype_file_name, header = T, stringsAsFactors = F) %>% distinct(snp_ID,.keep_all=T) %>% column_to_rownames(var="snp_ID")
  # str(samples);str(gt_df)
  gt_df <- gt_df[,samples] %>% t() %>% as.data.frame()
  effect_allele_freqs <- colMeans(gt_df) / 2
  gt_df <- gt_df[,which((effect_allele_freqs >= maf) | (effect_allele_freqs <= 1 - maf))]
  gt_df
}


##Value check means that I'm checking nothing
get_gene_annotation <- function(gene_annot_file_name, chrom, gene_types=c('protein_coding',"aptamer", 'pseudogene', 'lincRNA',"aptamer","VALUE")){
  gene_df <- read.table(gene_annot_file_name, header = TRUE, stringsAsFactors = FALSE) %>%
    filter((chr == chrom) & gene_type %in% gene_types) %>% mutate(gene_id=gsub("\\.[0-9]+","",gene_id))
  gene_df
}



get_gene_type <- function(gene_annot, gene) {
  filter(gene_annot, gene_id == gene)$gene_type
}

get_gene_expression <- function(gene_expression_file_name, gene_annot) { #row are obs, columns are features
  expr_df <- as.data.frame(t(read.table(gene_expression_file_name, header = T, stringsAsFactors = F, row.names = 1)))
  names<-gsub("\\.[0-9]+","",colnames(expr_df))
  colnames(expr_df)<-names
  expr_df <- expr_df %>% select(one_of(intersect(gene_annot$gene_id, colnames(expr_df)))) #%>% mutate(id=gsub("\\.[0-9]+","",id))
  expr_df
}

get_protein_expression<- function(protein_expression_file_name,expr_df) { #row are obs, columns are features
  protein_df <- as.data.frame(t(read.table(protein_expression_file_name, header = T, stringsAsFactors = F, row.names = 1)))
  names<-gsub("\\.[0-9]+","",colnames(protein_df))
  colnames(protein_df)<-names
  #protein_df <- protein_df %>% select(one_of(intersect(colnames(protein_df), colnames(expr_df))))
  protein_df
}

extract_gene_expression<-function(gene_df,gene_name){
  # str(gene_name)
  # print(colnames(gene_df))
  gene_df<-gene_df %>% select(contains(gene_name))
  gene_df
}

extract_protein_expression<-function(protein_df,gene_name){
  # str(gene_name)
  # print(colnames(protein_df))
  protein_df<-protein_df %>% select(contains(gene_name))
  protein_df
}

get_gene_coords <- function(gene_annot, gene) {
  row <- gene_annot[which(gene_annot$gene_id == gene),]
  c(row$start, row$end)
}

get_cis_genotype <- function(gt_df, snp_annot, coords, cis_window) {
  # str(snp_annot)
  snp_info <- snp_annot %>% filter((pos >= (coords[1] - cis_window) & !is.na(snp_ID)) & (pos <= (coords[2] + cis_window)))
  if (nrow(snp_info) == 0)
    return(NA)
  cis_gt <- gt_df %>% select(one_of(intersect(snp_info$snp_ID, colnames(gt_df))))
  column_labels <- colnames(cis_gt)
  row_labels <- rownames(cis_gt)
  #str(cis_gt)
  #str(row_labels)
  #str(column_labels)
  # Convert cis_gt to a matrix for glmnet
  cis_gt <- matrix(as.matrix(cis_gt), ncol=ncol(cis_gt)) # R is such a bad language
  colnames(cis_gt) <- column_labels
  rownames(cis_gt) <- row_labels
  cis_gt
}

get_covariates <- function(covariate_file_name, samples) {
  cov_df <- read.table(covariate_file_name, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  cov_df <- cov_df[,-1:-3] %>% as.data.frame()
  cov_df
}

generate_fold_ids <- function(n_samples, n_folds=10) {
  n <- ceiling(n_samples / n_folds)
  fold_ids <- rep(1:n_folds, n)
  sample(fold_ids[1:n_samples])
}

adjust_for_covariates <- function(expression_vec, cov_df) {
  combined_df <- cbind(expression_vec, cov_df)
  expr_resid <- summary(lm(expression_vec ~ ., data=combined_df))$residuals
  expr_resid <- scale(expr_resid, center = TRUE, scale = TRUE)
  expr_resid
}

calc_R2 <- function(y, y_pred) {
  # cat("str of obs:");print(str(y))
  # cat("str of pred:");print(str(y_pred));print(dim(y_pred)[2])
  y_pred<-matrix(y_pred,nrow=dim(y_pred)[1],ncol=dim(y_pred)[2])
  if( ncol(y) != ncol(y_pred)) return(NA)
  R2list<-list()
  for (i in 1:ncol(y))
  {
    tss <- sum(y[,i]**2)
    rss <- sum((y[,i] - y_pred[,i])**2)
    R2list[i]<-1 - rss/tss
  }
  # cat("str of R2 before assignment:");str(R2list)
  return(R2list)
}

calc_corr <- function(y_pred,y_test) {
  y_pred<-matrix(y_pred,nrow=dim(y_pred)[1],ncol=dim(y_pred)[2])
  if( ncol(y_test) != ncol(y_pred)) return(NA)
  corrlist<-list()
  # str(y_test);str(y_pred)
  for (i in 1:ncol(y_test))
  {
    corrlist[i]<-cor(y_pred[,i], y_test[,i])
  }
  # cat("str of corr before assignment:");str(corrlist)
  return(corrlist)
}

calc_corr_pval <- function(y_pred,y_test) {
  y_pred<-matrix(y_pred,nrow=dim(y_pred)[1],ncol=dim(y_pred)[2])
  if( ncol(y_test) != ncol(y_pred)) return(NA)
  corrplist<-list()
  for (i in 1:ncol(y_test))
  {
    corrplist[i]<-cor.test(y_pred[,i], y_test[,i])$p.value
  }
  # cat("str of corr p before assignment:");str(corrplist)
  return(corrplist)
}

combine_pvals<-function(pval_df,n_folds){
  plist<-list()
  for (i in 1:ncol(pval_df))
  {
    # Fisher's method for combining p-values: https://en.wikipedia.org/wiki/Fisher%27s_method
    plist[i]<-pchisq(-2 * sum(log(pval_df[,i])), 2*n_folds, lower.tail = F)
  }
  plist[ncol(pval_df)+1]<-pchisq(-2 * sum(log(pval_df)), 2*n_folds, lower.tail = F)
  # cat("str of corr p before assignment:");str(corrplist)
  return(plist)
}

nested_cv_elastic_net_perf <- function(x, y, n_samples, n_train_test_folds, n_k_folds, alpha) {
  # Gets performance estimates for k-fold cross-validated elastic-net models.
  # Splits data into n_train_test_folds disjoint folds, roughly equal in size,
  # and for each fold, calculates a n_k_folds cross-validated elastic net model. Lambda parameter is
  # cross validated. Then get performance measures for how the model predicts on the hold-out
  # fold. Get the coefficient of determination, R^2, and a p-value, where the null hypothesis
  # is there is no correlation between prediction and observed.
  #
  # The mean and standard deviation of R^2 over all folds is then reported, and the p-values
  # are combined using Fisher's method.
  R2_folds <- data.frame(y1=rep(0, n_train_test_folds),y2=rep(0, n_train_test_folds))
  corr_folds <- data.frame(y1=rep(0, n_train_test_folds),y2=rep(0, n_train_test_folds))
  zscore_folds <- data.frame(y1=rep(0, n_train_test_folds),y2=rep(0, n_train_test_folds))
  pval_folds <- data.frame(y1=rep(0, n_train_test_folds),y2=rep(0, n_train_test_folds))
  
  # Outer-loop split into training and test set.
  train_test_fold_ids <- generate_fold_ids(n_samples, n_folds=n_train_test_folds)
  for (test_fold in 1:n_train_test_folds) {
    train_idxs <- which(train_test_fold_ids != test_fold)
    test_idxs <- which(train_test_fold_ids == test_fold)
    
    # str(train_idxs)
    # str(test_idxs)    
    x_train <- x[train_idxs, ]
    y_train <- y[train_idxs,]
    x_test <- x[test_idxs, ]
    y_test <- y[test_idxs,]
    # str(y_train)
    # Inner-loop - split up training set for cross-validation to choose lambda.
    cv_fold_ids <- generate_fold_ids(nrow(y_train), n_k_folds)
    y_pred <- tryCatch({
      # Fit model with training data.
      fit <- cv.glmnet(x_train, y_train, nfolds = n_k_folds, alpha = alpha, type.measure='mse',family="mgaussian", foldid = cv_fold_ids)
      # Predict test data using model that had minimal mean-squared error in cross validation.
      predict(fit, x_test, s = 'lambda.min')},
      # if the elastic-net model did not converge, predict the mean of the y_train (same as all non-intercept coef=0)
      error = function(cond) do.call(rbind, replicate(nrow(y_test), colMeans(y_train), simplify=FALSE)))# rep(mean(y_train), nrow(y_test)))
    # print(str(fit))
    R2_folds[test_fold,] <- calc_R2(y_test, y_pred)
    # cat("R2 str is :");print(str(R2_folds))
    # Get p-value for correlation test between predicted y and actual y.
    # If there was no model, y_pred will have var=0, so cor.test will yield NA.
    # In that case, give a random number from uniform distribution, which is what would
    # usually happen under the null.
    corr_folds[test_fold,] <- ifelse(all(diag(var(y_pred)) != 0), calc_corr(y_pred, y_test), 0)
    zscore_folds[test_fold,] <- atanh(corr_folds[test_fold,])*sqrt(nrow(y_test) - 3) # Fisher transformation
    pval_folds[test_fold,] <- ifelse(all(diag(var(y_pred)) != 0), calc_corr_pval(y_pred, y_test), runif(ncol(y_test)))
    # cat("corr str is :");print(str(corr_folds))
    # cat("zscore str is :");print(str(zscore_folds))
    # cat("corr pval str is :");print(str(pval_folds))
  }
  # 
  str(corr_folds)
  str(pval_folds)
  R2_avg <- colMeans(R2_folds)
  R2_sd <- sqrt(diag(var(R2_folds)))
  rho_avg <- colMeans(corr_folds)
  rho_se <- sqrt(diag(var(corr_folds)))
  rho_avg_squared <- rho_avg**2
  # Stouffer's method for combining z scores.
  zscore_est <- colSums(zscore_folds) / sqrt(n_train_test_folds)
  zscore_pval <- 2*pnorm(abs(zscore_est), lower.tail = FALSE)
  # Fisher's method for combining p-values: https://en.wikipedia.org/wiki/Fisher%27s_method
  pval_est <- combine_pvals(pval_folds,n_train_test_folds)
  list(R2_avg=R2_avg, R2_sd=R2_sd, pval_est=pval_est, rho_avg=rho_avg, rho_se=rho_se, rho_zscore=zscore_est, rho_avg_squared=rho_avg_squared, zscore_pval=zscore_pval)
}

do_covariance <- function(gene_id, cis_gt, rsids, varIDs) {
  model_gt <- cis_gt[,rsids, drop=FALSE]
  colnames(model_gt) <- rsids
  geno_cov <- cov(model_gt)
  geno_cov[lower.tri(geno_cov)] <- NA
  cov_df <- melt(geno_cov, varnames = c("rsid1", "rsid2"), na.rm = TRUE) %>%
    mutate(gene=gene_id) %>%
    select(GENE=gene, RSID1=rsid1, RSID2=rsid2, VALUE=value) %>%
    arrange(GENE, RSID1, RSID2)
  cov_df
}

main <- function(snp_annot_file, gene_annot_file, genotype_file, expression_file, 
                 covariates_file, chrom, prefix, maf=0.01, n_folds=10, n_train_test_folds=5,
                 seed=NA, cis_window=1e6, alpha=0.5, null_testing=FALSE,protein_file) {
  cat("starting\n")
  gene_annot <- get_gene_annotation(gene_annot_file, chrom)
  # cat("here\n")
  expr_df <- get_gene_expression(expression_file, gene_annot)
  protein_df<-get_protein_expression(protein_file,expr_df)
  samples <- rownames(expr_df)
  # str(expr_df)
  n_samples <- length(samples)
  genes <- colnames(expr_df)
  n_genes <- length(genes)
  snp_annot <- get_filtered_snp_annot(snp_annot_file)
  gt_df <- get_maf_filtered_genotype(genotype_file, maf, samples)
  # covariates_df <- get_covariates(covariates_file, samples)
  # cat("here\n")
  # Set seed----
  seed<-1234
  seed <- ifelse(is.na(seed), sample(1:1000000, 1), seed)
  set.seed(seed)
  
  # Prepare output data----
  model_summary_file <- prefix %&% '_chr' %&% chrom %&% '_model_summaries.txt'
  model_summary_cols <- c('gene_id', 'gene_name', 'gene_type', 'alpha', 'n_snps_in_window', 'n_snps_in_model', 'lambda_min_mse',
                          'test_R2_avg', 'test_R2_sd', 'cv_R2_avg', 'cv_R2_sd', 'in_sample_R2',
                          'nested_cv_fisher_pval', 'rho_avg', 'rho_se', 'rho_zscore', 'rho_avg_squared', 'zscore_pval',
                          'cv_rho_avg', 'cv_rho_se', 'cv_rho_avg_squared', 'cv_zscore_est', 'cv_zscore_pval', 'cv_pval_est')
  write(model_summary_cols, file = model_summary_file, ncol = 24, sep = '\t')
  
  weights_file <- prefix %&% '_chr' %&% chrom %&% '_weights.txt'
  weights_col <- c('gene_id', 'rsid', 'varID', 'ref', 'alt', 'beta')
  write(weights_col, file = weights_file, ncol = 6, sep = '\t')
  
  tiss_chr_summ_f <- prefix %&% '_chr' %&% chrom %&% '_pop_chr_summary.txt'
  tiss_chr_summ_col <- c('n_samples', 'chrom', 'cv_seed', 'n_genes')
  tiss_chr_summ <- data.frame(n_samples, chrom, seed, n_genes)
  colnames(tiss_chr_summ) <- tiss_chr_summ_col
  write.table(tiss_chr_summ, file = tiss_chr_summ_f, quote = FALSE, row.names = FALSE, sep = '\t')
  
  covariance_file <-  prefix %&% '_chr' %&% chrom %&% '_covariances.txt'
  covariance_col <- c('GENE', 'RSID1', 'RSID2', 'VALUE')
  write(covariance_col, file = covariance_file, ncol = 4, sep = ' ')
  
  # Attempt to build model for each gene----
  for (i in 1:n_genes) {
    cat(i, "/", n_genes, "\n")
    gene <- unlist(genes[i])
    gene_name <- gene_annot$gene_name[gene_annot$gene_id == gene]
    gene_type <- get_gene_type(gene_annot, gene)
    coords <- get_gene_coords(gene_annot, gene)
    cis_gt <- get_cis_genotype(gt_df, snp_annot, coords, cis_window)
    str(gene)
    str(gene_type)
    str(coords)
    str(cis_gt)
    if (all(is.na(cis_gt))) {
      # No snps within window for gene.
      model_summary <- c(gene, gene_name, gene_type, alpha, 0, 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      write(model_summary, file = model_summary_file, append = TRUE, ncol = 24, sep = '\t')
      next
    }
    model_summary <- c(gene, gene_name, gene_type, alpha, ncol(cis_gt), 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
    protein_model_summary <- c(gene %&% "protein", gene_name, gene_type, alpha, ncol(cis_gt), 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
    
    if (ncol(cis_gt) >= 2) {
      
      expression_vec <- extract_gene_expression(expr_df,gene)
      # str(expression_vec)
      protein_vec<-extract_protein_expression(protein_df,gene)
      protein<-colnames(protein_vec)
      protein_model_summary <- c(protein, gene_name, gene_type, alpha, ncol(cis_gt), 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      
      joint_expression<-cbind.data.frame(expression_vec,protein_vec)
      joint_expression<-matrix(as.matrix(joint_expression), ncol=ncol(joint_expression))
      # str(joint_expression)
      #quit()
      perf_measures <- nested_cv_elastic_net_perf(cis_gt, joint_expression, n_samples, n_train_test_folds, n_folds, alpha)
      # str(perf_measures)
      #quit()
      R2_avg <- perf_measures$R2_avg
      R2_sd <- perf_measures$R2_sd
      pval_est <- perf_measures$pval_est
      rho_avg <- perf_measures$rho_avg
      rho_se <- perf_measures$rho_se
      rho_zscore <- perf_measures$rho_zscore
      rho_avg_squared <- perf_measures$rho_avg_squared
      zscore_pval <- perf_measures$zscore_pval
      # Fit on all data
      cat("here now\n")
      cv_fold_ids <- generate_fold_ids(nrow(joint_expression), n_folds)
      # str(cv_fold_ids)
      # str(joint_expression)
      cat("fitting model\n")
      fit <- tryCatch(cv.glmnet(cis_gt, joint_expression, nfolds = n_folds, alpha = 0.5, family="mgaussian", foldid = cv_fold_ids,type.measure='mse',  keep = TRUE), #
                      error = function(cond) {message('Error'); message(geterrmessage()); list()})
      # print(str(fit))
      cat("here\n")
      if (length(fit) > 0) {
        cv_R2_folds <- data.frame(y1=rep(0, n_folds),y2=rep(0, n_folds))
        cv_corr_folds <- data.frame(y1=rep(0, n_folds),y2=rep(0, n_folds))
        cv_zscore_folds <- data.frame(y1=rep(0, n_folds),y2=rep(0, n_folds))
        cv_pval_folds <- data.frame(y1=rep(0, n_folds),y2=rep(0, n_folds))
        best_lam_ind <- which.min(fit$cvm)
        for (j in 1:n_folds) {
          fold_idxs <- which(cv_fold_ids == j)
          adj_expr_fold_pred <- fit$fit.preval[fold_idxs,,best_lam_ind];# print(str(joint_expression));print(str(adj_expr_fold_pred))#str(fit$fit.preval[fold_idxs,,best_lam_ind]);str(adj_expr_fold_pred)
          cv_R2_folds[j] <- calc_R2(joint_expression[fold_idxs,], adj_expr_fold_pred);
          cv_corr_folds[j,] <- ifelse(all(diag(var(adj_expr_fold_pred))) != 0, calc_corr(adj_expr_fold_pred, joint_expression[fold_idxs,]), 0)
          cv_zscore_folds[j,] <- atanh(cv_corr_folds[j,])*sqrt(nrow(joint_expression[fold_idxs,]) - 3) # Fisher transformation
          cv_pval_folds[j,] <- ifelse(all(diag(var(adj_expr_fold_pred))) != 0, calc_corr_pval(adj_expr_fold_pred, joint_expression[fold_idxs,]), runif(1))
        }
        #all(diag(var(y_pred)) != 0)
        cv_R2_avg <- colMeans(cv_R2_folds)
        cv_R2_sd <- var(cv_R2_folds)
        adj_expr_pred <- predict(fit, as.matrix(cis_gt), s = 'lambda.min')
        training_R2 <- calc_R2(joint_expression, adj_expr_pred)
        
        cv_rho_avg <- colMeans(cv_corr_folds)
        cv_rho_se <- var(cv_corr_folds)
        cv_rho_avg_squared <- cv_rho_avg**2
        # Stouffer's method for combining z scores.
        cv_zscore_est <- sum(cv_zscore_folds) / sqrt(n_folds)
        cv_zscore_pval <- 2*pnorm(abs(cv_zscore_est), lower.tail = FALSE)
        cv_pval_est <- pchisq(-2 * sum(log(cv_pval_folds)), 2*n_folds, lower.tail = F)
        cat("here again\n")
        str(fit$nzero)
        str(fit$nzero[best_lam_ind])
        if (length(fit$glmnet.fit$beta[[2]][which(fit$glmnet.fit$beta[[2]][,best_lam_ind] != 0), best_lam_ind]) > 0) {
          cat("here a fifth time\n")
          weights <- fit$glmnet.fit$beta[[2]][which(fit$glmnet.fit$beta[[2]][,best_lam_ind] != 0), best_lam_ind]; str(weights)
          weighted_snps <- names(fit$glmnet.fit$beta[[2]][,best_lam_ind])[which(fit$glmnet.fit$beta[[2]][,best_lam_ind] != 0)]; #str(weighted_snps); #str(snp_annot)
          weighted_snps_info <- snp_annot %>% filter(snp_ID %in% weighted_snps) %>% select(snp_ID, varid, ref, alt)
          weighted_snps_info$protein <- protein
          weighted_snps_info <- weighted_snps_info %>%
            merge(data.frame(weights = weights, snp_ID=weighted_snps), by = 'snp_ID') %>%
            select(protein, snp_ID, varid, ref, alt, weights);#str(weighted_snps)
          write.table(weighted_snps_info, file = weights_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
          covariance_df <- do_covariance(protein, cis_gt, weighted_snps_info$snp_ID, weighted_snps_info$varid)
          write.table(covariance_df, file = covariance_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")
          protein_summary <- c(protein, gene_name, gene_type, alpha, ncol(cis_gt), fit$nzero[best_lam_ind], fit$lambda[best_lam_ind], unlist(R2_avg[2]), unlist(R2_sd[2]), unlist(cv_R2_avg[2]), unlist(cv_R2_sd[2]), unlist(training_R2[2]), unlist(pval_est[2]),
                               unlist(rho_avg[2]), unlist(rho_se[2]), unlist(rho_zscore[2]), unlist(rho_avg_squared[2]), unlist(zscore_pval[2]), unlist(cv_rho_avg[2]), unlist(cv_rho_se[2]), unlist(cv_rho_avg_squared[2]), unlist(cv_zscore_est[2]), unlist(cv_zscore_pval[2]), unlist(cv_pval_est[2])) 
          
          # print(str(fit$glmnet.fit$beta[[1]]))
          # # print(dim(fit$glmnet.fit$beta[[]]))
          # print(dim(fit$glmnet.fit$beta[[1]]))
          # print(length(fit$glmnet.fit$beta[[1]]))          
          # str(fit$glmnet.fit$beta[[1]]); #str(fit$glmnet.fit$beta[[2]])
          weights <- fit$glmnet.fit$beta[[1]][which(fit$glmnet.fit$beta[[1]][,best_lam_ind] != 0), best_lam_ind];#str(weights)#str(snp_annot)
          weighted_snps <- names(fit$glmnet.fit$beta[[1]][,best_lam_ind])[which(fit$glmnet.fit$beta[[1]][,best_lam_ind] != 0)]; #str(weighted_snps)#cat("here a third time\n")#str(weighted_snps); #str(snp_annot)
          weighted_snps_info <- snp_annot %>% filter(snp_ID %in% weighted_snps) %>% select(snp_ID, varid, ref, alt);#str(weighted_snps_info);
          #cat("here a third time\n")
          weighted_snps_info$gene <- gene;  cat("here a third time\n")
          weighted_snps_info <- weighted_snps_info %>%
            merge(data.frame(weights = weights, snp_ID=weighted_snps), by = 'snp_ID') %>%
            select(gene, snp_ID, varid, ref, alt, weights);#str(weighted_snps)
          write.table(weighted_snps_info, file = weights_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
          covariance_df <- do_covariance(gene, cis_gt, weighted_snps_info$snp_ID, weighted_snps_info$varid);  cat("here a fourth time\n")
          write.table(covariance_df, file = covariance_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")
          str(fit$lambda[best_lam_ind])
          model_summary <- c(gene, gene_name, gene_type, alpha, ncol(cis_gt), fit$nzero[best_lam_ind], fit$lambda[best_lam_ind], unlist(R2_avg[1]), unlist(R2_sd[1]), unlist(cv_R2_avg[1]), unlist(cv_R2_sd[1]), unlist(training_R2[1]), unlist(pval_est[1]),
                               unlist(rho_avg[1]), unlist(rho_se[1]), unlist(rho_zscore[1]), unlist(rho_avg_squared[1]), unlist(zscore_pval[1]), unlist(cv_rho_avg[1]), unlist(cv_rho_se[1]), unlist(cv_rho_avg_squared[1]), unlist(cv_zscore_est[1]), unlist(cv_zscore_pval[1]), unlist(cv_pval_est[1])) 
          # print(str(model_summary))
          
 } else {
          cat("length 0\n")
          protein_summary <- c(protein, gene_name, gene_type, alpha, ncol(cis_gt), 0, fit$lambda[best_lam_ind], unlist(R2_avg[2]), unlist(R2_sd[2]), unlist(cv_R2_avg[2]), unlist(cv_R2_sd[2]), unlist(training_R2[2]), unlist(pval_est[2]),
                        unlist(rho_avg[2]), unlist(rho_se[2]), unlist(rho_zscore[2]), unlist(rho_avg_squared[2]), unlist(zscore_pval[2]), unlist(cv_rho_avg[2]), unlist(cv_rho_se[2]), unlist(cv_rho_avg_squared[2]), unlist(cv_zscore_est[2]), unlist(cv_zscore_pval[2]), unlist(cv_pval_est[2])) 
          
          model_summary <- c(gene, gene_name, gene_type, alpha, ncol(cis_gt), 0, fit$lambda[best_lam_ind], unlist(R2_avg[1]), unlist(R2_sd[1]), unlist(cv_R2_avg[1]), unlist(cv_R2_sd[1]), unlist(training_R2[1]), unlist(pval_est[1]),
                        unlist(rho_avg[1]), unlist(rho_se[1]), unlist(rho_zscore[1]), unlist(rho_avg_squared[1]), unlist(zscore_pval[1]), unlist(cv_rho_avg[1]), unlist(cv_rho_se[1]), unlist(cv_rho_avg_squared[1]), unlist(cv_zscore_est[1]), unlist(cv_zscore_pval[1]), unlist(cv_pval_est[1])) 
   
          # model_summary <- c(gene, gene_name, gene_type, alpha, ncol(cis_gt), 0, fit$lambda[best_lam_ind], R2_avg, R2_sd,
          #                    cv_R2_avg, cv_R2_sd, training_R2, pval_est, rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval,
          #                    cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)
          # protein_summary <- c(protein, gene_name, gene_type, alpha, ncol(cis_gt), 0, fit$lambda[best_lam_ind], R2_avg, R2_sd,
          #                    cv_R2_avg, cv_R2_sd, training_R2, pval_est, rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval,
          #                    cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)
           
        }
      } else {
        cat("no snps found\n")
        # str(R2_avg[1])
        # str(R2_sd[1,1])
        model_summary <- c(gene, gene_name, gene_type, alpha, ncol(cis_gt), 0, NA, unlist(R2_avg[1]), R2_sd[1,1], NA, NA, NA, unlist(pval_est[1]), unlist(rho_avg[1]), rho_se[1,1], unlist(rho_zscore[1]), unlist(rho_avg_squared[1]), unlist(zscore_pval[1]),
                           NA, NA, NA, NA, NA, NA)
        protein_summary <- c(protein, gene_name, gene_type, alpha, ncol(cis_gt), 0, NA, unlist(R2_avg[2]), R2_sd[2,2], NA, NA, NA, unlist(pval_est[2]), unlist(rho_avg[2]), rho_se[2,2], unlist(rho_zscore[2]), unlist(rho_avg_squared[2]), unlist(zscore_pval[2]),
                             NA, NA, NA, NA, NA, NA)
      }
    }
    # str(model_summary);str(protein_summary)
    write(model_summary, file = model_summary_file, append = TRUE, ncol = 24, sep = '\t')
    write(protein_summary, file = model_summary_file, append = TRUE, ncol = 24, sep = '\t')
  }
}
