setwd("/home/ryan/topmed/multiomic_modeling/modeling_scripts/")
source("multiomic_nested_cv_elnet.R")
"%&%" <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly = TRUE)
chrom <- argv[1]

snp_annot_file <- "/home/ryan/topmed/multiomic_modeling/predictor_matrices/preddb_input/chr" %&% chrom %&% "_anno.txt.gz"
gene_annot_file <- "/home/ryan/topmed/multiomic_modeling/response_matrices/gene_annotation/gene_annotation_hg19.txt"
genotype_file <- "/home/ryan/topmed/multiomic_modeling/predictor_matrices/preddb_input/chr" %&% chrom %&% "_geno.txt.gz"
expression_file <- "/home/ryan/topmed/multiomic_modeling/response_matrices/TRUE_intersect_CAU_PBMC_RNASeq_levels.txt"
protein_file <- "/home/ryan/topmed/multiomic_modeling/response_matrices/TRUE_intersect_CAU_PBMC_Protein_levels.txt"
prefix <- "/home/ryan/topmed/multiomic_modeling/output/CAU_PBMC_chr" %&% chrom %&% "_multiomic_models.txt"
#prefix <- "/home/ryan/topmed/proteome/dapg_net/redo_AFA/PIP_" %&% pip %&% "_" %&% suffix %&% "/" %&% pop %&% "/" %&% pop %&% "_chr" %&% chrom %&% "_" %&% "PIP_" %&% pip %&% "_" %&% suffix

main(snp_annot_file=snp_annot_file, 
	gene_annot_file, 
	genotype_file=genotype_file, 
	expression_file=expression_file, 
	chrom=as.numeric(chrom), 
	prefix=prefix, 
	null_testing=FALSE,
	protein_file=protein_file)


