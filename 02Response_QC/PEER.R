"%&%" = function(a,b) paste(a,b,sep="")
library(peer)
library(tibble)
library(dplyr)
source('/home/ryan/software/QTPipe/eQTL/rntransform.R')
source('/home/ryan/software/QTPipe/eQTL/ztransform.R')
library(preprocessCore)
library(argparse)
library(methods)
library(data.table)
parser <- ArgumentParser()
parser$add_argument("-e", "--expression", help="file path of the gene expression file")
parser$add_argument("-o", "--outputdir", help="output directory", type="character")
parser$add_argument("-t", "--tag", help="file tag for this run of samples")
parser$add_argument("--cov", help="file of covariates")
parser$add_argument("--pilot", help="if present, only tests one level of peer factors instead of iterating.",action="store_true")
parser$add_argument("--margin",help="1 for rows are individuals, 2 for columns are individuals")
args <- parser$parse_args()

check.fread<-function(file_name){
  zipped<-grepl(".gz$",file_name)
  if(zipped == T){
    file<-fread('zcat ' %&% file_name, header = T)
    return(file)
  } else{
    file<-fread(file_name, header = T)
    return(file)
  }
}

expr = check.fread(file=args$expression)

if ( args$margin==1 ){
	cat("treating rows as individuals\n")
	rows<- colnames(expr)[-1]
	cols<- unname(unlist(expr[,1]))
	expmat = as.matrix(expr[,-1]) #margin 1 keep as is
} else {
	cat("treating columns as individuals\n")
	rows<- unname(unlist(expr[,1]))
        cols<- colnames(expr)[-1]
        expmat = as.matrix(t(expr[,-1])) #margin 2 transpose, genes in cols, ids in rows
}
cat("Quantile Normalizing\n")
qn.expmat <- as.matrix(normalize.quantiles(expmat)) ##quantile normalize
qn<-as.matrix(t(qn.expmat))
rownames(qn) <- rows
colnames(qn) <- cols
write.table(file = args$outputdir %&% "/" %&% args$tag %&% "QN.txt", x = qn, quote = F, sep = '\t')

str(qn.expmat)

cat("Rank Normalizing\n")
rn <- apply(qn.expmat,1,"rntransform") ##rank transform to normality & transposes##
rn.qn.expmat<-as.matrix(t(rn))
rownames(rn) <- rows
colnames(rn) <- cols
write.table(file =  args$outputdir %&% "/" %&% args$tag %&% "RNQN.txt", x = rn, quote = F, sep = '\t')

if (!is.null(args$cov)){
  covs = as.matrix(read.table(args$cov,header=T, row.names = 1))
  covs<-as.matrix(as.numeric(covs))
 }
  ## i is number of peer factors to calculate, recommend 25% of sample size, but no more than 100
if ( args$pilot==T ) {
	pf_vec<-c(10)
	cat("Performing test analysis at PF 10\n")
} else {
	pf_vec<-c(10,15,20)
}

for(i in pf_vec){
    model = PEER()
    PEER_setPhenoMean(model, rn.qn.expmat)
    if (!is.null(args$cov)){
      PEER_setCovariates(model, as.matrix(covs))
    }  
    PEER_setNk(model,i)
    PEER_update(model)
    factors = PEER_getX(model)
    weights = PEER_getW(model)
    precision = PEER_getAlpha(model)
    residuals = PEER_getResiduals(model)
    print(dim(residuals))
    cat(args$expression, "Nk =", i)
    print(PEER_plotModel(model))
    adjexp <- t(residuals)
    print(dim(adjexp))
    print(nrow(adjexp))
#    str(cols)
#    str(rows)
    rownames(adjexp) <- rows
    colnames(adjexp) <- cols
    adjexp <- as.data.frame(adjexp) %>% rownames_to_column('gene_id')
    png(filename =args$outputdir %&% "/" %&% args$tag %&% "_" %&% i %&% '_peer.png')
    factors = PEER_getX(model)
    rownames(factors) <- cols
    write.table(factors,file=args$outputdir %&% "/" %&% args$tag %&% "_TOPMED_expression" %&% i %&% ".PEER.factors" %&% ".txt", quote=F,col.names=F)
    PEER_plotModel(model)
    dev.off()
    write.table(adjexp, file=args$outputdir %&% "/" %&% args$tag %&% "_TOPMED_expression" %&% i %&% '_peer_factor_adjusted.txt' , quote=F, row.names=F, sep='\t')
}


