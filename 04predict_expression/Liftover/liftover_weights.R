library(dplyr)
library(tidyr)
library(ggplot2)
library(argparse)
library(data.table)

parser <- ArgumentParser()
parser$add_argument("--lifted", help="directory where all the hapmap files are written")
parser$add_argument("--weights", help="directory where all the hapmap files are written")
parser$add_argument("--out", help="directory where you would like to output your plots")
args <- parser$parse_args()
"%&%" = function(a,b) paste (a,b,sep="")


weights<-fread(args$weights,header = F)
colnames(weights)<-c("gene","id","varID","a1","a2","weight")
lifted<-fread(args$lifted,header = F)
colnames(lifted)<-c("chr","pos1","pos2","varID37","cpos","varID38")

final<-inner_join(weights,lifted,by=c("varID"="varID37")) %>% 
  select(gene,cpos,varID38,a1,a2,weight)
fwrite(final,args$out,col.names = F,sep="\t")