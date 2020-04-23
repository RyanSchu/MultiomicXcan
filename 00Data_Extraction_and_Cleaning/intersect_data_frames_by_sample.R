library(argparse)
library(dplyr)
library(tidyr)
library(data.table)

parser <- ArgumentParser()
parser$add_argument("--df1", help="file containing square matrix with observations and all features")
parser$add_argument("--df2", help="file containing square matrix with observations and all features")
parser$add_argument("--out1",help="directory to output as")
parser$add_argument("--out2",help="directory to output as")
args <- parser$parse_args()

"%&%" = function(a,b) paste(a,b,sep="")

check.fread<-function(file_name,select=NULL,header=T){
  zipped<-grepl(".gz$",file_name)
  if(zipped == T){
    file<-fread('zcat ' %&% file_name, header = header,select=select)
    return(file)
  } else{
    file<-fread(file_name, header = header,select=select)
    return(file)
  }
}

#?read.table

read.samples<-function(filname){
  samples<-read.table(filname,nrows=1,header=F,stringsAsFactors = F) %>% unlist() %>% unname()
  return(samples)
}


samples_1<-read.samples(args$df1)
samples_2<-read.samples(args$df2)

str(samples_1)
str(samples_2)


intersect_samples<-base::intersect(samples_1,samples_2)

df1<-check.fread(args$df1) %>% select(one_of(intersect_samples))
df2<-check.fread(args$df2) %>% select(one_of(intersect_samples))

fwrite(df1,args$out1,col.names = T,sep="\t")
fwrite(df2,args$out2,col.names = T,sep="\t")