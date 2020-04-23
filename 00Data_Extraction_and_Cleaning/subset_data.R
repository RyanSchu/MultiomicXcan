library(argparse)
library(dplyr)
library(tidyr)
library(data.table)

parser <- ArgumentParser()
parser$add_argument("--superset", help="file containing square matrix with observations and all features")
parser$add_argument("--margin", help="a margin of 1 indicates that rows of the superset are individuals/observations with columns as features. Margin of 2 implies columns are obs/inds")
parser$add_argument("--subset", help="a text file containing n lines, where each line is the id of the observations to extract")
parser$add_argument("--features", help="a text file containing p lines, where each line is the id of the features to extract")
parser$add_argument("--sanitize", help="sanitize the ENSG ids",action="store_true")
parser$add_argument("--out",help="name of the file you would like to output as, optionally including full path")
args <- parser$parse_args()

"%&%" = function(a,b) paste(a,b,sep="")

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

extract.subset<-function(superset,margin,features,observations){
  if (margin == 1){
    
    reduced<-superset %>% 
      rename_at(colnames(superset)[1] , ~ "id") %>% 
      select(one_of(c("id",features))) %>% 
      filter(id %in% observations)
    return( reduced )
  } else {
    reduced<-superset %>% 
      rename_at(colnames(superset)[1] , ~ "id") %>% 
      select(one_of(c("id",observations))) %>% 
      filter(id %in% features)
    return( reduced )
  }
}

sanitize.ENSG<-function(df,margin){
  if (margin == 1){
    names<-gsub("\\.[0-9]+","",colnames(df))
    colnames(df)<-names
  } else {
    # str(df)
    df<-df %>%  rename_at(colnames(superset)[1] , ~ "id") %>% mutate(id=gsub("\\.[0-9]+","",id))
    return(df)
  }
}

superset<-check.fread(args$superset)
superset<-sanitize.ENSG(superset,margin=args$margin)

features<-check.fread(args$features) %>% unlist() %>% unname() %>% gsub("\\.[0-9]+","",.)


subset<-check.fread(args$subset) %>% unlist() %>% unname()


reduced<-extract.subset(superset = superset,features = features,observations = subset,margin = args$margin)

fwrite(reduced,args$out,sep='\t',col.names = T)






