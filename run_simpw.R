library("tidyverse")
library("parallel")
library("topconfects")
library("edgeR")
library("DESeq2")
library("limma")
library("ABSSeq")
library("stringi")

source("simpw_func.R")

# a is orig expression data
a<-read.table("https://raw.githubusercontent.com/markziemann/simde/master/start_data/ERR2539161/ERR2539161.se.tsv")

# GI is the gene information. It has gene name information for mapping accessions to gene symbols in GMT files
gi<-read.table("https://raw.githubusercontent.com/markziemann/simde/master/start_data/ERR2539161/GeneInfo.tsv",header=T,row.names=1)

# merge gene names
aa<-merge(a,gi,by=0)
aa<-aa[,-c(4:7)]
aa<-aggregate(. ~ GeneSymbol,aa,function(x) sum(as.numeric(as.character(x))))
aa$Row.names=NULL
rownames(aa)<-aa$GeneSymbol
aa$GeneSymbol=NULL

a<-aa[which(aa$ERR2539161>=10),,drop=F]

# generate some gene sets
gsets<-sapply( rep(100,1000) , function(x) {list(as.character(sample(rownames(a),x))) } )
#gsets<-sapply( rep(100,1000) , function(x) {list(as.character(sample(head(rownames(a),1000),x))) } )
names(gsets)<-stri_rand_strings(length(gsets), 15, pattern = "[A-Za-z]")

# Get the reactome gene set GMT
#download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip",destfile="ReactomePathways.gmt.zip")
#unzip("ReactomePathways.gmt.zip")
#GMT="ReactomePathways.gmt"

###############################################
# 10M reads with edger classic
###############################################
SIMS=10
unlink("simpw_res_running.tsv")
res=NULL
for ( FRAC_DE in c(0.2)) {
  PDFNAME=paste(FRAC_DE,"_pw.pdf",sep="")
  pdf(file=PDFNAME,width=11.7,height=6.9)
  for (FC in c(1)) {
    par(mfrow=c(3,3))
    for (N_REPS in c(3,5,10)) {
      for (DGE_FUNC in c("deseq2")) {
        for ( SUM_COUNT in c(10000000,40000000,100000000)) {
          for  ( VARIANCE in c(0,0.2,0.3,0.4,0.5)) {
            x<-agg_dge(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,gsets)
            write.table(x,file="simpw_res_running.tsv",quote=F,sep='\t',append=T)
            res=rbind(res,x)
          }
        }

      }
    }
  }
}

