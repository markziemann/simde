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

for ( FRAC_DE in c(0.05,0.1,0.5)) {
  PDFNAME=paste(FRAC_DE,"_pw.pdf",sep="")
  pdf(file=PDFNAME,width=11.7,height=6.9)
  for (FC in c(0.584,1,1.584,2)) {
    par(mfrow=c(3,3))
    for (N_REPS in c(3,5,10)) {
      res=NULL
      for (DGE_FUNC in c("deseq2","edger","limma")) {
#      for (DGE_FUNC in c("edger","edger_ql","deseq","limma","absseq")) {
        for ( SUM_COUNT in c(10000000,40000000,100000000)) {
          for  ( VARIANCE in c(0,0.2,0.3,0.4,0.5)) {
            res_new<-agg_dge(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,gsets)
            res=rbind(res,res_new)
          }
        }
        res1<-res[which(res$DGE_FUNC==DGE_FUNC),]
        res1<-res1[which(res1$FRAC_DE==FRAC_DE),]

        #points
        res1_1e7<-res1[which(res1$SUM_COUNT==1E+7),]
        res1_1e7_v0<-res1_1e7[which(res1_1e7$VARIANCE==0),]
        res1_1e7_v2<-res1_1e7[which(res1_1e7$VARIANCE==0.2),]
        res1_1e7_v3<-res1_1e7[which(res1_1e7$VARIANCE==0.3),]
        res1_1e7_v4<-res1_1e7[which(res1_1e7$VARIANCE==0.4),]
        res1_1e7_v5<-res1_1e7[which(res1_1e7$VARIANCE==0.5),]

        res1_4e7<-res1[which(res1$SUM_COUNT==4E+7),]
        res1_4e7_v0<-res1_4e7[which(res1_4e7$VARIANCE==0),]
        res1_4e7_v2<-res1_4e7[which(res1_4e7$VARIANCE==0.2),]
        res1_4e7_v3<-res1_4e7[which(res1_4e7$VARIANCE==0.3),]
        res1_4e7_v4<-res1_4e7[which(res1_4e7$VARIANCE==0.4),]
        res1_4e7_v5<-res1_4e7[which(res1_4e7$VARIANCE==0.5),]

        res1_1e8<-res1[which(res1$SUM_COUNT==1E+8),]
        res1_1e8_v0<-res1_1e8[which(res1_1e8$VARIANCE==0),]
        res1_1e8_v2<-res1_1e8[which(res1_1e8$VARIANCE==0.2),]
        res1_1e8_v3<-res1_1e8[which(res1_1e8$VARIANCE==0.3),]
        res1_1e8_v4<-res1_1e8[which(res1_1e8$VARIANCE==0.4),]
        res1_1e8_v5<-res1_1e8[which(res1_1e8$VARIANCE==0.5),]

        par(mar=c(4,4,3,2))
        plot(res1_1e7_v0$r,res1_1e7_v0$p,xlab="recall",ylab="precision",pch=15,col="red",xlim=c(0,1),ylim=c(0,1))
        points(res1_1e7_v2$r,res1_1e7_v2$p,xlab="recall",ylab="precision",pch=16,col="red",xlim=c(0,1),ylim=c(0,1))
        points(res1_1e7_v3$r,res1_1e7_v3$p,xlab="recall",ylab="precision",pch=17,col="red",xlim=c(0,1),ylim=c(0,1))
        points(res1_1e7_v4$r,res1_1e7_v4$p,xlab="recall",ylab="precision",pch=0,col="red",xlim=c(0,1),ylim=c(0,1))
        points(res1_1e7_v5$r,res1_1e7_v5$p,xlab="recall",ylab="precision",pch=1,col="red",xlim=c(0,1),ylim=c(0,1))

        points(res1_4e7_v0$r,res1_4e7_v0$p,xlab="recall",ylab="precision",pch=15,col="blue",xlim=c(0,1),ylim=c(0,1))
        points(res1_4e7_v2$r,res1_4e7_v2$p,xlab="recall",ylab="precision",pch=16,col="blue",xlim=c(0,1),ylim=c(0,1))
        points(res1_4e7_v3$r,res1_4e7_v3$p,xlab="recall",ylab="precision",pch=17,col="blue",xlim=c(0,1),ylim=c(0,1))
        points(res1_4e7_v4$r,res1_4e7_v4$p,xlab="recall",ylab="precision",pch=0,col="blue",xlim=c(0,1),ylim=c(0,1))
        points(res1_4e7_v5$r,res1_4e7_v5$p,xlab="recall",ylab="precision",pch=1,col="blue",xlim=c(0,1),ylim=c(0,1))

        points(res1_1e8_v0$r,res1_1e8_v0$p,xlab="recall",ylab="precision",pch=15,col="dark gray",xlim=c(0,1),ylim=c(0,1))
        points(res1_1e8_v2$r,res1_1e8_v2$p,xlab="recall",ylab="precision",pch=16,col="dark gray",xlim=c(0,1),ylim=c(0,1))
        points(res1_1e8_v3$r,res1_1e8_v3$p,xlab="recall",ylab="precision",pch=17,col="dark gray",xlim=c(0,1),ylim=c(0,1))
        points(res1_1e8_v4$r,res1_1e8_v4$p,xlab="recall",ylab="precision",pch=0,col="dark gray",xlim=c(0,1),ylim=c(0,1))
        points(res1_1e8_v5$r,res1_1e8_v5$p,xlab="recall",ylab="precision",pch=1,col="dark gray",xlim=c(0,1),ylim=c(0,1))

        legend(0.4,0.3,legend=c("10M","40M","100M"),col=c("red", "blue","dark gray") ,pch=19,cex=0.6,title="read depth")
        legend(0.75,0.3,legend=c("0","0.2","0.3","0.4","0.5"),col=c("dark gray") ,pch=c(15:17,0,1),cex=0.6,title="added variance")
        mtext(paste(DGE_FUNC,N_REPS,"reps, 10% DEG") ,cex=0.8); grid()

      }
    }
  }
  dev.off()
}
write.table(res,file="simpw_res.tsv",quote=F,sep='\t')

