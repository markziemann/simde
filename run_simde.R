library("tidyverse")
library("parallel")
library("topconfects")
library("edgeR")
library("DESeq2")
library("limma")
library("ABSSeq")

source("simde_func.R")

#a is orig expression data
a<-read.table("https://raw.githubusercontent.com/markziemann/simde/master/start_data/ERR2539161/ERR2539161.se.tsv")

###############################################
#sanity checks to make sure the simulation is working as expected
###############################################
pdf(file="sanity_check.pdf")
x1<-simrna(a,5,10000000,0,0,0)
colnames(x1$x)<-paste("x1",colnames(x1$x),sep="_")
x2<-simrna(a,5,10000000,1,0,0)
colnames(x2$x)<-paste("x2",colnames(x2$x),sep="_")
xx<-merge(x1$x,x2$x,by=0)
heatmap(cor(xx[,2:ncol(xx)]),scale="none",cex.main=0.9,main="check1")

x1<-simrna(a,5,10000000,0,0,0)
colnames(x1$x)<-paste("x1",colnames(x1$x),sep="_")
x2<-simrna(a,5,10000000,0,0.1,1)
colnames(x2$x)<-paste("x2",colnames(x2$x),sep="_")
xx<-merge(x1$x,x2$x,by=0)
heatmap(cor(xx[,2:ncol(xx)]),scale="none",cex.main=0.9,main="check2")

x1<-simrna(a,5,10000000,0,0,0)
colnames(x1$x)<-paste("x1",colnames(x1$x),sep="_")
x2<-simrna(a,5,10000000,1,0.1,1)
colnames(x2$x)<-paste("x2",colnames(x2$x),sep="_")
xx<-merge(x1$x,x2$x,by=0)
heatmap(cor(xx[,2:ncol(xx)]),scale="none",cex.main=0.9,main="check3")
dev.off()


###############################################
# 10M reads with edger classic
###############################################
SIMS=10
for ( FRAC_DE in c(0.01,0.05,0.1,0.25)) {
  PDFNAME=paste(FRAC_DE,"_pr.pdf",sep="")
  pdf(file=PDFNAME,width=11.7,height=6.9)
  for (FC in c(0.584,1,1.584,2)) {
    par(mfrow=c(3,5))
    for (N_REPS in c(3,5,10)) {
      res=NULL
      for (DGE_FUNC in c("edger","edger_ql","deseq","limma","absseq")) {
        for ( SUM_COUNT in c(10000000,40000000,100000000)) {
          for  ( VARIANCE in c(0,0.2,0.3,0.4,0.5)) {
            res_new<-agg_dge(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC)
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

        legend(0.4,0.2,legend=c("10M","40M","100M"),col=c("red", "blue","dark gray") ,pch=19,cex=0.6,title="read depth")
        legend(0.75,0.2,legend=c("0","0.2","0.3","0.4","0.5"),col=c("dark gray") ,pch=c(15:17,0,1),cex=0.6,title="added variance")
        mtext(paste(DGE_FUNC,N_REPS,"reps, ",FRAC_DE,"DEG") ,cex=0.8); grid()

      }
    }
  }
  dev.off()
}

