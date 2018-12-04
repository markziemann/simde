library("limma")
library("edgeR")
library("DESeq2")
library("tidyverse")
library("parallel")

#a is orig expression data
a<-read.table("start_data/ERR2539161/ERR2539161.se.tsv")

########################################
# simulate some gene exression data
########################################
#inputs=(a,N_REPS,SUM_COUNT,
simrna<-function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC) {
#N_REPS=5 ; SUM_COUNT=10000000 ; VARIANCE=0.1 ; FRAC_DE=0.1 ; FC=3
df = NULL
for (k in paste0("data",1:(N_REPS*2)))  {
        b<-thinCounts(a,target.size=SUM_COUNT)
        colnames(b)=k
        df = cbind(df,b)
     }
#Number of differential genes
NDIF=round(nrow(df)*FRAC_DE)

if (VARIANCE>0) {
  #create some random values centred around 1 with some% error
  rand<-matrix(log2(rnorm(nrow(a)*N_REPS*2 , 2, VARIANCE)),ncol=N_REPS*2)
  #incorporate the noise
  df<-round(df*rand)
  #set any negative counts to zero
  df<-apply(df, 2, function(x) {ifelse(x < 0, 0, x)})
} 

if (NDIF>0) {
  message("prep fold changes")
  #Make even
  if ( NDIF%%2==1 ) { print("odd") ; NDIF=NDIF-1 }
  #sample DE genes
  DE_LIST<-sample(rownames(df) , NDIF)
  #split list in half
  UP_DE<-sample(DE_LIST,length(DE_LIST)/2)
  DN_DE<-setdiff(DE_LIST,UP_DE)
  #reformat as df and add fold change
  UP_DE<-as.data.frame(UP_DE)
  UP_DE$V1<-2^FC
  colnames(UP_DE)=c("Gene","FC")
  DN_DE<-as.data.frame(DN_DE)
  DN_DE$V1<-2^-FC
  colnames(DN_DE)=c("Gene","FC")
  ALL_DE<-rbind(DN_DE,UP_DE)
  #Go back to list for downstream work
  UP_DE<-UP_DE$Gene
  DN_DE<-DN_DE$Gene
  NON_DE<-as.data.frame(setdiff(rownames(df),ALL_DE$Gene))
  colnames(NON_DE)="Gene"
  NON_DE$FC=1
  ALL_DE<-rbind(ALL_DE,NON_DE)
  ALL_DE<-ALL_DE[ order(as.vector(ALL_DE$Gene)) , ]
  message("incorporate changes")
  df <- df[ order(row.names(df)), ]
  df2<-cbind(df,ALL_DE)
  df2$Gene=NULL
} else {
  df2<-as.data.frame( df )
  df2$FC<- 1
  UP_DE=NULL
  DN_DE=NULL
}
ODD_COLS=(1:(ncol(df2)-1))[c(TRUE,FALSE)]
EVEN_COLS=(1:(ncol(df2)-1))[c(FALSE,TRUE)]
controls<-df2[,ODD_COLS]
colnames(controls)=paste0( "ctrl_" ,1:ncol(controls) )
treatments<-round(df2[,EVEN_COLS]*df2$FC)
colnames(treatments)=paste0( "trt_" ,1:ncol(treatments) )
x<-cbind(controls,treatments)
rownames(x)=rownames(df2)
#filter out genes that are not expressed
x<- x[which(rowSums(x)/ncol(x)>10),]
UP_DE<-intersect(UP_DE,rownames(x))
DN_DE<-intersect(DN_DE,rownames(x))
xx <- list("x" = x, "UP_DE" = UP_DE, "DN_DE" = DN_DE )
xx
}
#simrna(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC)
#working examples
#xx<-simrna(a,5,10000000,0.2,20)  
#xxx<-lapply(list(a,a,a), simrna, 5, 10000000, 0.2, 20)
#xxx<-replicate(2,simrna(a,5,10000000,0.2,20))

#Thanks Gray Calhoun gcalhoun@iastate.edu for the following function
RepParallel <- function(n, expr, simplify = "array",...) {
      answer <-
        mclapply(integer(n), eval.parent(substitute(function(...) expr)),...)
      if (!identical(simplify, FALSE) && length(answer)) 
        return(simplify2array(answer, higher = (simplify == "array")))
      else return(answer)
    }
#xxx<-RepParallel(10,simrna(a,5,10000000,0.2,20), simplify=F, mc.cores = detectCores() )

#################################################
# define edgeR classic function
##################################################
edger<-function(y) {
res=NULL
label="simulate"
samplesheet<-as.data.frame(colnames(y))
colnames(samplesheet)="sample"
samplesheet$trt<-as.numeric(grepl("trt",colnames(y)))
design<-model.matrix(~samplesheet$trt)
rownames(design)=samplesheet$sample
y<-y[which(rowSums(y)/ncol(y)>=(10)),]
z<-DGEList(counts=y)
z<-calcNormFactors(z)
z<-estimateDisp(z, design,robust=TRUE,prior.df=1)
fit<-glmFit(z, design)
lrt<-glmLRT(fit)
dge<-as.data.frame(topTags(lrt,n=Inf))
dge$dispersion<-lrt$dispersion
dge<-merge(dge,lrt$fitted.values,by='row.names')
rownames(dge)=dge$Row.names
dge$Row.names=NULL
dge<-merge(dge,z$counts,by='row.names')
dge<-dge[order(dge$PValue),]
dge2<-subset(dge,FDR<0.05)

up_de<-dge2[which(dge2$logFC>1),1]
dn_de<-dge2[which(dge2$logFC<1),1]

res <- list("dge" = dge, "up_de" = up_de, "dn_de" = dn_de)
res
}

#################################################
# define edgeR QL function
##################################################
edger_ql<-function(y) {
res=NULL
label="simulate"
samplesheet<-as.data.frame(colnames(y))
colnames(samplesheet)="sample"
samplesheet$trt<-as.numeric(grepl("trt",colnames(y)))
design<-model.matrix(~samplesheet$trt)
rownames(design)=samplesheet$sample
y<-y[which(rowSums(y)/ncol(y)>=(10)),]
z<-DGEList(counts=y)
z<-calcNormFactors(z)
z<-estimateDisp(z, design,robust=TRUE,prior.df=1)
fit<-glmQLFit(z, design)
lrt<-glmQLFTest(fit)
dge<-as.data.frame(topTags(lrt,n=Inf))
dge$dispersion<-lrt$dispersion
dge<-merge(dge,lrt$fitted.values,by='row.names')
rownames(dge)=dge$Row.names
dge$Row.names=NULL
dge<-merge(dge,z$counts,by='row.names')
dge<-dge[order(dge$PValue),]
dge2<-subset(dge,FDR<0.05)

up_de<-dge2[which(dge2$logFC>1),1]
dn_de<-dge2[which(dge2$logFC<1),1]

res <- list("dge" = dge, "up_de" = up_de, "dn_de" = dn_de)
res
}

#################################################
# define DESeq2 function
##################################################
deseq<-function(y) {
res=NULL
label="simulate"
samplesheet<-as.data.frame(colnames(y))
colnames(samplesheet)="sample"
samplesheet$trt<-factor(as.numeric(grepl("trt",colnames(y))))
dds <- DESeqDataSetFromMatrix(countData = y, colData = samplesheet, design = ~ trt )
res <- DESeq(dds)
z<- results(res)
vsd <- vst(dds, blind=FALSE)
zz<-cbind(z,assay(vsd))
zz<-as.data.frame(zz[order(zz$padj),])
dge2<-subset(zz,padj<0.05)
up_de<-rownames(dge2[which(dge2$log2FoldChange>1),])
dn_de<-rownames(dge2[which(dge2$log2FoldChange<1),])
res <- list("dge" = zz, "up_de" = up_de, "dn_de" = dn_de)
res
}


#################################################
# define limma function
##################################################
limma<-function(y) {
res=NULL
label="simulate"
samplesheet<-as.data.frame(colnames(y))
colnames(samplesheet)="sample"
samplesheet$trt<-as.numeric(grepl("trt",colnames(y)))
design<-model.matrix(~samplesheet$trt)
rownames(design)=samplesheet$sample
y<-y[which(rowSums(y)/ncol(y)>=(10)),]
z<-DGEList(counts=y)
z <- calcNormFactors(z)
v <- voom(z,design,plot=F)
fit <- lmFit(v, design)
fit.de <- eBayes(fit, robust=TRUE)
dge<-topTable(fit.de,n=Inf)
dge<-dge[order(dge$adj.P.Val),]
dge2<-subset(dge,adj.P.Val<0.05)
up_de<-rownames(dge2[which(dge2$logFC>1),])
dn_de<-rownames(dge2[which(dge2$logFC<1),])
res <- list("dge" = dge, "up_de" = up_de, "dn_de" = dn_de)
res
}


##################################
# aggregate script
##################################
agg_dge<-function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC) {
#N_REPS=5 ; SUM_COUNT=10000000 ; VARIANCE=0.3 ; FRAC_DE=0.2 ; FC=1 ; SIMS=10 ; DGE_FUNC="edger"
xxx<-RepParallel(SIMS,simrna(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC), simplify=F, mc.cores = detectCores() )
dge<-mclapply( sapply(xxx,"[",1) , DGE_FUNC , mc.cores = detectCores() )
ups<-sapply(xxx,"[",2)
dns<-sapply(xxx,"[",3)
ups_edger<-sapply(dge,"[",2)
dns_edger<-sapply(dge,"[",3)
true_pos_up<-as.numeric(mapply( function(x,y) length(intersect(x,y)) , ups ,  ups_edger ))
true_pos_dn<-as.numeric(mapply( function(x,y) length(intersect(x,y)) , dns ,  dns_edger ))
false_pos_up<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , ups_edger ,  ups ))
false_pos_dn<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , dns_edger , dns ))
false_neg_up<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , ups ,  ups_edger ))
false_neg_dn<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , dns ,  dns_edger ))

true_pos<-mean(true_pos_up+true_pos_dn)
false_pos<-mean(false_pos_up+false_pos_dn)
false_neg<-mean(false_neg_up+false_neg_dn)
nrows<-as.numeric(lapply( sapply(xxx,"[",1 ), nrow))
true_neg<-mean(nrows-(true_pos+false_pos+false_neg))

p<-true_pos/(true_pos+false_pos)
r<-true_pos/(true_pos+false_neg)
f<-2*p*r/(p+r)

dge_res<-data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,true_pos,false_pos,true_neg,false_neg,p,r,f)
dge_res
}
#res1<-agg_dge(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS)




#sanity checks
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
N_REPS=5
SIMS=10
res=NULL
for (DGE_FUNC in c("edger","edger_ql","deseq","limma")) {
  for ( SUM_COUNT in c(10000000,40000000,100000000)) {
    for ( FRAC_DE in c(0.01,0.05,0.1,0.25)) {
      for (FC in c(1)) {
        for  ( VARIANCE in c(0,0.2,0.4,0.6)) {
          res_new<-agg_dge(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,10,DGE_FUNC)
          res=rbind(res,res_new)
        }
      }
    }
  }
}

rownames(res)<-paste(res$DGE_FUNC,res$SUM_COUNT,res$FRAC_DE,res$VARIANCE)
save.image(file="simde.RData")

res_edger<-res[which(res$DGE_FUNC=="edger"),]
res_edger<-res_edger[which(res_edger$FRAC_DE==0.1),]
res_edgerql<-res[which(res$DGE_FUNC=="edger_ql"),]
res_edgerql<-res_edgerql[which(res_edgerql$FRAC_DE==0.1),]
res_deseq<-res[which(res$DGE_FUNC=="deseq"),]
res_deseq<-res_deseq[which(res_deseq$FRAC_DE==0.1),]
res_limma<-res[which(res$DGE_FUNC=="limma"),]
res_limma<-res_limma[which(res_limma$FRAC_DE==0.1),]

pdf(file="barplots.pdf") ; par(mfrow=c(2,2)) ; par(mar=c(4,8,3,2))

barplot(prop.table(as.matrix(t(apply(res_edger[,8:11],2,rev))),2),las=1,horiz=T,xlab="Gene proportion", 
 cex.names=0.5,cex.main=0.9,main="edgeR 5 reps varied variance and DEG%",
 legend.text = TRUE,args.legend=list(x="topright",bty="n",inset=c(-0.05, 0)))

barplot(prop.table(as.matrix(t(apply(res_edgerql[,8:11],2,rev))),2),las=1,horiz=T,xlab="Gene proportion", 
 cex.names=0.5,cex.main=0.9,main="edgeR QL 5 reps varied variance and DEG%",
 legend.text = F,args.legend=list(x="topright",bty="n",inset=c(-0.05, 0)))

barplot(prop.table(as.matrix(t(apply(res_deseq[,8:11],2,rev))),2),las=1,horiz=T,xlab="Gene proportion", 
 cex.names=0.5,cex.main=0.9,main="DESeq 5 reps varied variance and DEG%",
 legend.text = F,args.legend=list(x="topright",bty="n",inset=c(-0.05, 0)))

barplot(prop.table(as.matrix(t(apply(res_limma[,8:11],2,rev))),2),las=1,horiz=T,xlab="Gene proportion", 
 cex.names=0.5,cex.main=0.9,main="limma 5 reps varied variance and DEG%",
 legend.text = F,args.legend=list(x="topright",bty="n",inset=c(-0.05, 0)))

dev.off()

res_edger_1e7<-res_edger[which(res_edger$SUM_COUNT==1E+7),]
res_edger_1e7_v0<-res_edger_1e7[which(res_edger_1e7$VARIANCE==0),]
res_edger_1e7_v2<-res_edger_1e7[which(res_edger_1e7$VARIANCE==0.2),]
res_edger_1e7_v4<-res_edger_1e7[which(res_edger_1e7$VARIANCE==0.4),]
res_edger_1e7_v6<-res_edger_1e7[which(res_edger_1e7$VARIANCE==0.6),]

res_edger_4e7<-res_edger[which(res_edger$SUM_COUNT==4E+7),]
res_edger_4e7_v0<-res_edger_4e7[which(res_edger_4e7$VARIANCE==0),]
res_edger_4e7_v2<-res_edger_4e7[which(res_edger_4e7$VARIANCE==0.2),]
res_edger_4e7_v4<-res_edger_4e7[which(res_edger_4e7$VARIANCE==0.4),]
res_edger_4e7_v6<-res_edger_4e7[which(res_edger_4e7$VARIANCE==0.6),]

res_edger_1e8<-res_edger[which(res_edger$SUM_COUNT==1E+8),]
res_edger_1e8_v0<-res_edger_1e8[which(res_edger_1e8$VARIANCE==0),]
res_edger_1e8_v2<-res_edger_1e8[which(res_edger_1e8$VARIANCE==0.2),]
res_edger_1e8_v4<-res_edger_1e8[which(res_edger_1e8$VARIANCE==0.4),]
res_edger_1e8_v6<-res_edger_1e8[which(res_edger_1e8$VARIANCE==0.6),]

res_edgerql_1e7<-res_edgerql[which(res_edgerql$SUM_COUNT==1E+7),]
res_edgerql_1e7_v0<-res_edgerql_1e7[which(res_edgerql_1e7$VARIANCE==0),]
res_edgerql_1e7_v2<-res_edgerql_1e7[which(res_edgerql_1e7$VARIANCE==0.2),]
res_edgerql_1e7_v4<-res_edgerql_1e7[which(res_edgerql_1e7$VARIANCE==0.4),]
res_edgerql_1e7_v6<-res_edgerql_1e7[which(res_edgerql_1e7$VARIANCE==0.6),]

res_edgerql_4e7<-res_edgerql[which(res_edgerql$SUM_COUNT==4E+7),]
res_edgerql_4e7_v0<-res_edgerql_4e7[which(res_edgerql_4e7$VARIANCE==0),]
res_edgerql_4e7_v2<-res_edgerql_4e7[which(res_edgerql_4e7$VARIANCE==0.2),]
res_edgerql_4e7_v4<-res_edgerql_4e7[which(res_edgerql_4e7$VARIANCE==0.4),]
res_edgerql_4e7_v6<-res_edgerql_4e7[which(res_edgerql_4e7$VARIANCE==0.6),]

res_edgerql_1e8<-res_edgerql[which(res_edgerql$SUM_COUNT==1E+8),]
res_edgerql_1e8_v0<-res_edgerql_1e8[which(res_edgerql_1e8$VARIANCE==0),]
res_edgerql_1e8_v2<-res_edgerql_1e8[which(res_edgerql_1e8$VARIANCE==0.2),]
res_edgerql_1e8_v4<-res_edgerql_1e8[which(res_edgerql_1e8$VARIANCE==0.4),]
res_edgerql_1e8_v6<-res_edgerql_1e8[which(res_edgerql_1e8$VARIANCE==0.6),]

res_deseq_1e7<-res_deseq[which(res_deseq$SUM_COUNT==1E+7),]
res_deseq_1e7_v0<-res_deseq_1e7[which(res_deseq_1e7$VARIANCE==0),]
res_deseq_1e7_v2<-res_deseq_1e7[which(res_deseq_1e7$VARIANCE==0.2),]
res_deseq_1e7_v4<-res_deseq_1e7[which(res_deseq_1e7$VARIANCE==0.4),]
res_deseq_1e7_v6<-res_deseq_1e7[which(res_deseq_1e7$VARIANCE==0.6),]

res_deseq_4e7<-res_deseq[which(res_deseq$SUM_COUNT==4E+7),]
res_deseq_4e7_v0<-res_deseq_4e7[which(res_deseq_4e7$VARIANCE==0),]
res_deseq_4e7_v2<-res_deseq_4e7[which(res_deseq_4e7$VARIANCE==0.2),]
res_deseq_4e7_v4<-res_deseq_4e7[which(res_deseq_4e7$VARIANCE==0.4),]
res_deseq_4e7_v6<-res_deseq_4e7[which(res_deseq_4e7$VARIANCE==0.6),]

res_deseq_1e8<-res_deseq[which(res_deseq$SUM_COUNT==1E+8),]
res_deseq_1e8_v0<-res_deseq_1e8[which(res_deseq_1e8$VARIANCE==0),]
res_deseq_1e8_v2<-res_deseq_1e8[which(res_deseq_1e8$VARIANCE==0.2),]
res_deseq_1e8_v4<-res_deseq_1e8[which(res_deseq_1e8$VARIANCE==0.4),]
res_deseq_1e8_v6<-res_deseq_1e8[which(res_deseq_1e8$VARIANCE==0.6),]

res_limma_1e7<-res_limma[which(res_limma$SUM_COUNT==1E+7),]
res_limma_1e7_v0<-res_limma_1e7[which(res_limma_1e7$VARIANCE==0),]
res_limma_1e7_v2<-res_limma_1e7[which(res_limma_1e7$VARIANCE==0.2),]
res_limma_1e7_v4<-res_limma_1e7[which(res_limma_1e7$VARIANCE==0.4),]
res_limma_1e7_v6<-res_limma_1e7[which(res_limma_1e7$VARIANCE==0.6),]

res_limma_4e7<-res_limma[which(res_limma$SUM_COUNT==4E+7),]
res_limma_4e7_v0<-res_limma_4e7[which(res_limma_4e7$VARIANCE==0),]
res_limma_4e7_v2<-res_limma_4e7[which(res_limma_4e7$VARIANCE==0.2),]
res_limma_4e7_v4<-res_limma_4e7[which(res_limma_4e7$VARIANCE==0.4),]
res_limma_4e7_v6<-res_limma_4e7[which(res_limma_4e7$VARIANCE==0.6),]

res_limma_1e8<-res_limma[which(res_limma$SUM_COUNT==1E+8),]
res_limma_1e8_v0<-res_limma_1e8[which(res_limma_1e8$VARIANCE==0),]
res_limma_1e8_v2<-res_limma_1e8[which(res_limma_1e8$VARIANCE==0.2),]
res_limma_1e8_v4<-res_limma_1e8[which(res_limma_1e8$VARIANCE==0.4),]
res_limma_1e8_v6<-res_limma_1e8[which(res_limma_1e8$VARIANCE==0.6),]

pdf(file="pr.pdf")
par(mfrow=c(2,2))

par(mar=c(4,4,3,2))
plot(res_edger_1e7_v0$r,res_edger_1e7_v0$p,xlab="recall",ylab="precision",pch=15,col="red",xlim=c(0,1),ylim=c(0,1))
points(res_edger_1e7_v2$r,res_edger_1e7_v0$p,xlab="recall",ylab="precision",pch=16,col="red",xlim=c(0,1),ylim=c(0,1))
points(res_edger_1e7_v4$r,res_edger_1e7_v4$p,xlab="recall",ylab="precision",pch=17,col="red",xlim=c(0,1),ylim=c(0,1))
points(res_edger_1e7_v6$r,res_edger_1e7_v6$p,xlab="recall",ylab="precision",pch=18,col="red",xlim=c(0,1),ylim=c(0,1))

points(res_edger_4e7_v0$r,res_edger_4e7_v0$p,xlab="recall",ylab="precision",pch=15,col="blue",xlim=c(0,1),ylim=c(0,1))
points(res_edger_4e7_v2$r,res_edger_4e7_v0$p,xlab="recall",ylab="precision",pch=16,col="blue",xlim=c(0,1),ylim=c(0,1))
points(res_edger_4e7_v4$r,res_edger_4e7_v4$p,xlab="recall",ylab="precision",pch=17,col="blue",xlim=c(0,1),ylim=c(0,1))
points(res_edger_4e7_v6$r,res_edger_4e7_v6$p,xlab="recall",ylab="precision",pch=18,col="blue",xlim=c(0,1),ylim=c(0,1))

points(res_edger_1e8_v0$r,res_edger_1e8_v0$p,xlab="recall",ylab="precision",pch=15,col="black",xlim=c(0,1),ylim=c(0,1))
points(res_edger_1e8_v2$r,res_edger_1e8_v2$p,xlab="recall",ylab="precision",pch=16,col="black",xlim=c(0,1),ylim=c(0,1))
points(res_edger_1e8_v4$r,res_edger_1e8_v4$p,xlab="recall",ylab="precision",pch=17,col="black",xlim=c(0,1),ylim=c(0,1))
points(res_edger_1e8_v6$r,res_edger_1e8_v6$p,xlab="recall",ylab="precision",pch=18,col="black",xlim=c(0,1),ylim=c(0,1))

legend(0.80,1,legend=c("10M","40M","100M"),col=c("red", "blue","black") ,pch=19,cex=0.6,title="read depth")
legend(0.80,0.7,legend=c("0","0.2","0.4","0.6"),col=c("black") ,pch=15:19,cex=0.6,title="added variance")
mtext("edgeR 5 reps, 10% DEG" ,cex=0.8); grid()

par(mar=c(4,4,3,2))
plot(res_edgerql_1e7_v0$r,res_edgerql_1e7_v0$p,xlab="recall",ylab="precision",pch=15,col="red",xlim=c(0,1),ylim=c(0,1))
points(res_edgerql_1e7_v2$r,res_edgerql_1e7_v0$p,xlab="recall",ylab="precision",pch=16,col="red",xlim=c(0,1),ylim=c(0,1))
points(res_edgerql_1e7_v4$r,res_edgerql_1e7_v4$p,xlab="recall",ylab="precision",pch=17,col="red",xlim=c(0,1),ylim=c(0,1))
points(res_edgerql_1e7_v6$r,res_edgerql_1e7_v6$p,xlab="recall",ylab="precision",pch=18,col="red",xlim=c(0,1),ylim=c(0,1))

points(res_edgerql_4e7_v0$r,res_edgerql_4e7_v0$p,xlab="recall",ylab="precision",pch=15,col="blue",xlim=c(0,1),ylim=c(0,1))
points(res_edgerql_4e7_v2$r,res_edgerql_4e7_v0$p,xlab="recall",ylab="precision",pch=16,col="blue",xlim=c(0,1),ylim=c(0,1))
points(res_edgerql_4e7_v4$r,res_edgerql_4e7_v4$p,xlab="recall",ylab="precision",pch=17,col="blue",xlim=c(0,1),ylim=c(0,1))
points(res_edgerql_4e7_v6$r,res_edgerql_4e7_v6$p,xlab="recall",ylab="precision",pch=18,col="blue",xlim=c(0,1),ylim=c(0,1))

points(res_edgerql_1e8_v0$r,res_edgerql_1e8_v0$p,xlab="recall",ylab="precision",pch=15,col="black",xlim=c(0,1),ylim=c(0,1))
points(res_edgerql_1e8_v2$r,res_edgerql_1e8_v2$p,xlab="recall",ylab="precision",pch=16,col="black",xlim=c(0,1),ylim=c(0,1))
points(res_edgerql_1e8_v4$r,res_edgerql_1e8_v4$p,xlab="recall",ylab="precision",pch=17,col="black",xlim=c(0,1),ylim=c(0,1))
points(res_edgerql_1e8_v6$r,res_edgerql_1e8_v6$p,xlab="recall",ylab="precision",pch=18,col="black",xlim=c(0,1),ylim=c(0,1))

legend(0.80,1,legend=c("10M","40M","100M"),col=c("red", "blue","black") ,pch=19,cex=0.6,title="read depth")
legend(0.80,0.7,legend=c("0","0.2","0.4","0.6"),col=c("black") ,pch=15:19,cex=0.6,title="added variance")
mtext("edgeR QL 5 reps, 10% DEG" ,cex=0.8); grid()

par(mar=c(4,4,3,2))
plot(res_deseq_1e7_v0$r,res_deseq_1e7_v0$p,xlab="recall",ylab="precision",pch=15,col="red",xlim=c(0,1),ylim=c(0,1))
points(res_deseq_1e7_v2$r,res_deseq_1e7_v0$p,xlab="recall",ylab="precision",pch=16,col="red",xlim=c(0,1),ylim=c(0,1))
points(res_deseq_1e7_v4$r,res_deseq_1e7_v4$p,xlab="recall",ylab="precision",pch=17,col="red",xlim=c(0,1),ylim=c(0,1))
points(res_deseq_1e7_v6$r,res_deseq_1e7_v6$p,xlab="recall",ylab="precision",pch=18,col="red",xlim=c(0,1),ylim=c(0,1))

points(res_deseq_4e7_v0$r,res_deseq_4e7_v0$p,xlab="recall",ylab="precision",pch=15,col="blue",xlim=c(0,1),ylim=c(0,1))
points(res_deseq_4e7_v2$r,res_deseq_4e7_v0$p,xlab="recall",ylab="precision",pch=16,col="blue",xlim=c(0,1),ylim=c(0,1))
points(res_deseq_4e7_v4$r,res_deseq_4e7_v4$p,xlab="recall",ylab="precision",pch=17,col="blue",xlim=c(0,1),ylim=c(0,1))
points(res_deseq_4e7_v6$r,res_deseq_4e7_v6$p,xlab="recall",ylab="precision",pch=18,col="blue",xlim=c(0,1),ylim=c(0,1))

points(res_deseq_1e8_v0$r,res_deseq_1e8_v0$p,xlab="recall",ylab="precision",pch=15,col="black",xlim=c(0,1),ylim=c(0,1))
points(res_deseq_1e8_v2$r,res_deseq_1e8_v2$p,xlab="recall",ylab="precision",pch=16,col="black",xlim=c(0,1),ylim=c(0,1))
points(res_deseq_1e8_v4$r,res_deseq_1e8_v4$p,xlab="recall",ylab="precision",pch=17,col="black",xlim=c(0,1),ylim=c(0,1))
points(res_deseq_1e8_v6$r,res_deseq_1e8_v6$p,xlab="recall",ylab="precision",pch=18,col="black",xlim=c(0,1),ylim=c(0,1))

legend(0.80,1,legend=c("10M","40M","100M"),col=c("red", "blue","black") ,pch=19,cex=0.6,title="read depth")
legend(0.80,0.7,legend=c("0","0.2","0.4","0.6"),col=c("black") ,pch=15:19,cex=0.6,title="added variance")
mtext("DESeq 5 reps, 10% DEG" ,cex=0.8) ; grid()

par(mar=c(4,4,3,2))
plot(res_limma_1e7_v0$r,res_limma_1e7_v0$p,xlab="recall",ylab="precision",pch=15,col="red",xlim=c(0,1),ylim=c(0,1))
points(res_limma_1e7_v2$r,res_limma_1e7_v0$p,xlab="recall",ylab="precision",pch=16,col="red",xlim=c(0,1),ylim=c(0,1))
points(res_limma_1e7_v4$r,res_limma_1e7_v4$p,xlab="recall",ylab="precision",pch=17,col="red",xlim=c(0,1),ylim=c(0,1))
points(res_limma_1e7_v6$r,res_limma_1e7_v6$p,xlab="recall",ylab="precision",pch=18,col="red",xlim=c(0,1),ylim=c(0,1))

points(res_limma_4e7_v0$r,res_limma_4e7_v0$p,xlab="recall",ylab="precision",pch=15,col="blue",xlim=c(0,1),ylim=c(0,1))
points(res_limma_4e7_v2$r,res_limma_4e7_v0$p,xlab="recall",ylab="precision",pch=16,col="blue",xlim=c(0,1),ylim=c(0,1))
points(res_limma_4e7_v4$r,res_limma_4e7_v4$p,xlab="recall",ylab="precision",pch=17,col="blue",xlim=c(0,1),ylim=c(0,1))
points(res_limma_4e7_v6$r,res_limma_4e7_v6$p,xlab="recall",ylab="precision",pch=18,col="blue",xlim=c(0,1),ylim=c(0,1))

points(res_limma_1e8_v0$r,res_limma_1e8_v0$p,xlab="recall",ylab="precision",pch=15,col="black",xlim=c(0,1),ylim=c(0,1))
points(res_limma_1e8_v2$r,res_limma_1e8_v2$p,xlab="recall",ylab="precision",pch=16,col="black",xlim=c(0,1),ylim=c(0,1))
points(res_limma_1e8_v4$r,res_limma_1e8_v4$p,xlab="recall",ylab="precision",pch=17,col="black",xlim=c(0,1),ylim=c(0,1))
points(res_limma_1e8_v6$r,res_limma_1e8_v6$p,xlab="recall",ylab="precision",pch=18,col="black",xlim=c(0,1),ylim=c(0,1))

legend(0.80,1,legend=c("10M","40M","100M"),col=c("red", "blue","black") ,pch=19,cex=0.6,title="read depth")
legend(0.80,0.7,legend=c("0","0.2","0.4","0.6"),col=c("black") ,pch=15:19,cex=0.6,title="added variance")
mtext("Limma 5 reps, 10% DEG" ,cex=0.8); grid()

dev.off()

q()








res
smm_10m <- list("smm_0de"=smm_0de,"smm_1de"=smm_1de,"smm_2de"=smm_2de,"smm_3de"=smm_3de)
save.image(file="simde.RData")

pdf(file="10Mreads.pdf")
par(mfrow=c(2,4))

barplot(prop.table(as.matrix(t(apply(smm_0de[,1:4],2,rev))),2),las=1,
 ylab="Added noise",horiz=T,xlab="Gene proportion",cex.main=0.9,main="edger 10M reads no DGEs",
 legend.text = TRUE,args.legend=list(x="topright",bty="n",inset=c(-0.05, 0)))

barplot(prop.table(as.matrix(t(apply(smm_1de[,1:4],2,rev))),2),las=1,
 ylab="Added noise",horiz=T,xlab="Gene proportion",cex.main=0.9,main="edger 10M reads 5% DGEs",
 legend.text = FALSE,args.legend=list(x="topright",bty="n",inset=c(-0.05, 0)))

dotchart(rev(smm_0de$f),pch=19,labels=rev(rownames(smm)),xlim=c(0,1),cex.main=0.9,main="edger 10M reads no DGEs")
legend("topright", pch = c(19,19,19), col = c("black", "red", "blue"),legend = c("F1","p","r"))
grid()
points(rev(smm_0de$p),1:5,pch=19,col="red")
points(rev(smm_0de$r),1:5,pch=19,col="blue")
points(rev(smm_0de$f),1:5,pch=19,col="black")

dotchart(rev(smm_1de$f),pch=19,labels=rev(rownames(smm)),xlim=c(0,1),cex.main=0.9,main="edger 10M reads 5% DGEs")
#legend("topright", pch = c(19,19,19), col = c("black", "red", "blue"),legend = c("F1","p","r"))
grid()
points(rev(smm_1de$p),1:5,pch=19,col="red")
points(rev(smm_1de$r),1:5,pch=19,col="blue")
points(rev(smm_1de$f),1:5,pch=19,col="black")

