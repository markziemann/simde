library(edgeR)
library(tidyverse)
library(parallel)

#a is orig expression data
a<-read.table("start_data/ERR2539161/ERR2539161.se.tsv")

########################################
# simulate some gene exression data
########################################
#inputs=(a,N_REPS,SUM_COUNT,
simrna<-function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC) {
#N_REPS=5 ; SUM_COUNT=10000000 ; VARIANCE=0.5 ; FRAC_DE=0 ; FC=0
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
  message("incorporate changes")
  df2<-merge(df,ALL_DE,by.x=0,by.y="Gene",all.x=T)
  df2[is.na(df2)] <- 1
} else {
  df2<-as.data.frame( df )
  df2$FC<- 1
  UP_DE=NULL
  DN_DE=NULL
}
ODD_COLS=(2:(ncol(df2)-1))[c(TRUE,FALSE)]
EVEN_COLS=(2:(ncol(df2)-1))[c(FALSE,TRUE)]
controls<-df2[,ODD_COLS]
colnames(controls)=paste0( "ctrl_" ,1:ncol(controls) )
treatments<-round(df2[,EVEN_COLS]*df2$FC)
colnames(treatments)=paste0( "trt_" ,1:ncol(treatments) )
x<-cbind(controls,treatments)
rownames(x)=rownames(df)
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
# define edgeR function
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
write.table(dge,file=paste(label,"_edger.tsv",sep=""),sep="\t",quote=F,row.names=F)
##output rank file
rnk<-as.data.frame(-log10(dge$PValue+1E-300)/sign(dge$logFC))
colnames(rnk)="Score"
rownames(rnk)=dge$Row.names
write.table(rnk,file=paste(label,".rnk",sep=""),sep="\t",quote=F,row.names=F)
#plots
plotname=paste(label,"_plots.pdf",sep="")
pdf(plotname)
#MA plot
dge2<-subset(dge,FDR<0.05)
up_de<-dge2[which(dge2$logFC>1),1]
dn_de<-dge2[which(dge2$logFC<1),1]
cntdge=nrow(dge2)
cntup=nrow(subset(dge2,logFC>0))
cntdn=nrow(subset(dge2,logFC<0))
HEADING=paste(label, cntdge, "DGEs", cntup, "up", cntdn, "dn")
plot(dge$logCPM,dge$logFC,main=HEADING,col="gray",pch=19,cex=0.5,xlab="log2(CPM)",ylab="log2(Fold Change)")
points(dge2$logCPM,dge2$logFC,col="red",pch=19,cex=0.5)
text(dge2$logCPM+1,dge2$logFC,labels=dge2$Row.names,cex=0.6)
#volcano plot
dge2<-subset(dge,FDR<0.05)
cntdge=nrow(dge2)
cntup=nrow(subset(dge2,logFC>0))
cntdn=nrow(subset(dge2,logFC<0))
HEADING=paste(label, cntdge, "DGEs", cntup, "up", cntdn, "dn")
XMAX=max(dge$logFC+0.5)
XMIN=min(dge$logFC)
plot(dge$logFC,-log2(dge$PValue),xlim=c(XMIN,XMAX),main=HEADING,col="gray",pch=19,cex=0.5,xlab="log2(Fold Change)",ylab="log2(nom P value)")
points(dge2$logFC,-log2(dge2$PValue),col="red",pch=19,cex=0.5)
text(dge2$logFC+0.25,-log2(dge2$PValue),labels=dge2$Row.names,cex=0.6)
#Heatnap
cols=paste(grep(".x", names(dge), value = TRUE))
zz<-as.data.frame(dge[,cols])
rownames(zz)=dge$Row.names
zz<-zz[1:50,]
heatmap(as.matrix(zz),margins = c(10,16))
#Gene Barchart
rownames(dge)=dge$Row.names
par(mfrow=c(3,3))
for (i in 1:9) {
 COLS=grep(".x",colnames(dge))
 b<-t(dge[i,COLS])
 logFC=dge[i,grep("logFC",colnames(dge)) ]
 Pval=dge[i,grep("PVal",colnames(dge)) ]
 FDR=dge[i,grep("FDR",colnames(dge)) ]
 HEADER=colnames(b)
 SUBHEAD=paste("logFC=",signif(logFC,3)," Pval=",signif(Pval,3)," FDR=",signif(FDR,3),sep="")
 par(mar=c(10,8,4,5))
 barplot(b,beside=T,main=HEADER,names.arg=sub(".x","",rownames(b)),las=2,ylab="CPM",cex.lab=0.8, cex.axis=0.8, cex.main=1, cex.sub=0.9,cex.names=0.6)
 mtext(SUBHEAD,cex=0.5)
}
dev.off()
res <- list("dge" = dge, "up_de" = up_de, "dn_de" = dn_de, "rnk" = rnk)
res
}

##################################
# aggregate script
##################################
agg_edger<-function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS) {
#N_REPS=5 ; SUM_COUNT=10000000 ; VARIANCE=0.3 ; FRAC_DE=0.2 ; FC=1 ; SIMS=10
xxx<-RepParallel(SIMS,simrna(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC), simplify=F, mc.cores = detectCores() )
dge<-mclapply( sapply(xxx,"[",1) , edger, mc.cores = detectCores() )
ups<-sapply(xxx,"[",2)
dns<-sapply(xxx,"[",3)
ups_edger<-sapply(dge,"[",2)
dns_edger<-sapply(dge,"[",3)
true_pos_up<-as.numeric(mapply( function(x,y) length(intersect(x,y)) , ups ,  ups_edger ))
true_pos_dn<-as.numeric(mapply( function(x,y) length(intersect(x,y)) , dns ,  dns_edger ))
true_pos<-true_pos_up+true_pos_dn
false_pos_up<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , ups_edger ,  ups ))
false_pos_dn<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , dns_edger , dns ))
false_pos<-false_pos_up+false_pos_dn
false_neg_up<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , ups ,  ups_edger ))
false_neg_dn<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , dns ,  dns_edger ))
false_neg<-false_neg_up+false_neg_dn
nrows<-as.numeric(lapply( sapply(xxx,"[",1 ), nrow))
true_neg<-nrows-(true_pos+false_pos+false_neg)
edger_res<-cbind(true_pos,false_pos,true_neg,false_neg)
edger_res
}
#res1<-agg_edger(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS)


#sanity checks
pdf(file="sanity_check.pdf")
x1<-simrna(a,5,10000000,0,0,0)
colnames(x1$x)<-paste("x1",colnames(x1$x),sep="_")
x2<-simrna(a,5,10000000,1,0,0)
colnames(x2$x)<-paste("x2",colnames(x2$x),sep="_")
xx<-merge(x1$x,x2$x,by=0)
heatmap(cor(xx[,2:ncol(xx)]),scale="none",main="check1")

x1<-simrna(a,5,10000000,0,0,0)
colnames(x1$x)<-paste("x1",colnames(x1$x),sep="_")
x2<-simrna(a,5,10000000,0,0.1,1)
colnames(x2$x)<-paste("x2",colnames(x2$x),sep="_")
xx<-merge(x1$x,x2$x,by=0)
heatmap(cor(xx[,2:ncol(xx)]),scale="none",main="check2")

x1<-simrna(a,5,10000000,0,0,0)
colnames(x1$x)<-paste("x1",colnames(x1$x),sep="_")
x2<-simrna(a,5,10000000,1,0.1,1)
colnames(x2$x)<-paste("x2",colnames(x2$x),sep="_")
xx<-merge(x1$x,x2$x,by=0)
heatmap(cor(xx[,2:ncol(xx)]),scale="none",main="check3")

dev.off()


# No DE just adding noise
res0<-agg_edger(a,5,10000000,0,0,0,10)
res1<-agg_edger(a,5,10000000,0.2,0,0,10)
res2<-agg_edger(a,5,10000000,0.4,0,0,10)
res3<-agg_edger(a,5,10000000,0.6,0,0,10)
res4<-agg_edger(a,5,10000000,0.8,0,0,10)
res5<-agg_edger(a,5,10000000,1,0,0,10)
res6<-agg_edger(a,5,10000000,1.5,0,0,10)
res7<-agg_edger(a,5,10000000,2,0,0,10)
res_0de <- list("v0"=res0,"v0.2"=res1,"v0.4"=res2,"v0.6"=res3,"v0.8"=res4,"v1"=res5,"v1.5"=res6,"v2"=res7)
lapply( res_0de , colMeans)

sum_0de<-t(matrix(unlist(lapply( res_0de , colMeans)),nrow=4))
colnames(sum_0de)=c("true_pos","false_pos","true_neg","false_neg")
rownames(sum_0de)=c("v0","v0.2","v0.4","v0.6","v0.8","v1","v1.5","v2")



res0<-agg_edger(a,5,10000000,0,0.05,1,10)
res1<-agg_edger(a,5,10000000,0.2,0.05,1,10)
res2<-agg_edger(a,5,10000000,0.4,0.05,1,10)
res3<-agg_edger(a,5,10000000,0.6,0.05,1,10)
res4<-agg_edger(a,5,10000000,0.8,0.05,1,10)
res5<-agg_edger(a,5,10000000,1,0.05,1,10)
res6<-agg_edger(a,5,10000000,1.5,0.05,1,10)
res7<-agg_edger(a,5,10000000,2,0.05,1,10)
res_1de <- list("v0"=res0,"v0.2"=res1,"v0.4"=res2,"v0.6"=res3,"v0.8"=res4,"v1"=res5,"v1.5"=res6,"v2"=res7)
lapply( res_1de , colMeans)

sum_1de<-t(matrix(unlist(lapply( res_1de , colMeans)),nrow=4))
colnames(sum_1de)=c("true_pos","false_pos","true_neg","false_neg")
rownames(sum_1de)=c("v0","v0.2","v0.4","v0.6","v0.8","v1","v1.5","v2")

res0<-agg_edger(a,5,10000000,0,0.1,1,10)
res1<-agg_edger(a,5,10000000,0.2,0.1,1,10)
res2<-agg_edger(a,5,10000000,0.4,0.1,1,10)
res3<-agg_edger(a,5,10000000,0.6,0.1,1,10)
res4<-agg_edger(a,5,10000000,0.8,0.1,1,10)
res5<-agg_edger(a,5,10000000,1,0.1,1,10)
res6<-agg_edger(a,5,10000000,1.5,0.1,1,10)
res7<-agg_edger(a,5,10000000,2,0.1,1,10)
res_2de <- list("v0"=res0,"v0.2"=res1,"v0.4"=res2,"v0.6"=res3,"v0.8"=res4,"v1"=res5,"v1.5"=res6,"v2"=res7)
lapply( res_1de , colMeans)

sum_2de<-t(matrix(unlist(lapply( res_2de , colMeans)),nrow=4))
colnames(sum_2de)=c("true_pos","false_pos","true_neg","false_neg")
rownames(sum_2de)=c("v0","v0.2","v0.4","v0.6","v0.8","v1","v1.5","v2")


res0<-agg_edger(a,5,10000000,0,0.25,1,10)
res1<-agg_edger(a,5,10000000,0.2,0.25,1,10)
res2<-agg_edger(a,5,10000000,0.4,0.25,1,10)
res3<-agg_edger(a,5,10000000,0.6,0.25,1,10)
res4<-agg_edger(a,5,10000000,0.8,0.25,1,10)
res5<-agg_edger(a,5,10000000,1,0.25,1,10)
res6<-agg_edger(a,5,10000000,1.5,0.25,1,10)
res7<-agg_edger(a,5,10000000,2,0.25,1,10)
res_3de <- list("v0"=res0,"v0.2"=res1,"v0.4"=res2,"v0.6"=res3,"v0.8"=res4,"v1"=res5,"v1.5"=res6,"v2"=res7)
lapply( res_2de , colMeans)

sum_3de<-t(matrix(unlist(lapply( res_3de , colMeans)),nrow=4))
colnames(sum_3de)=c("true_pos","false_pos","true_neg","false_neg")
colnames(sum_3de)=c("v0","v0.2","v0.4","v0.6","v0.8","v1","v1.5","v2")


#

#simrna(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC)



#pie(colSums(edger_res), paste(colnames(edger_res),colSums(edger_res)) ,main="edgeR result")
save.image(file="de_eval.RData")
