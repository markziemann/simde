library(edgeR)
library(tidyverse)
library(fgsea)
library(parallel)

#############################
# Todo - dont incorporate gene names at all.
#############################
#a is orig 
a<-read.table("start_data/ERR2539161/ERR2539161.se.tsv")

########################################
# simulate some gene exression data
########################################
#inputs=(a,N_REPS,SUM_COUNT,
simrna<-function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC) {
#N_REPS=5 #SUM_COUNT=10000000 #VARIANCE=0.2 #FRAC_DE=0.05 #FC=1
df = NULL
for (k in paste0("data",1:(N_REPS*2)))  {
        b<-thinCounts(a,target.size=SUM_COUNT)
        colnames(b)=k
        df = cbind(df,b)
     }
#create some random values centred around 1 with some% error
rand<-matrix(log2(rnorm(nrow(a)*N_REPS*2 , 2, VARIANCE)),ncol=N_REPS*2)
#incorporate the noise
df2<-round(df*rand)
#set any negative counts to zero
df2<-apply(df2, 2, function(x) {ifelse(x < 0, 0, x)}) 
message("prep fold changes")
#Number of differential genes
NDIF=round(nrow(df)*FRAC_DE)
#Make even
if ( NDIF%%2==1 ) { print("odd") ; NDIF=NDIF-1 }
#sample DE genes
DE_LIST<-sample(rownames(df2) , NDIF)
#split list in half
UP_DE<-sample(DE_LIST,length(DE_LIST)/2)
DN_DE<-setdiff(DE_LIST,UP_DE)
#reformat as df and add fold change
UP_DE<-as.data.frame(UP_DE)
colnames(UP_DE)="Gene"
UP_DE$V1<-2^FC
DN_DE<-as.data.frame(DN_DE)
colnames(DN_DE)="Gene"
DN_DE$V1<-2^-FC
ALL_DE<-rbind(DN_DE,UP_DE)
#Go back to list for downstream work
UP_DE<-UP_DE$Gene
DN_DE<-DN_DE$Gene
message("incorporate changes")
df2<-merge(df,ALL_DE,by.x=0,by.y="Gene",all.x=T)
df2[is.na(df2)] <- 1
ODD_COLS=(2:(ncol(df2)-1))[c(TRUE,FALSE)]
EVEN_COLS=(2:(ncol(df2)-1))[c(FALSE,TRUE)]
controls<-df2[,ODD_COLS]
colnames(controls)=paste0( "ctrl_" ,1:ncol(controls) )
treatments<-round(df2[,EVEN_COLS]*df2$V1)
colnames(treatments)=paste0( "trt_" ,1:ncol(treatments) )
x<-cbind(controls,treatments)
rownames(x)=df2$Row.names
#filter out genes that are not expressed
x<- x[which(rowSums(x)/ncol(x)>10),]
UP_DE<-intersect(UP_DE,rownames(x))
DN_DE<-intersect(DN_DE,rownames(x))
xx <- list("x" = x, "UP_DE" = UP_DE, "DN_DE" = DN_DE )
xx
}
#simrna(a,N_REPS,SUM_COUNT,VARIANCE,N_DIFF_PW,FC)
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
message("run edgeR")
##################################
#res<-edger(y)
#simrna(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC)
xxx<-RepParallel(10,simrna(a,5,10000000,0.3,0.2,1), simplify=F, mc.cores = detectCores() )
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

pie(colSums(edger_res), paste(colnames(edger_res),colSums(edger_res)) ,main="edgeR result")
#get numrows
#lapply( sapply(xxx,"[",1 ), nrow)

save.image(file="de_eval.RData")
