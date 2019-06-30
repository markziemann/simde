library("tidyverse")
library("parallel")
library("topconfects")
library("edgeR")
library("DESeq2")
library("limma")
library("ABSSeq")

########################################
# simulate some gene exression data
########################################
#inputs=(a,N_REPS,SUM_COUNT,
simrna<-function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC) {
library("edgeR")

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
# RepParallel usage
#xxx<-RepParallel(10,simrna(a,5,10000000,0.2,20), simplify=F, mc.cores = detectCores() )

#################################################
# define edgeR classic function
##################################################
edger<-function(y) {
library("limma")
library("edgeR")
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

up_de<-dge2[which(dge2$logFC>0),1]
dn_de<-dge2[which(dge2$logFC<0),1]

res <- list("dge" = dge, "up_de" = up_de, "dn_de" = dn_de)
res
}

#################################################
# define edgeR QL function
##################################################
edger_ql<-function(y) {
library("limma")
library("edgeR")
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

up_de<-dge2[which(dge2$logFC>0),1]
dn_de<-dge2[which(dge2$logFC<0),1]

res <- list("dge" = dge, "up_de" = up_de, "dn_de" = dn_de)
res
}

#################################################
# define DESeq2 function
##################################################
deseq<-function(y) {
library("DESeq2")
res=NULL
label="simulate"
samplesheet<-as.data.frame(colnames(y))
colnames(samplesheet)="sample"
samplesheet$trt<-factor(as.numeric(grepl("trt",colnames(y))))
dds <- DESeqDataSetFromMatrix(countData = y, colData = samplesheet, design = ~ trt )
res <- DESeq(dds)
z<- DESeq2::results(res)
vsd <- vst(dds, blind=FALSE)
zz<-cbind(z,assay(vsd))
zz<-as.data.frame(zz[order(zz$padj),])
dge2<-subset(zz,padj<0.05)
up_de<-rownames(dge2[which(dge2$log2FoldChange>0),])
dn_de<-rownames(dge2[which(dge2$log2FoldChange<0),])
res <- list("dge" = zz, "up_de" = up_de, "dn_de" = dn_de)
res
}


#################################################
# define limma function
##################################################
limma<-function(y) {
library("limma")
library("edgeR")
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
up_de<-rownames(dge2[which(dge2$logFC>0),])
dn_de<-rownames(dge2[which(dge2$logFC<0),])
res <- list("dge" = dge, "up_de" = up_de, "dn_de" = dn_de)
res
}


#################################################
# define absseq function
##################################################
absseq<-function(y) {
library("ABSSeq")
res=NULL
label="simulate"
samplesheet<-as.data.frame(colnames(y))
colnames(samplesheet)="sample"
samplesheet$trt<-as.numeric(grepl("trt",colnames(y)))
obj<-ABSDataSet(y, factor(samplesheet$trt))  #default normalisation is qtotal
obj<-ABSSeq(obj)
dge<- as.data.frame(cbind(obj$Amean,obj$Bmean,obj$foldChange,obj$pvalue,obj$adj.pvalue))
colnames(dge)=c("Amean","Bmean","logFC","PValue","FDR")
dge<-dge[order(dge$PValue),]
dge2<-subset(dge,FDR<0.05)
up_de<-rownames(dge2[which(dge2$logFC>0),])
dn_de<-rownames(dge2[which(dge2$logFC<0),])
res <- list("dge" = dge, "up_de" = up_de, "dn_de" = dn_de)
res
}

## here is how to run topconfects from an edger analysis
## slow but works
#econfects <- edger_confects(fit, coef=2, fdr=0.05)
#econfects <- edger_confects(fit, coef=2, fdr=0.5)
# confects_plot(econfects)

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


