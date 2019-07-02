library("tidyverse")
library("parallel")
library("topconfects")
library("edgeR")
library("DESeq2")
library("limma")
library("ABSSeq")
library("stringi")
library("mitch")


########################################
# simulate some gene expression data
########################################
simrna<-function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,gsets) {

# N_REPS=5 ; SUM_COUNT=10000000 ; VARIANCE=0.2 ; FRAC_DE=0.05 ; FC=1 ; GMT="ReactomePathways.gmt"

library("edgeR")

df = NULL
for (k in paste0("data",1:(N_REPS*2)))  {
        b<-thinCounts(a,target.size=SUM_COUNT)
        colnames(b)=k
        df = cbind(df,b)
     }

# now need to only include gsets with 10 members in the 
gsets_sub<-which(unlist( lapply(gsets,function(x) { length(which(rownames(a) %in% as.character(unlist(x)))) >10 }  ) ) )
gsets<-gsets[which(names(gsets) %in% names(gsets_sub))]

#Number of differential genes
NDIF=round(length(gsets)*FRAC_DE)

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

  # sample some pathways to fiddle with
  DE_LIST<-sample(gsets , NDIF)

  # divide the list in 2 with half up and half down
  UP_LIST=sample(DE_LIST , NDIF/2)
  DN_LIST<-DE_LIST[!(DE_LIST %in% UP_LIST)]

  # now find a list of genes inside the pathways
  UP_DE<-unique(unlist(unname(UP_LIST)))
  # select the ones that are also in the profile
  UP_DE<-UP_DE[which(UP_DE %in% row.names(df))]

  # same for down genes
  DN_DE<-unique(unlist(unname(DN_LIST)))
  DN_DE<-DN_DE[which(DN_DE %in% row.names(df))]


  ITX<-intersect(UP_DE,DN_DE)
  # need to eliminate the overlapping ones for simplicity
  UP_DE<-setdiff(UP_DE,ITX)
  DN_DE<-setdiff(DN_DE,ITX)

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
  UP_LIST=NULL
  DN_LIST=NULL
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
xx <- list("x" = x,"UP_DE"=UP_DE,"DN_DE"=DN_DE,"UP_LIST"=UP_LIST,"DN_LIST"=DN_LIST)
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
de<-as.data.frame(topTags(lrt,n=Inf))
de$dispersion<-lrt$dispersion
de<-de[order(de$PValue),]
x[[6]]<-de[order(de$PValue),]
sig<-subset(de,FDR<0.05)
x[[7]]<-sig[which(sig$logFC>0),1]
x[[8]]<-sig[which(sig$logFC<0),1]
x 
}

#################################################
# define edgeR QL function
##################################################
edger_ql<-function(y) {
library("limma")
library("edgeR")
label="simulate"
y<-x[[1]]
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
de<-as.data.frame(topTags(lrt,n=Inf))
de$dispersion<-lrt$dispersion
de<-de[order(de$PValue),]
x[[6]]<-de[order(de$PValue),]
sig<-subset(de,FDR<0.05)
x[[7]]<-sig[which(sig$logFC>0),1]
x[[8]]<-sig[which(sig$logFC<0),1]
x
}

#################################################
# define DESeq2 function
##################################################
deseq2<-function(x) {
library("DESeq2")
label="simulate"
y<-x[[1]]
samplesheet<-as.data.frame(colnames(y))
colnames(samplesheet)="sample"
samplesheet$trt<-factor(as.numeric(grepl("trt",colnames(y))))
dds <- DESeqDataSetFromMatrix(countData = y, colData = samplesheet, design = ~ trt )
res <- DESeq(dds)
z<- DESeq2::results(res)
vsd <- vst(dds, blind=FALSE)
zz<-cbind(z,assay(vsd))
x[[6]]<-as.data.frame(zz[order(zz$padj),])
sig<-subset(zz,padj<0.05)
x[[7]]<-rownames(sig[which(sig$log2FoldChange>0),])
x[[8]]<-rownames(sig[which(sig$log2FoldChange<0),])
x
}


#################################################
# define limma function
##################################################
limma<-function(y) {
library("limma")
library("edgeR")
label="simulate"
y<-x[[1]]
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
x[[6]]<-dge[order(dge$adj.P.Val),]
sig<-subset(dge,adj.P.Val<0.05)
x[[7]]<-rownames(sig[which(sig$logFC>0),])
x[[8]]<-rownames(sig[which(sig$logFC<0),])
x
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
x[[6]]<-dge[order(dge$PValue),]
sig<-subset(dge,FDR<0.05)
x[[7]]<-rownames(dge2[which(dge2$logFC>0),])
x[[8]]<-rownames(dge2[which(dge2$logFC<0),])
x
}

## here is how to run topconfects from an edger analysis
## slow but works
#econfects <- edger_confects(fit, coef=2, fdr=0.05)
#econfects <- edger_confects(fit, coef=2, fdr=0.5)
# confects_plot(econfects)

#################################################
# define mitch function
##################################################
run_mitch<-function(y,DGE_FUNC,gsets, N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS) {
library("mitch")
dge<-sapply(y,"[",6)
names(dge)<-paste0("x",1:length(dge),sep="")

w<-mitch_import(dge, DGE_FUNC )
res<-mitch_calc(w,gsets,priority="significance")
padj<-apply( res$manova_result[,grep("^p.x",colnames(res$manova_result))] , 2 , p.adjust, method="BH")
colnames(padj)<-paste0("adj.",colnames(padj),sep="")
res$manova_result<-data.frame(res$manova_result,padj)


# now to extract the observed changes
for (N in 1:length(grep("^p.x",colnames(res$manova_result)))) {
  SCOL=3+N
  NDIMS=length(grep("^p.x",colnames(res$manova_result)))
  FDRCOL=ncol(res$manova_result)-NDIMS+N
  y[[N]][[9]]<-res$manova_result[which(res$manova_result[,SCOL]>0 & res$manova_result[,FDRCOL]<0.05 ),1]
  y[[N]][[10]]<-res$manova_result[which(res$manova_result[,SCOL]<0 & res$manova_result[,FDRCOL]<0.05 ),1]
}

obs_up<-sapply(y,"[",9)
obs_dn<-sapply(y,"[",10)

gt_up<-sapply(y,"[",4)
gt_up<-lapply( gt_up , names)
gt_dn<-sapply(y,"[",5)
gt_dn<-lapply( gt_dn , names)

true_pos_up<-as.numeric(mapply( function(x,y) length(intersect(x,y)) , obs_up ,  gt_up ))
true_pos_dn<-as.numeric(mapply( function(x,y) length(intersect(x,y)) , obs_dn ,  gt_dn ))
false_pos_up<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , obs_up ,  gt_up ))
false_pos_dn<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , obs_dn , gt_dn ))
false_neg_up<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , gt_up ,  obs_up ))
false_neg_dn<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , gt_dn ,  obs_dn ))

true_pos<-mean(true_pos_up+true_pos_dn)
false_pos<-mean(false_pos_up+false_pos_dn)
false_neg<-mean(false_neg_up+false_neg_dn)
nrows<-as.numeric(lapply( sapply(y,"[",1 ), nrow))
true_neg<-mean(nrows-(true_pos+false_pos+false_neg))

p<-true_pos/(true_pos+false_pos)
r<-true_pos/(true_pos+false_neg)
f<-2*p*r/(p+r)

attr(y,'mitch_res') <-data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,true_pos,false_pos,true_neg,false_neg,p,r,f)
y
}


##################################
# hypergeometric test function (limited to deseq2)
##################################
run_hypergeometric<-function(x,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS){

dge<-sapply(x,"[",6)

ups<-lapply(dge, function(x) { rownames(subset(x, padj<0.05 & log2FoldChange > 0)) } )
dns<-lapply(dge, function(x) { rownames(subset(x, padj<0.05 & log2FoldChange < 0)) } )

l_ups<-sapply(ups,length)
l_dns<-sapply(dns,length)
geneset_sizes<-sapply( gsets , length )

# calculate number of genes in sets that are up and downregulated
n_dns=n_ups=p_ups=p_dns=obs_up=obs_dn=NULL
n_dns=n_ups=p_ups=p_dns=obs_up=obs_dn=list()

for (d in 1:length(dge)) {
  universe=length(rownames(dge[[d]]))
  n_ups[[d]]<-sapply( 1:length(gsets),function(x){length(which(gsets[[x]] %in% ups[[d]] ))} )
  n_dns[[d]]<-sapply( 1:length(gsets),function(x){length(which(gsets[[x]] %in% dns[[d]] ))} )

  p_ups[[d]]<-sapply( 1:length(gsets),function(x){phyper((n_ups[[d]][[x]]-1),l_ups[[d]],universe-geneset_sizes[[x]],geneset_sizes[[x]],lower.tail=FALSE,log.p=FALSE)})
  p_dns[[d]]<-sapply( 1:length(gsets),function(x){phyper((n_dns[[d]][[x]]-1),l_dns[[d]],universe-geneset_sizes[[x]],geneset_sizes[[x]],lower.tail=FALSE,log.p=FALSE)})

  x[[d]][[11]]<-names(gsets[which(p.adjust(p_ups[[d]],method="fdr")<0.05)])
  x[[d]][[12]]<-names(gsets[which(p.adjust(p_dns[[d]],method="fdr")<0.05)])
}

obs_up<-sapply(x,"[",11)
obs_dn<-sapply(x,"[",12)

gt_up<-sapply(x,"[",4)
gt_up<-lapply( gt_up , names)
gt_dn<-sapply(x,"[",5)
gt_dn<-lapply( gt_dn , names)

true_pos_up<-as.numeric(mapply( function(x,y) length(intersect(x,y)) , obs_up ,  gt_up ))
true_pos_dn<-as.numeric(mapply( function(x,y) length(intersect(x,y)) , obs_dn ,  gt_dn ))
false_pos_up<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , obs_up ,  gt_up ))
false_pos_dn<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , obs_dn , gt_dn ))
false_neg_up<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , gt_up ,  obs_up ))
false_neg_dn<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , gt_dn ,  obs_dn ))

true_pos<-mean(true_pos_up+true_pos_dn)
false_pos<-mean(false_pos_up+false_pos_dn)
false_neg<-mean(false_neg_up+false_neg_dn)
nrows<-as.numeric(lapply( sapply(x,"[",1 ), nrow))
true_neg<-mean(nrows-(true_pos+false_pos+false_neg))

p<-true_pos/(true_pos+false_pos)
r<-true_pos/(true_pos+false_neg)
f<-2*p*r/(p+r)

attr(x,'phyper_res') <-data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,true_pos,false_pos,true_neg,false_neg,p,r,f)
x

}


##################################
# aggregate function
##################################
agg_dge<-function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,gsets) {
# N_REPS=3 ; SUM_COUNT=40000000 ; VARIANCE=0.4 ; FRAC_DE=0.5 ; FC=1 ; SIMS=10 ; DGE_FUNC="deseq2" ; gsets=gsets
library("mitch")
xxx<-RepParallel(SIMS,simrna(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,gsets), simplify=F, mc.cores = detectCores() )

# run deseq2
xxx<-mclapply(xxx , DGE_FUNC , mc.cores = detectCores() )

# run mitch
xxx<-run_mitch(xxx,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS)

# run phyper
xxx<-run_hypergeometric(xxx,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS)

# return the result
g=list()
for (f in 1:length(attributes(xxx))) {
 PWAY_FUNC<-names(attributes(xxx)[f])
 PWAY_FUNC<-as.data.frame(PWAY_FUNC)
 g[[f]]<-cbind(unname(attributes(xxx)[f]),PWAY_FUNC)
}

g
}
# x<-agg_dge(a,10,40000000,0.4,0.2,1,10,"deseq2",gsets) 

