library("tidyverse")
library("parallel")
library("topconfects")
library("edgeR")
library("DESeq2")
library("limma")
library("ABSSeq")


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


# Get the reactome gene set GMT
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip",destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
GMT="ReactomePathways.gmt"
#' gmt_import
#'
#' This function imports GMT files into a list of character vectors for mitch analysis. GMT files are a commonly used
#' format for lists of genes used in pathway enrichment analysis. GMT files can be obtained from Reactome, MSigDB, etc.
#' @param gmtfile a gmt file.
#' @return a list of gene sets.
#' @keywords import genesets
#' @export
#' @examples
#' # Import some gene sets
#' genesets<-gmt_import("MyGeneSets.gmt")

gmt_import<-function(gmtfile){
    genesetLines <- strsplit(readLines(gmtfile), "\t")
    genesets <- lapply(genesetLines, utils::tail, -2)
    names(genesets) <- sapply(genesetLines, head, 1)
    attributes(genesets)$originfile<-gmtfile
    genesets
}

########################################
# simulate some gene expression data
########################################
simrna<-function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,GMT) {

# N_REPS=5 ; SUM_COUNT=10000000 ; VARIANCE=0.2 ; FRAC_DE=0.05 ; FC=1 ; GMT="ReactomePathways.gmt"

gset<-gmt_import(GMT)

library("edgeR")

df = NULL
for (k in paste0("data",1:(N_REPS*2)))  {
        b<-thinCounts(a,target.size=SUM_COUNT)
        colnames(b)=k
        df = cbind(df,b)
     }

# now need to only include gsets with 10 members in the 
gset_sub<-which(unlist( lapply(gset,function(x) { length(which(rownames(a) %in% as.character(unlist(x)))) >10 }  ) ) )
gset<-gset[which(names(gset) %in% names(gset_sub))]

#Number of differential genes
NDIF=round(length(gset)*FRAC_DE)

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
  DE_LIST<-sample(gset , NDIF)

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
xx <- list("x" = x,"UP_LIST"=UP_LIST,"DN_LIST"=DN_LIST,"UP_DE"=UP_DE,"DN_DE"=DN_DE)
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
z<- results(res)
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
agg_dge<-function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,GMT) {
#N_REPS=5 ; SUM_COUNT=10000000 ; VARIANCE=0.3 ; FRAC_DE=0.2 ; FC=1 ; SIMS=10 ; DGE_FUNC="edger"
xxx<-RepParallel(SIMS,simrna(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,GMT), simplify=F, mc.cores = detectCores() )
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



###############################################
#sanity checks to make sure the simulation is working as expected
###############################################
pdf(file="sanity_check.pdf")
x1<-simrna(a,5,10000000,0,0,0,GMT)
colnames(x1$x)<-paste("x1",colnames(x1$x),sep="_")
x2<-simrna(a,5,10000000,1,0,0,GMT)
colnames(x2$x)<-paste("x2",colnames(x2$x),sep="_")
xx<-merge(x1$x,x2$x,by=0)
heatmap(cor(xx[,2:ncol(xx)]),scale="none",cex.main=0.9,main="check1")

x1<-simrna(a,5,10000000,0,0,0,GMT)
colnames(x1$x)<-paste("x1",colnames(x1$x),sep="_")
x2<-simrna(a,5,10000000,0,0.1,1,GMT)
colnames(x2$x)<-paste("x2",colnames(x2$x),sep="_")
xx<-merge(x1$x,x2$x,by=0)
heatmap(cor(xx[,2:ncol(xx)]),scale="none",cex.main=0.9,main="check2")

x1<-simrna(a,5,10000000,0,0,0,GMT)
colnames(x1$x)<-paste("x1",colnames(x1$x),sep="_")
x2<-simrna(a,5,10000000,1,0.1,1,GMT)
colnames(x2$x)<-paste("x2",colnames(x2$x),sep="_")
xx<-merge(x1$x,x2$x,by=0)
heatmap(cor(xx[,2:ncol(xx)]),scale="none",cex.main=0.9,main="check3")
dev.off()


###############################################
# 10M reads with edger classic
###############################################
SIMS=10
for ( FRAC_DE in c(0.05)) {
  PDFNAME=paste(FRAC_DE,"_pw.pdf",sep="")
  pdf(file=PDFNAME,width=11.7,height=6.9)
  for (FC in c(1)) {
    par(mfrow=c(3,5))
    for (N_REPS in c(3,5,10)) {
      res=NULL
      for (DGE_FUNC in c("edger","edger_ql","deseq","limma","absseq")) {
        for ( SUM_COUNT in c(10000000,40000000,100000000)) {
          for  ( VARIANCE in c(0,0.2,0.3,0.4,0.5)) {
            res_new<-agg_dge(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,GMT)
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

        points(res1_1e8_v0$r,res1_1e8_v0$p,xlab="recall",ylab="precision",pch=15,col="gray",xlim=c(0,1),ylim=c(0,1))
        points(res1_1e8_v2$r,res1_1e8_v2$p,xlab="recall",ylab="precision",pch=16,col="gray",xlim=c(0,1),ylim=c(0,1))
        points(res1_1e8_v3$r,res1_1e8_v3$p,xlab="recall",ylab="precision",pch=17,col="gray",xlim=c(0,1),ylim=c(0,1))
        points(res1_1e8_v4$r,res1_1e8_v4$p,xlab="recall",ylab="precision",pch=0,col="gray",xlim=c(0,1),ylim=c(0,1))
        points(res1_1e8_v5$r,res1_1e8_v5$p,xlab="recall",ylab="precision",pch=1,col="gray",xlim=c(0,1),ylim=c(0,1))

        legend(0.4,0.2,legend=c("10M","40M","100M"),col=c("red", "blue","gray") ,pch=19,cex=0.6,title="read depth")
        legend(0.75,0.2,legend=c("0","0.2","0.3","0.4","0.5"),col=c("gray") ,pch=c(15:17,0,1),cex=0.6,title="added variance")
        mtext(paste(DGE_FUNC,N_REPS,"reps, 10% DEG") ,cex=0.8); grid()

      }
    }
  }
  dev.off()
}


q()



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
res_absseq<-res[which(res$DGE_FUNC=="absseq"),]
res_absseq<-res_absseq[which(res_absseq$FRAC_DE==0.1),]

pdf(file="barplots.pdf") ; par(mfrow=c(1,5)) ; par(mar=c(4,8,3,2))

barplot(prop.table(as.matrix(t(apply(res_edger[,8:11],2,rev))),2),las=1,horiz=T,xlab="Gene proportion", 
 cex.names=0.5,cex.main=0.9,main="edgeR 5 reps varied variance and DEG%",
 legend.text = TRUE,args.legend=list(x="topright",bty="n",inset=c(-0.05, 0)))

barplot(prop.table(as.matrix(t(apply(res_edgerql[,8:11],2,rev))),2),las=1,horiz=T,xlab="Gene proportion", 
 cex.names=0.5,cex.main=0.9,main="edgeR QL 5 reps varied variance and DEG%",
 legend.text = F,args.legend=list(x="topright",bty="n",inset=c(-0.05, 0)))

barplot(prop.table(as.matrix(t(apply(res_deseq[,8:11],2,rev))),2),las=1,horiz=T,xlab="Gene proportion", 
 cex.names=0.5,cex.main=0.9,main="DESeq2 5 reps varied variance and DEG%",
 legend.text = F,args.legend=list(x="topright",bty="n",inset=c(-0.05, 0)))

barplot(prop.table(as.matrix(t(apply(res_limma[,8:11],2,rev))),2),las=1,horiz=T,xlab="Gene proportion", 
 cex.names=0.5,cex.main=0.9,main="Limma 5 reps varied variance and DEG%",
 legend.text = F,args.legend=list(x="topright",bty="n",inset=c(-0.05, 0)))

barplot(prop.table(as.matrix(t(apply(res_absseq[,8:11],2,rev))),2),las=1,horiz=T,xlab="Gene proportion",
 cex.names=0.5,cex.main=0.9,main="ABSseq 5 reps varied variance and DEG%",
 legend.text = F,args.legend=list(x="topright",bty="n",inset=c(-0.05, 0)))


barplot(prop.table(as.matrix(t(apply(smm_0de[,1:4],2,rev))),2),las=1,
 ylab="Added noise",horiz=T,xlab="Gene proportion",cex.main=0.9,main="edger 10M reads no DGEs",
 legend.text = TRUE,args.legend=list(x="topright",bty="n",inset=c(-0.05, 0)))

barplot(prop.table(as.matrix(t(apply(smm_1de[,1:4],2,rev))),2),las=1,
 ylab="Added noise",horiz=T,xlab="Gene proportion",cex.main=0.9,main="edger 10M reads 5% DGEs",
 legend.text = FALSE,args.legend=list(x="topright",bty="n",inset=c(-0.05, 0)))

dotchart(rev(smm_0de$f),pch=19,labels=rev(rownames(smm)),xlim=c(0,1),cex.main=0.9,main="edger 10M reads no DGEs")
legend("topright", pch = c(19,19,19), col = c("gray", "red", "blue"),legend = c("F1","p","r"))
grid()
points(rev(smm_0de$p),1:5,pch=19,col="red")
points(rev(smm_0de$r),1:5,pch=19,col="blue")
points(rev(smm_0de$f),1:5,pch=19,col="gray")

dotchart(rev(smm_1de$f),pch=19,labels=rev(rownames(smm)),xlim=c(0,1),cex.main=0.9,main="edger 10M reads 5% DGEs")
#legend("topright", pch = c(19,19,19), col = c("gray", "red", "blue"),legend = c("F1","p","r"))
grid()
points(rev(smm_1de$p),1:5,pch=19,col="red")
points(rev(smm_1de$r),1:5,pch=19,col="blue")
points(rev(smm_1de$f),1:5,pch=19,col="gray")

