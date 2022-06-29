#### Projection to CPTAC LUAD

# set a working directory
setwd("./src/Subtype projection codes")

#############################$
# Load CPTAC LUAD TPM data
#############################$

CPTAC.LUAD <- readRDS(file = "CPTAC.LUAD.rds") # Replace this with the expression matrix of interest for subtype projection

######### projection to CPTAC LUAD

# Load TCGA LUAD expression data
TUMOR.GEXP.log2.mat <- readRDS(file = "TUMOR.GEXP.log2.mat.rds")

OUTPUT_CPTAC <- "./input_data/OUTPUT_CPTAC/"

threshold.NA <- 0.1
cut.bottom <- 0.1 ## filter out bottom 10% genes in terms of mean expression
cut.sd <- 0.75 ## include genes with top 25% in terms of SD
cut.gene <- 100 ### number of markers in each subtype is restricted by cut.gene
cut.W <- 0.50 ### only consider genes with the normalized association to clusters >= 0.50
cut.fold <- rep(0.50,length(unique(g.Bayes))) ### select markers with mean difference of log2(fold changes) >= 0.50
cut.NA <- 0.1 ### remove genes with NA values across samples > 0.1 

tumor.type <- "LUAD"

source("./src/get.marker_for_fold.ver2.WR.TMP.R")
marker0 <- sapply(marker0,function(x) strsplit(x,"\\|")[[1]][1])
marker0.original <- marker0

cohort <- "CPTAC-LUAD-Tumor"
expr <- CPTAC.LUAD

index.NA <- rowSums(is.na(expr))<= round(cut.NA*ncol(expr))
expr <- expr[index.NA,]
expr.fold <- t(apply(expr,1,function(x) x-median(x,na.rm=T)))
expr.fold[is.na(expr.fold)] <- 0
expr.fold.up <- apply(expr.fold,2,function(x) ifelse(x>0,x,0))
rownames(W1) <- sapply(rownames(W1),function(x) strsplit(x,"[|]")[[1]][1])
for (i in 1:1) {
  gene.expr <- sapply(rownames(expr.fold),function(x) strsplit(x,"_")[[1]][2])
  marker0 <- rownames(expr.fold)[gene.expr%in%marker0.original] 
  gene.marker0 <- sapply(marker0,function(x) strsplit(x,"_")[[1]][2]) 
  W1.new <- W1[match(gene.marker0,rownames(W1),nomatch=0),]
  rownames(W1.new) <- marker0
  expr.fold.marker <- expr.fold[match(marker0,rownames(expr.fold),nomatch=0),]
  x <- get.SSEC.fold.short(expr,expr.fold,expr.fold.up,cohort,marker0,W1.new)
  g.Bayes.new <- x[[1]]
  H.Bayes.new <- x[[2]]
  colnames(H.Bayes.new) <- names(g.Bayes.new)
  save(H.Bayes.new,file=paste(OUTPUT_CPTAC,paste(cohort,"H.Bayes.RData",sep="."),sep=""))
}

# Normalized H matrix
H.Bayes.new.norm <- apply(H.Bayes.new,2,function(x) x/sum(x))

