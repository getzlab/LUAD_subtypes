############################################################################################
############################################################################################
#### Copyright (c) 2017, Broad Institute, Inc. All rights reserved.
#### Redistribution and use in source and binary forms, with or without
#### modification, are permitted provided that the following conditions are
#### met:
####     Redistributions of source code must retain the above copyright
####     notice, this list of conditions and the following disclaimer.
####     Redistributions in binary form must reproduce the above copyright
####     notice, this list of conditions and the following disclaimer in
####     the documentation and/or other materials provided with the
####     distribution.
####     Neither the name of the Broad Institute nor the names of its
####     contributors may be used to endorse or promote products derived
####     from this software without specific prior written permission.
#### THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#### "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#### LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#### A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#### HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#### SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#### LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#### DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#### THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#### (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#### OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
############################################################################################
############################################################################################

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

