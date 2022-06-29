NMF.W <- function(X,W,tol,K) {
  n.run <- 1
  n.iter <- 1000000
  eps <- 1.e-50
  N <- dim(X)[1]
  M <- dim(X)[2]
  meanX <- mean(X,na.rm=T)
  eps <- 1.e-50
  for (j in 1:n.run) {
    H <- matrix(runif(K * M)*meanX,ncol=M)
    X.ap <- W %*% H
    error.EU <- sum((X-X.ap)^2)
    error.KL <- sum(X*log((X+eps)/(X.ap+eps))+X.ap-X)
    del <- 1
    count <- 1
    while (del >= tol & count < n.iter) {
      H <- H * (t(W) %*% X) / (t(W)%*%(W%*%H) + eps)
      X.ap <- W %*% (H)
      del <- abs(error.EU-sum((X-X.ap)^2))
      error.EU <- sum((X-X.ap)^2)
      error.KL <- sum(X*log((X+eps)/(X.ap+eps))+X.ap-X)
      if (count %% 100 == 0) cat(count,error.EU,error.KL,del,'\n')
      count <- count+1
    }
  }
  return(list(W,H))
}

NMF.H <- function(X,H,tol,K) {
  n.run <- 1
  n.iter <- 1000000
  eps <- 1.e-50
  N <- dim(X)[1]
  M <- dim(X)[2]
  meanX <- mean(X,na.rm=T)
  eps <- 1.e-50
  for (j in 1:n.run) {
    W <- matrix(runif(N * K)*meanX,ncol=K)
    X.ap <- W %*% H
    error.EU <- sum((X-X.ap)^2)
    error.KL <- sum(X*log((X+eps)/(X.ap+eps))+X.ap-X)
    del <- 1
    count <- 1
    while (del >= tol & count < n.iter) {
      W <- W * (X %*% t(H)) / ((W%*%H) %*% t(H)+eps)
      X.ap <- W %*% H
      del <- abs(error.EU-sum((X-X.ap)^2))
      error.EU <- sum((X-X.ap)^2)
      if (count %% 100 == 0) cat(count,error.EU,del,'\n')
      count <- count+1
    }
  }
  return(list(W,H))
}

get.SSEC.fold.short <- function(expr,expr.fold,expr.fold.up,cohort,marker0,W1) {
  n.sample <- ncol(expr)
  comm <- intersect(rownames(expr),marker0)
  W0.tmp <- as.matrix(W1[match(comm,rownames(W1),nomatch=0),])
  W0.tmp.norm <- t(apply(W0.tmp,1,function(x) x/sum(x)))
  K <- ncol(W1)
  H.tmp <- array(0,dim=c(K,n.sample))
  for (i in 1:n.sample) {
    X0.tmp <- as.matrix(expr.fold.up[match(comm,rownames(expr.fold.up),nomatch=0),i])
    x <- get.single.NMF(X0.tmp,W0.tmp,1.e-07,K)
    H.tmp[,i] <- x[[1]]
  }
  rownames(H.tmp) <- colnames(W0.tmp)
  colnames(H.tmp) <- colnames(expr.fold)
  H.tmp.norm <- apply(H.tmp,2,function(x) x/sum(x))
  g.tmp <- apply(H.tmp.norm,2,function(x) which.max(x))
  #rownames(H.tmp.norm) <- c("Luminal","Luminal-infiltrated","Basal-squamous","Neuronal","Luminal-papiilary")
  rownames(H.tmp.norm) <- paste("G",seq(1:ncol(W1)),sep="")
  x <- list(g.tmp,H.tmp)
  return(x)
}

get.single.NMF <- function(X1,W0,tol,K) {
  X1[is.na(X1)] <- 0
  res <- NMF.W(X1,W0,tol,K)
  H1 <- res[[2]]
  g1 <- apply(H1,2,function(x) which.max(x))
  return(list(H1,g1))
}

get.marker.fold <- function(W1.sig,W1.norm.sig,cut.W,cut.fold,expression) {
  
  mRNA.fold <- expression
  index.max <- apply(W1.norm.sig,1,function(x) which.max(x))
  
  summary <- list()
  
  for(i in 1:ncol(W1.sig)){ # the total number of clusters

    eval(parse(text=paste0("gene",i," <- unique(names(index.max)[index.max==",i,"])")))
    eval(parse(text=paste0("index.sample",i," <- colnames(mRNA.fold)%in%sample",i)))
    eval(parse(text=paste0("index.gene",i," <- rownames(mRNA.fold)%in%gene",i)))
    eval(parse(text=paste0("index.sig",i," <- match(gene",i,",rownames(W1.sig),nomatch=0)")))
    eval(parse(text=paste0("summary.gene",i," <- data.frame(gene",i,",W1.sig[index.sig",i,",",i,"],W1.norm.sig[index.sig",i,",],rowMeans(mRNA.fold[index.gene",i,",index.sample",i,"],na.rm=T),rowMeans(mRNA.fold[index.gene",i,",!index.sample",i,"],na.rm=T))")))
    
    temp <- ""
    for(j in 1:ncol(W1.sig)){
      temp <- paste0(temp,"\"W.norm.sig",j,"\",")
    }
    eval(parse(text=paste0("colnames(summary.gene",i,") <- c(\"gene\",\"W.sig\",",temp,"\"mean1\",\"mean2\"",")")))
    
    eval(parse(text=paste0("summary.gene",i,"[,\"fold\"] <- summary.gene",i,"$mean1-summary.gene",i,"$mean2")))
    eval(parse(text=paste0("summary.gene",i," <- summary.gene",i,"[order(summary.gene",i,"$fold,decreasing=T),]")))
    eval(parse(text=paste0("gene",i," <- summary.gene",i,"$gene[(summary.gene",i,"$mean1-summary.gene",i,"$mean2)>cut.fold[",i,"]]")))
    eval(parse(text=paste0("summary[[",i,"]] <- summary.gene",i)))
    
  }
  
  temp <- ""
  for(j in 1:ncol(W1.sig)){
    temp <- paste0(temp,"gene",j,",")
  }
  eval(parse(text=paste0("return(list(",temp,"summary))")))
  
}

########## loading expression data
##########

mRNA.comm <- TUMOR.GEXP.log2.mat

index.gene <- rowSums(is.na(mRNA.comm))<round(cut.NA*ncol(mRNA.comm))
mRNA.comm <- mRNA.comm[index.gene,]
mRNA.fold <- t(apply(mRNA.comm,1,function(x) x-median(x,na.rm=T)))
mRNA.fold1 <- mRNA.fold
mRNA.fold2 <- mRNA.fold
mRNA.fold1[is.na(mRNA.fold1)] <- 0
mRNA.fold2[is.na(mRNA.fold2)] <- 0
mRNA.fold1 <- apply(mRNA.fold1,2,function(x) ifelse(x>0,x,0))
mRNA.fold2 <- apply(mRNA.fold2,2,function(x) ifelse(x>0,0,-x))
rownames(mRNA.fold1) <- paste(rownames(mRNA.fold1),"up",sep=".")
rownames(mRNA.fold2) <- paste(rownames(mRNA.fold2),"dn",sep=".")
mRNA.fold0 <- rbind(mRNA.fold1,mRNA.fold2)

########## loading clustering data
######################

# Retrieve df information  

OUTPUT_TPM <- paste0("/Users/wroh/Google Drive/Work/Projects/Expression_Subtyping/data/GDAN_TMP/BayesNMF/TMP_Tier3/OUTPUT_",tumor.type,"/")
#OUTPUT_TPM <- paste0("/Users/wroh/Google Drive/Work/Projects/Expression_Subtyping/data/BayesNMF/TPM_No_sizefactor_norm/OUTPUT_",tumor.type,"/")
#OUTPUT_TPM <- paste0("/Users/wroh/Google Drive/Work/Projects/Expression_Subtyping/data/BayesNMF/RSEM/OUTPUT_",tumor.type,"/")

n.iter <- 50 ## # of iterations 

tmpK <- rep(0,n.iter)
tmpE <- rep(0,n.iter)
for (i in 1:n.iter) {
  load(file=paste(OUTPUT_TPM,paste("res.L1EU.Bayes",i,"RData",sep="."),sep=""))
  lambda <- res.Bayes[[5]]
  lambda <- unlist(lambda[length(lambda)])
  lambda <- lambda-min(lambda)
  cat(lambda,sum(lambda!=0),'\n')
  tmpK[i] <- sum(lambda > 0)
  tmpE[i] <- res.Bayes[[4]][length(res.Bayes[[4]])]
}

################ summary of BayesNMF runs
#### df has all info for BayesNMF runs and please choose the run # with the lowest "evid" for given K (== runK)
df <- data.frame(seq(1:n.iter),tmpK,unlist(tmpE))
colnames(df) <- c("run","K","evid")
df <- df[order(df$evid,decreasing=T),]

run.K = min(subset(df,evid==min(subset(df,K==max(as.numeric(names(table(df$K))[table(df$K)==max(table(df$K))])))$evid))$run)
#run.K = 28  # manual selection of run.K (for K=5 in BCLA RSEM)

load(file=paste(OUTPUT_TPM,paste("res.L1EU.Bayes",run.K,"RData",sep="."),sep=""))

#load(file=paste(paste("res.L1EU.Bayes.2.RData",sep="."),sep="")) #### BayesNMF ouput used for the de-novo expression subtyping
res <- res.Bayes
W <- res[[1]]
H <- res[[2]]
W <- W[,colSums(W)!=0]
H <- H[rowSums(H)!=0,]
rownames(H) <- paste("G",seq(1:nrow(H)),sep="")
H.norm <- apply(H,2,function(x) x/sum(x))
g.Bayes <- apply(H.norm,2,function(x) which.max(x))

for(i in 1:length(unique(g.Bayes))){
  eval(parse(text=paste0("sample",i," <- names(g.Bayes)[g.Bayes==",i,"]")))
}

# sample1 <- names(g.Bayes)[g.Bayes==1]
# sample2 <- names(g.Bayes)[g.Bayes==2]
# sample3 <- names(g.Bayes)[g.Bayes==3]
# sample4 <- names(g.Bayes)[g.Bayes==4]
# sample5 <- names(g.Bayes)[g.Bayes==5]

############ determine W for tumor.type fold expression and H.norm

H1 <- H
H1.norm <- H.norm

set.seed(10)
if (!file.exists(paste0(OUTPUT_TPM,paste("res.NMF2",cut.NA,"RData",sep=".")))) {
  res.NMF2 <- NMF.H(mRNA.fold0,H1.norm,1.e-05,nrow(H1))
  save(res.NMF2,file=paste0(OUTPUT_TPM,paste("res.NMF2",cut.NA,"RData",sep=".")))
} else {
  load(file=paste0(OUTPUT_TPM,paste("res.NMF2",cut.NA,"RData",sep=".")))
}

W1 <- res.NMF2[[1]]
W1 <- W1[index.gene,]
W1.up <- W1[grep("up",rownames(W1)),]
W1.dn <- W1[grep("dn",rownames(W1)),]
W1.norm <- t(apply(W1,1,function(x) x/sum(x)))
W1.norm.up <-  W1.norm[grep("up",rownames(W1.norm)),]
W1.norm.dn <-  W1.norm[grep("dn",rownames(W1.norm)),]

############ compute the maximum association to clusters across genes
############ up-regulated genes
max.up <- apply(W1.norm.up,1,function(x) max(x))
cut.W1 <- quantile(max.up,prob=cut.W,na.rm=T)
gene.sig.up <- rownames(W1.norm.up)[max.up >= cut.W1]
W1.sig.up <- W1.up[match(gene.sig.up,rownames(W1.up),nomatch=0),]
W1.norm.sig.up <- W1.norm.up[match(gene.sig.up,rownames(W1.norm.up),nomatch=0),]
x <- get.marker.fold(W1.sig.up,W1.norm.sig.up,cut.W,cut.fold,mRNA.fold1) 

for(i in 1:ncol(W1.sig.up)){
  eval(parse(text=paste0("gene",i,".up <- x[[",i,"]]")))
}
summary.gene.up <- x[[ncol(W1.sig.up)+1]]

# gene1.up <- x[[1]]
# gene2.up <- x[[2]]
# gene3.up <- x[[3]]
# gene4.up <- x[[4]]
# gene5.up <- x[[5]]
# summary.gene.up <- x[[6]]

############ down-regulated genes, we will not use this information
# cut.W1 <- quantile(W1.norm.dn,prob=cut.W,na.rm=T)
# W1.norm.index <- W1.norm.dn > cut.W1
# gene.sig.dn <- unique(c(rownames(W1.norm.dn)[rowSums(W1.norm.index)>0])) #"PEG10","PADI3"))
# W1.sig.dn <- W1.dn[match(gene.sig.dn,rownames(W1.dn),nomatch=0),]
# W1.norm.sig.dn <- W1.norm.dn[match(gene.sig.dn,rownames(W1.norm.dn),nomatch=0),]
# x <- get.marker.fold(W1.sig.dn,W1.norm.sig.dn,cut.W,cut.fold,mRNA.fold2)
# gene1.dn <- x[[1]]
# gene2.dn <- x[[2]]
# gene3.dn <- x[[3]]
# gene4.dn <- x[[4]]
# gene5.dn <- x[[5]]
# summary.gene.dn <- x[[6]]

for(i in 1:ncol(W1.sig.up)){
  eval(parse(text=paste0("gene",i," <- gsub(\".up\",\"\",gene",i,".up)")))
  eval(parse(text=paste0("if (length(gene",i,")>cut.gene) gene",i," <- gene",i,"[1:cut.gene]")))
  eval(parse(text=paste0("gene",i," <- gene",i,"[gene",i,"%in%rownames(mRNA.comm)]")))
}

########### adding manual genes which you may interestrd in
# if (EXTRA) {
#   gene1.extra <- c("PPARG","GATA3","FOXA1")
#   gene2.extra <- c("ZEB1","ZEB2","SNAI1","TWIST1","CDH2","CLDN3","CLDN4","CLDN7")
#   gene2.extra <- c("ZEB1","ZEB2","SNAI1","TWIST1","CHD2")
#   gene3.extra <- c("CD44")
#   gene4.extra <- c("CHGA","SCG2","ENO2","SYP")
#   gene5.extra <- c("FGFR3","SHH","BMP5")
#   gene1 <- unique(c(as.character(gene1),gene1.extra))
#   gene2 <- unique(c(as.character(gene2),gene2.extra))
#   gene3 <- unique(c(as.character(gene3),gene3.extra))
#   gene4 <- unique(c(as.character(gene4),gene4.extra))
#   gene5 <- unique(c(as.character(gene5),gene5.extra))
# }
# 

############################
############################
W1 <- W1.up
rownames(W1) <- gsub(".up","",rownames(W1))

for(i in 1:ncol(W1.up)){
  eval(parse(text=paste0("W.marker",i," <- W1[match(as.character(gene",i,"),rownames(W1),nomatch=0),]")))
  eval(parse(text=paste0("W.marker",i,".norm <- t(apply(W.marker",i,",1,function(x) x/sum(x)))")))
  eval(parse(text=paste0("gene",i," <- gene",i,"[order(W.marker",i,".norm[,",i,"],decreasing=T)]")))
}

####

# W.marker1 <- W1[match(as.character(gene1),rownames(W1),nomatch=0),]
# W.marker2 <- W1[match(as.character(gene2),rownames(W1),nomatch=0),]
# W.marker3 <- W1[match(as.character(gene3),rownames(W1),nomatch=0),]
# W.marker4 <- W1[match(as.character(gene4),rownames(W1),nomatch=0),]
# W.marker5 <- W1[match(as.character(gene5),rownames(W1),nomatch=0),]
# W.marker1.norm <- t(apply(W.marker1,1,function(x) x/sum(x)))
# W.marker2.norm <- t(apply(W.marker2,1,function(x) x/sum(x)))
# W.marker3.norm <- t(apply(W.marker3,1,function(x) x/sum(x)))
# W.marker4.norm <- t(apply(W.marker4,1,function(x) x/sum(x)))
# W.marker5.norm <- t(apply(W.marker5,1,function(x) x/sum(x)))
# gene1 <- gene1[order(W.marker1.norm[,1],decreasing=T)]
# gene2 <- gene2[order(W.marker2.norm[,2],decreasing=T)]
# gene3 <- gene3[order(W.marker3.norm[,3],decreasing=T)]
# gene4 <- gene4[order(W.marker4.norm[,4],decreasing=T)]
# gene5 <- gene5[order(W.marker5.norm[,5],decreasing=T)]

temp <- "marker0 <- unique(c("
for(j in 1:ncol(W1.sig.up)){
  ifelse(j!=ncol(W1.sig.up), temp <- paste0(temp,"as.character(gene",j,"),"),
                             temp <- paste0(temp,"as.character(gene",j,")))")      )
}
eval(parse(text=temp))
             
#marker0 <- unique(c(as.character(gene1),as.character(gene2),as.character(gene3),as.character(gene4),as.character(gene5)))
W.marker0 <- W1[match(marker0,rownames(W1),nomatch=0),]
W.marker0.norm <- t(apply(W.marker0,1,function(x) x/sum(x)))