library(gplots)
library(RColorBrewer)
library(clusterSim)
library(ggplot2)
library(reshape)
library(reshape2)

get.consensus.clustering <- function(res,consensus.norm) {
        W <- res[[1]]
        H <- res[[2]]
        W <- W[,colSums(W)!=0]
        H <- H[rowSums(H)!=0,]
        rownames(H) <- paste("G",seq(1:nrow(H)),sep="")
        H.norm <- apply(H,2,function(x) x/sum(x))
        g.Bayes <- apply(H.norm,2,function(x) which.max(x))
        K0 <- nrow(H.norm)
        #h.Bayes.HC <- hclust(dist(t(H.norm),method="euclidean"),method=finalLinkage)
        #g.Bayes.HC <- cutree(h.Bayes.HC,K0)

        #corr.consensus <- get.correlation.spearman(consensus.norm)
        #dt.consensus <- as.dist(1-corr.consensus)
        #h.HC <- hclust(dt.consensus,method=finalLinkage)
        #g.HC <- cutree(h.HC,K0)

        #s.Bayes <- silhouette(g.Bayes,dt.consensus)
        #s.Bayes.HC <- silhouette(g.Bayes.HC,dt.consensus)
        #s.HC <- silhouette(g.HC,dt.consensus)
        #x <- list(g.Bayes,h.Bayes.HC,g.Bayes.HC,s.Bayes.HC,h.HC,g.HC,s.HC,s.Bayes)
	x <- list(g.Bayes)
        return(x)
}

get.correlation.pearson <- function(mat) {
        n.gene <- nrow(mat)
        n.sample <- ncol(mat)
        corr <- array(0,dim=c(n.sample,n.sample))
        for (i in 1:n.sample) {
        for (j in i:n.sample) {
                corr[i,j] <- cor.test(mat[,i],mat[,j],method="pearson")$estimate
                corr[j,i] <- corr[i,j]
        }
        }
        return(corr)
}

get.correlation.spearman <- function(mat) {
        n.gene <- nrow(mat)
        n.sample <- ncol(mat)
        corr <- array(0,dim=c(n.sample,n.sample))
        for (i in 1:n.sample) {
        for (j in i:n.sample) {
                corr[i,j] <- cor.test(mat[,i],mat[,j],method="spearman")$estimate
                corr[j,i] <- corr[i,j]
        }
        cat(i,'\n')
        }
        return(corr)
}

BayesNMF.L1EU <- function(V0,n.iter,a0,tol,K,K0,phi) {
        eps <- 1.e-50
        del <- 1.0
        active_nodes <- colSums(V0) != 0
        V0 <- V0[,active_nodes] 
        V <- V0-min(V0) 
        Vmin <- min(V)
        Vmax <- max(V)
        N <- dim(V)[1]
        M <- dim(V)[2]

        W <- matrix(runif(N * K)*Vmax,ncol=K)
        H <- matrix(runif(M * K)*Vmax,ncol=M)
        V.ap <- W%*%H+eps
        I <- array(1,dim=c(N,M))

        phi <- sd(V)^2*phi
        C <- N+M+a0+1
        b0 <- sqrt((a0-1)*(a0-2)*mean(V,na.rm=T)/K0)
        lambda.bound <- b0/C
        lambda <- (colSums(W)+rowSums(H)+b0)/C
        lambda.cut <- 1.5*lambda.bound

        n.like <- list()
        n.evid <- list()
        n.error <- list()
        n.lambda <- list()
        n.lambda[[1]] <- lambda
        iter <- 2
        count <- 1
        while (del >= tol & iter < n.iter) {
                H <- H * (t(W) %*% V) / (t(W) %*% V.ap + phi * matrix(rep(1/lambda,M),ncol=M) + eps)
                V.ap <- W %*% H + eps
                W <- W * (V %*% t(H)) / (V.ap %*% t(H) + phi * t(matrix(rep(1/lambda,N),ncol=N)) + eps)
                V.ap <- W %*%H + eps
                lambda <- (colSums(W)+rowSums(H)+b0)/C
                del <- max(abs(lambda-n.lambda[[iter-1]])/n.lambda[[iter-1]])
                like <- sum((V-V.ap)^2)/2
                n.like[[iter]] <- like
                n.evid[[iter]] <- like + phi*sum((colSums(W)+rowSums(H)+b0)/lambda+C*log(lambda))
                n.lambda[[iter]] <- lambda
                n.error[[iter]] <- sum((V-V.ap)^2)
                if (iter %% 100 == 0) {
                        cat(iter,n.evid[[iter]],n.like[[iter]],n.error[[iter]],del,sum(colSums(W)!=0),sum(lambda>=lambda.cut),'\n')
                }
                iter <- iter+1
        }
        return(list(W,H,n.like,n.evid,n.lambda,n.error))
}

plot.heatmap.3 <- function(x,rowTF,colTF,main) {
        s1 <- 0.75
        s2 <- 1.0
        s3 <- 1.5
        mydist <- function(c) {dist(c,method="euclidean")}
        myclust <- function(c) {hclust(c,method="ward.D")}
        heatmap.2(as.matrix(x), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both",margins=c(5,5),
                Rowv=rowTF, Colv=colTF, symbreaks=F, key=TRUE, symkey=F,main=main,
                density.info="none", trace="none",labCol=colnames(x),labRow=rownames(x),col=greenred(40),cex.lab=s1,cexRow=0.50,cexCol=0.50,keysize=s1)
}

get.consensus.mat <- function(corr,max.K,pItem,innerLinkage,finalLinkage,n.base) {
        sample <- rownames(corr)
        n.sample <- ncol(corr)
        n.resample <- round(pItem*ncol(corr))
        consensus <- array(0,dim=c(n.sample,n.sample))
        n.count <- 0
        for (K in 2:max.K) {
        n.iter <- K*n.base
        for (i in 1:n.iter) {
                index <- sample.int(n.sample,n.resample)
                corr.resample <- corr[index,index]
                dt <- as.dist(1-corr.resample)
                h.inner <- hclust(dt,method=innerLinkage)
                g <- cutree(h.inner,K)
                for (j in 1:K) {
                        id <- match(names(g)[g==j],sample,nomatch=0)
                        consensus[id,id] <- consensus[id,id]+1
                }
                n.count <- n.count+1
        }
        cat(K,n.count,'\n')
        }
        colnames(consensus) <- rownames(corr)
        rownames(consensus) <- rownames(corr)
        consensus.norm <- consensus/n.count
        return(list(consensus,consensus.norm))
}

scale <- 0.8
.theme_ss <- theme_bw(base_size=14) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=10*scale, family="mono"),
        axis.text.y = element_text(hjust = 0.5,size=10*scale, family="mono"),
        axis.title.x = element_text(face="bold",colour="black",size=14*scale),
        axis.title.y = element_text(face="bold",colour="black",size=14*scale),
        axis.text = element_text(size = 16*scale, family = "mono"),
        strip.text = element_text(lineheight=0.5),
        strip.text.x = element_text(size=10*scale,face='bold',angle=00),
        strip.text.y = element_text(size=10*scale,face="bold"),
        strip.background = element_rect(colour="black",fill="gray85"),
        panel.margin = unit(0.20,"lines"),
        plot.title=element_text(lineheight=1.0,face="bold",size=12*scale))

get.sample.association.heatmap <- function(H,g.Bayes,scale0) {
        #g.ordering <- c("G5","G4","G3","G2","G1")
        commands <- c("g.ordering <- c(")
        for(i in nrow(H):1){
          if(i != 1){
            commands <- paste0(commands,"\"G",i,"\",")
          }
          else {
            commands <- paste0(commands,"\"G",i,"\")")
          }
        }
        eval(parse(text=commands))
  
        sample.ordering <- colnames(H)[order(g.Bayes,decreasing=F)]
        df <- t(H)
        df <- melt(df)
        colnames(df) <- c("sample","cluster","value")
        df$sample <- factor(df$sample,sample.ordering)
        df$cluster <- factor(df$cluster,g.ordering)
        df1 <- df
        df1[,"type"] <- "H matrix"

        H.norm <- apply(H,2,function(x) x/sum(x))
        df <- t(H.norm)
        df <- melt(df) 
        colnames(df) <- c("sample","cluster","value")
        df$sample <- factor(df$sample,sample.ordering)
        df$cluster <- factor(df$cluster,g.ordering)
        df2 <- df
        df2[,"type"] <- "Normalized H"

        p = ggplot(df1,aes(x=sample,y=cluster,fill=value))+geom_tile() #geom_tile(colour="yellow")
        #p = p + facet_grid(. ~ type, scale = "free_y")
        p = p + scale_fill_gradient2(low="blue",mid="white",high ="darkblue",name=paste("Activity",sep=""))
        p = p + .theme_ss
        p = p + ggtitle("H matrix")
        p = p + xlab("Sample") + ylab("mRNA Clusters")
        p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=14*scale))
        p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=14*scale))
        p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=5*scale, family="mono",face='bold',color=color.axis))
        p1 = p + theme(legend.position="top")

        p = ggplot(df2,aes(x=sample,y=cluster,fill=value))+geom_tile() #geom_tile(colour="yellow")
        #p = p + facet_grid(. ~ type, scale = "free_y")
        p = p + scale_fill_gradient2(low="blue",mid="white",high ="darkblue",name=paste("Activity",sep=""))
        p = p + .theme_ss
        p = p + ggtitle("Normalized H matrix")
        p = p + xlab("Sample") + ylab("mRNA Clusters")
        p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=14*scale))
        p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=14*scale))
        p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=5*scale, family="mono",face='bold',color=color.axis))
        p2 = p + theme(legend.position="top")

        df <- data.frame(t(H))
        colnames(df) <- paste("H",seq(1:nrow(H)),sep="")
        df[,"mRNA"] <- g.Bayes
        df$mRNA <- paste("G",df$mRNA,sep="")
        df <- melt(df,id="mRNA")
        colnames(df) <- c("Subtype","cluster","Association")
        p = ggplot(df,aes(x=Association)) 
        p = p + geom_histogram(color="black",fill="gray75",binwidth=0.05) 
        p = p + facet_grid(Subtype ~ cluster,scale='free_y')
        p = p + .theme_ss
        p = p + ggtitle("Sample Associations")
        p = p + xlab("Association") + ylab("Sample Counts")
        p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=14*scale))
        p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=14*scale))
        p = p + theme(axis.text.x = element_text(angle=0,vjust=0.0, size=12*scale, family="mono",face='bold',color=color.axis))
        p3 = p + theme(legend.position="top")
        return(list(p1,p2,p3))
}

