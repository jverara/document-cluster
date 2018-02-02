ICD10_cluster <- function(filename, H, reduce.dim.lsa=FALSE,reduce.dim.pca=FALSE, rows=FALSE){
 
library(readxl)
library(proxy)  
library(xlsx)  
   
  #Read in data 
  if(!rows[1]){
    data  <- read_xlsx(filename)[,] #1:160
  }else{
    data  <- read_xlsx(filename)[rows,]
  }
  #return(data)
  
  # In depends on the header -> needs to be checked 
  head1 <- names(data)[5] #SNPs
  head2 <- names(data)[4] #diagnose 
  
  #if(rows[1] & is.element(1,rows)){
  #  head1 <- NULL
  #  head2 <- NULL
  #}
  
  #return(data) 
  ll   <- dim(data)[1]
 
#Get the diagnosis (documents)  
diag <- data[,4]
diag <- as.character(diag)
diag <- strsplit(diag,",")[[1]]
#concatenate very first elements
diag <- c(head2, diag)
SNPs <- c(head1,data[,5][[1]])
#return(diag)
#Cut the strings into snippets 
SL <- vector("list",length(SNPs))  
for(xx in 1:length(SL)){
  str      <- as.character(SNPs[xx])
  str      <- strsplit(str," ")
  SL[[xx]] <- str
}
  
Rw_names      <- unique(unlist(SL))
MAT           <- matrix(0, length(Rw_names), length(SL) )
rownames(MAT) <- Rw_names

# Calculate the distance matrix 
for(xx in 1:length(Rw_names)){
  for(yy in 1:length(SL)){
    m1 <- Rw_names[xx]
    m2 <- SL[[yy]][[1]]
    MAT[xx,yy] <- sum(m1==m2) 
  }
}  

colnames(MAT) <- 1:length(SL)
MAT           <- t(MAT)

return(MAT)

########################################
# Calculate tfidf to normalize the data 
########################################
MAT <- calcTFIDF2(MAT)
#######################################

# Dimension Reduction ###################################
if(reduce.dim.pca){
 ####PCA#################################################
 #print("Reduce Dimension PCA")
 log.ir <- t(MAT)
 ir.pca <- prcomp(log.ir)#, center=TRUE, scale.=TRUE)
 #return(summary(ir.pca))
 #plot(ir.pca, type="l", col="orange")
 MAT <- ir.pca$rotation[,1:reduce.dim.pca]#works
 ########################################################
}

#LSA/SVD  Dimension reduction ##########################
if(reduce.dim.lsa){
 #print("Reduce Dimension LSA")
 X <- t(MAT)
 s <- svd(X) # Singulärwertzerlegung
 D <- diag(s$d) # Singulärwerte
 n_s  <- length(s$d)
 Xnew  <- diag(s$d[1:reduce.dim.lsa])%*%(t(s$v)[1:reduce.dim.lsa,]) # weil bei der transponierten Zeilen ...  
 MAT   <- t(Xnew)
}
########################################################

#HCLUS##################################################
#d  <- dist(scale(MAT), method="cosine")
d  <- dist(scale(MAT), method="cosine")

#return(d)
#the distance measure to be used. This must be one of "euclidean", 
#"maximum", "manhattan", "canberra", "binary" or "minkowski"
hc <- hclust(d, method="average") 
#return(hc)
#######################################################
# the agglomeration method to be used. This should be (an unambiguous abbreviation of) 
# one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), 
# "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#Create list of K clusters ############################## 
cl      <- cutree(hc, h=H)
#######################################################
#return(d)

#agnes #############
#hc  <- agnes(test)
#cl  <- cutree(hc, h=H)
#### -- same results as hclust 

##dbscan ##################################
#res <- dbscan(test, eps = H, minPts = 1)
#cl  <- res$cluster
###########################################

#kmeans################################################
#cl <- kmeans(d, H)
#cl <- cl$cluster
#######################################################

CLUSTER <- list()

n.cluster <- length(unique(cl))

for(xx in 1:n.cluster){
  CLUSTER[[xx]] <- diag[cl==xx]  
}
WORD.L <- NULL
WORD.L <- sapply(diag,nchar)
return(list(CLUSTER=CLUSTER,cl=cl,d=as.matrix(d)))#, WORD.L=WORD.L))
#####################################################

} # ENd of function 

AssignLabels <- function(CosMatrix, cluster, filename){
 
  # filename to get the already ICD10 coded diagnosis  
  data   <- read_xlsx(filename)
  head   <- names(data)[2]
  labels <- c(names(data)[2], data[,2][[1]])
  # 1 Clean labels (take just the first 3 identifier)
  labels <- sapply(strsplit(labels, split='.', fixed=TRUE), function(x) (x[1]))
  n.documents <- length(labels)
    
  labelsNEW   <- labels
  for(xx in 1:n.documents){
    
    if(!is.na(labels[xx])){next}# ignore alreay labeled docs
    vec    <- CosMatrix[xx,]
    # check the members of the cluster xx is in 
    cl     <- cluster[xx]
    memb   <- which(cl==cluster)
    memb[memb==xx]            <- NaN
    memb[is.na(labels[memb])] <- NaN
    #if(xx==176){
    #print(cl)
    #print(cl)
    #print(memb)
    #}
    id      <- which.min(vec[memb])
    minmemb <- memb[id]
    #print(minmemb)
    if(length(minmemb)==0){next}
    labelsNEW[xx] <- labels[minmemb]
    #print(labels[minmemb])
    #stop("Test")
  }
names(labelsNEW) <- 1:length(labels)  

# As a last step iterate through the documents and append .x to those 
# documents which are syntoms instead of diagnosis
snps <- c(names(data)[4], data[,4][[1]])
#return(snps)
for(xx in 1:n.documents){
 
  stp <- grep("st p",snps[xx])
  zn  <- grep("z n",snps[xx])
  
  if((length(stp)!=0) | (length(zn)!=0)){labelsNEW[xx] <- paste(labelsNEW[xx],".x", sep="")}
}
###################################################################

return(labelsNEW)  
}

ModCluster <- function(cluster,newICD10){
  # This function takes the ICD10 code inormation 
  # to post-process/modify the clustering 
  #ower cluster assignments 
  #ower ICD10 assignments 
  newcluster    <- cluster
  newnewcluster <- cluster  
  ICD           <- unique(newICD10)
  #print(ICD)
  for(xx in 1:length(ICD)){
    if(is.na(ICD[xx])){next}
    if(ICD[xx]=="NA.x"){
    idd        <-  newICD10=="NA.x"
    idd[is.na(idd)] <- FALSE
    newcluster[idd] <- newcluster[idd] + 10000;next
    }
    newcluster[newICD10==ICD[xx]] <- xx + 1000 
  }
  
  idd <- unique(newcluster)
  for(xx in 1:length(idd)){
    newnewcluster[newcluster==idd[xx]] <- xx
  }
  
return(newnewcluster)
}

### SUBFUNCTIONS ###################################################################

LinePlot <- function(df){
  
  ggplot(data=df, aes(x=(100/75)*id, y=Fmax, group=1)) +
         geom_line() +
         geom_point() +
         ylim(0,1) +
         geom_hline(yintercept = 0.7) +
    labs(title="",y="F1-measure", x = "Dimension (%)")
  
}

BarPlot <- function(){
  df <- data.frame(measure=c("Precision", "Recall", "F1-measure"),value=c(0.8165, 0.6223, 0.7063))
  p  <- ggplot(data=df, aes(x=measure, y=value)) + 
        geom_bar(stat="identity") +
        ylim(0,1) +
        geom_text(aes(label=value), vjust=-0.3, size=4) +
        labs(title="",x="", y = "") #+
        #theme_minimal()
  plot(p)
}

### write Output file ###
WriteOUT <- function(filename, cluster=FALSE, newICD10=FALSE, mod.cluster=FALSE){
  
  data  <- read_xlsx(filename)[,]
  h1    <- names(data)[1]
  h2    <- names(data)[2]
  h3    <- names(data)[3]
  h4    <- names(data)[4]
  
  col1  <- c(h1,unlist(data[,1]))
  col2  <- c(h2,unlist(data[,2]))
  col3  <- c(h3,unlist(data[,3]))
  col4  <- c(h4,unlist(data[,4]))
  
  
  DAT <- cbind(col1,col2,col3,col4,cluster,newICD10, mod.cluster)
  rownames(DAT) <- 1:dim(DAT)[1]
  
  DAT<-DAT[order(newICD10),]
  colnames(DAT) <- c("cluster-per-hand","ICD10","document","norm","aut.cluster","ass. ICD10","mod.cluster")
  return(DAT)
  
} 

### Calculate Elbow Method ###################################
calcElbow <- function(d, K){

  # MMM is the distance matrix 
  # Note, K is now a range of K values
  # Iterate of all K and calculate sum SSE
  SSE   <- numeric(length(K))
  count <- 1
  for(xx in K){
    hc  <- hclust(d, method="average")#method="ward")
    cl  <- cutree(hc, xx)
    # calculate SSE
    D   <- as.matrix(d)
    CL  <- as.matrix(cl)
    CL  <- cbind(as.numeric(rownames(CL)),CL)
    sse <- numeric(xx)
    for(yy in (1:xx)){
      ids     <- CL[,2]==yy
      CLsub   <- CL[ids,,drop=FALSE]
      el      <- apply(CLsub,1,function(x){return(D[x[1],x[2]])})
      sse[yy] <- (el-mean(el))^2
    }
  SSE[count] <- sum(sse)
  count      <- count + 1
  }
return(SSE)  
}
### Calculate tfidf ##############################
calcTFIDF <- function(MAT){
  #input is the unnormalized count matrix 
  #rows: douments 
  #cols: words
  #tf 
  n.words <- apply(MAT,1,sum)
  tf      <- MAT/n.words #TF(t) = (Number of times term t appears in a document) / (Total number of terms in the document).
  idf     <- apply(MAT,2,function(x){return(sum(x!=0))})
  idf     <- log(dim(MAT)[1]/idf) #IDF(t) = log_e(Total number of documents / Number of documents with term t in it). 
  idf     <- matrix(rep(idf,dim(MAT)[1]),dim(MAT)[1], dim(MAT)[2], byrow=TRUE)
  return(tf*idf)
}

calcTFIDF2 <- function(MAT){
  d     <- MAT
  tf    <- d
  idf   <- log(nrow(d)/colSums(d))
  #print(idf)
  tfidf <- d
  
  for(word in names(idf)){
    tfidf[,word] <- tf[,word] * idf[word]
  } 
  return(tfidf)
}

calcOkapi <- function(MAT, k1, b) {
  
  tf  <- MAT
  idf <- log(nrow(MAT)/colSUms(MAT))
  k1 <- 1.5
  b <- 0.75
  OkapiBM <- idf*((tf*(k1+1))/(tf+(k1*(1-b+b*(length(rowSums(tf))/mean(rowSums(tf)))))))
  
  return(OkapiBM)
  
}

EvalClust <- function(L,TL){
################################  
# L  is vector of cluster as returned from hclust 
# TL is the true vector of cluster as returned from hclust 
# for simplifications just numeric IDs within the list   
################################
docs   <- unlist(L)  
n.docs <- length(docs)  

EVAL   <<- matrix(0,2,2) 
rownames(EVAL) <- c("Same class","Different class")
colnames(EVAL) <- c("Same cluster","Different cluster")

# Calculate TP, FP, TN, TP
pairs <- combn(n.docs,2)

 apply(pairs,2, function(x){
 d1 <- x[1]
 d2 <- x[2]
 # Check L 
 EVAL[1,1] <<- EVAL[1,1] + ((TL[d1] == TL[d2]) & (L[d1] == L[d2]))
 EVAL[2,2] <<- EVAL[2,2] + ((TL[d1] != TL[d2]) & (L[d1] != L[d2]))
 EVAL[1,2] <<- EVAL[1,2] + ((TL[d1] == TL[d2]) & (L[d1] != L[d2]))
 EVAL[2,1] <<- EVAL[2,1] + ((TL[d1] != TL[d2]) & (L[d1] == L[d2]))   
})
  
TP <- EVAL[1,1]
TN <- EVAL[2,2]
FP <- EVAL[2,1]
FN <- EVAL[1,2]

cat("TP",TP,"\n")
cat("FP",FP,"\n")
cat("FN",FN,"\n")
cat("TN",TN,"\n")

P  <- TP/(TP+FP)
R  <- TP/(TP+FN)
F1 <- (2*P*R)/(P+R) 

RandI <- (TN+TP)/(TP+TN+FP+FN)

return(list(P=P,R=R,F1=F1, RandI=RandI)) 
}

calcHeatMap <- function(L,TL){
  
  ### Heatmap ##############################
  docs <- names(L)
  
  heat <<- matrix(0, length(docs), length(docs))
  rownames(heat) <<- docs
  colnames(heat) <<- docs
  
  pairs <- combn(length(docs),2)
  
  apply(pairs,2, function(x){
    d1 <- x[1]
    d2 <- x[2]
    # Check L 
    if(((TL[d1] == TL[d2]) & (L[d1] == L[d2]))){heat[d1,d2]<<-"TP"}
    if(((TL[d1] != TL[d2]) & (L[d1] != L[d2]))){heat[d1,d2]<<-"TN"}
    if(((TL[d1] == TL[d2]) & (L[d1] != L[d2]))){heat[d1,d2]<<-"FN"}
    if(((TL[d1] != TL[d2]) & (L[d1] == L[d2]))){heat[d1,d2]<<-"FP"}   
  })
  ##########################################
  dat                  <- data.frame(doc=as.numeric(docs), heat)
  colnames(dat)[-1]    <- as.numeric(docs)
  colnames(dat)[1]     <- "doc"

  return(dat)
  
}

plotHeatMap <- function(dat){
  library(ggplot2); library(reshape2)
  dat3 <- melt(dat, id.var = 'doc')
  ggplot(dat3, aes(variable, doc)) + 
    geom_tile(aes(fill = value), colour = "white") + 
    scale_fill_manual(values=c("red", "blue", "black","green","pink"))+
    scale_y_continuous(breaks = 1:dim(data)[1]) +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank() ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(text = element_text(size=5))
}

plotHeatMap_DIST <- function(dat){
  library(ggplot2); library(reshape2)
  mat.melted <- melt(dat)
  ggplot(data = mat.melted, aes(x=X1, y=X2, fill=value)) +
    scale_y_continuous(breaks = 1:dim(data)[1]) +
    scale_x_continuous(breaks = 1:dim(data)[1]) +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank() ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(text = element_text(size=5)) +
    geom_tile()
}

## Function to find cutoff 
calcF <- function(f, filename, reduce=FALSE, rws=FALSE){

cutoffs <- seq(0,1,by=0.01) #hclust
#cutoffs <- seq(10,50,by=1) #kmeans  

# f are the real cluster 
F1 <- matrix(0,length(reduce),length(cutoffs))
P  <- matrix(0,length(reduce),length(cutoffs))
R  <- matrix(0,length(reduce),length(cutoffs))

for(yy in 1:length(reduce)){
  cat(yy, " of ", length(reduce),"\n")
 for(xx in 1:length(cutoffs)){
 #cat("Run with Cluster: ",xx,"\n")  
 test <- ICD10_cluster(filename, cutoffs[xx], reduce.dim.lsa = reduce[yy], rows=rws)  
 res  <- EvalClust(test[[2]],f)
 F1[yy,xx] <- res$F1
 P[yy,xx]  <- res$P
 R[yy,xx]  <- res$R
 }  
}
colnames(F1) <- cutoffs
rownames(F1) <- reduce
return(F1)#list(F1,P,R))
}

### legend(-1, 1.9, c("F", "P", "R"), col = c("orange","black","blue"), text.col = "green4", lty = c(2, 2, 2))
##  pch = c(19, 19, 19)),merge = TRUE, bg = "gray90")
###
###