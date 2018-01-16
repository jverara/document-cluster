library(readxl)
NormalizationResults <- read_excel("~/Documentos/Markus/Text clustering/NormalizationResults_002_cut.xlsx", col_names = FALSE)
NormalizationResults_FullDocument <- read_excel("~/Documentos/Markus/Text clustering/NormalizationResults_003_cut.xlsx", col_names = FALSE)

# OkapiBM25
# strsplit(pangram, " ") --> To split a string into vector
# Create a vector with all the different snippets
snippet_dictionary <- c()
split_snippets_list <- c()
for (i in 1:nrow(NormalizationResults)) {
  text_snippets_split <- unlist(strsplit(as.character(NormalizationResults[i,2]), " "))
  split_snippets_list[[i]] <- text_snippets_split
  for (j in 1:length(text_snippets_split)) {
    if (!(text_snippets_split[j] %in% snippet_dictionary)) {
      snippet_dictionary <- c(snippet_dictionary, text_snippets_split[j])
    }
  }
}

idf <- c()
for  (i in 1:length(snippet_dictionary)) {
  docs_with_term = 0
  for (j in 1:length(split_snippets_list)) {
    if (snippet_dictionary[i] %in% split_snippets_list[[j]]) {
      docs_with_term <- docs_with_term + 1
    }
  }
  idf <- c(idf, (log(nrow(NormalizationResults))/(docs_with_term)))
}
names(idf) = snippet_dictionary

# Create matrix of vectors
total_lengths <-c()
for (i in 1:length(split_snippets_list)) {
  total_lengths <- c(total_lengths, length(split_snippets_list[[i]]))
}
avgdl <- mean(total_lengths) #avgdl = average documents' length

mapping_vector_matrix <- matrix(NA, nrow = length(snippet_dictionary), ncol = 0)
rownames(mapping_vector_matrix) = snippet_dictionary
mapping_vector_matrix_tf.idf <- matrix(NA, nrow = length(snippet_dictionary), ncol = 0)
rownames(mapping_vector_matrix_tf.idf) = snippet_dictionary
k1 = 1.5 # k1 in [1.2,2.0]
b = 0.75
for (i in 1:length(split_snippets_list)) {
  appending_vector <- integer(length(snippet_dictionary))
  for (j in 1:length(split_snippets_list[[i]])) {
    match_position <- match(split_snippets_list[[i]][j], snippet_dictionary)
    appending_vector[match_position] <- appending_vector[match_position] + 1
  }
  appending_vector_tf.idf <- integer(length(idf))
  for (k in 1:length(idf)) {
    #appending_vector_tf.idf[k] <- appending_vector[k]*idf[k]
    appending_vector_tf.idf[k] <- idf[k]*((appending_vector[k]*(k1+1))/(appending_vector[k]+(k1*(1-b+b*(length(split_snippets_list[[i]])/avgdl)))))
  }
  mapping_vector_matrix <- cbind(mapping_vector_matrix, appending_vector)
  mapping_vector_matrix_tf.idf <- cbind(mapping_vector_matrix_tf.idf, appending_vector_tf.idf)
  colnames(mapping_vector_matrix)[ncol(mapping_vector_matrix)] <- paste("Doc", ncol(mapping_vector_matrix), sep = " ")
  colnames(mapping_vector_matrix_tf.idf)[ncol(mapping_vector_matrix_tf.idf)] <- paste("Doc", ncol(mapping_vector_matrix_tf.idf), sep = " ")
}

mapping_vector_matrix <- scale(mapping_vector_matrix)
mapping_vector_matrix_tf.idf <- scale(mapping_vector_matrix_tf.idf)
## Clustering process ##

#distances <- dist(t(mapping_vector_matrix), method = "euclidian")
#distances_tf.idf <- dist(t(mapping_vector_matrix_tf.idf), method = "euclidian")
library(proxy)
distances_tf.idf_cos.sim <- dist(t(mapping_vector_matrix_tf.idf), method = "cosine")
#dendrogram<- hclust(t(distances), method="average")
#dendrogram_tf.idf<- hclust(t(distances_tf.idf), method="average")
dendrogram_tf.idf_cos.sim<- hclust(t(distances_tf.idf_cos.sim), method="average")
#plot(dendrogram, hang=-0.3, cex=0.4)
#plot(dendrogram_tf.idf, hang=-0.3, cex=0.4)
plot(dendrogram_tf.idf_cos.sim, hang = -1, cex=0.4)

#cor.pe <- cor(as.matrix(t(mapping_vector_matrix)), method=c("pearson")) 
#cor.pe.col <- cor(as.matrix(mapping_vector_matrix), method=c("pearson")) #cosine similarity
library(coop)
#cosine_cor <- cosine(as.matrix(mapping_vector_matrix))
cosine_cor_tf.idf <- cosine(as.matrix(mapping_vector_matrix_tf.idf))

# Cosine Similarity
#cos.sim <- function(ix) 
#{
#  A = mapping_vector_matrix[ix[1],]
#  B = mapping_vector_matrix[ix[2],]
#  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
#}   
#n <- nrow(mapping_vector_matrix) 
#cmb <- expand.grid(i=1:n, j=1:n) 
#cosine_cor <- matrix(apply(cmb,1,cos.sim),n,n)

#correlation heatmap
library(reshape2)
#melted_cor.pe.col <- melt(cor.pe.col)
#melted_cosine_cor <- melt(cosine_cor)
melted_cosine_cor_tf.idf <- melt(cosine_cor_tf.idf)

library(ggplot2)
#ggplot(data = melted_cor.pe.col, aes(x=Var1, y=Var2, fill=value)) + 
#  geom_tile()

#ggplot(data = melted_cosine_cor, aes(x=Var1, y=Var2, fill=value)) + 
#  geom_tile()

ggplot(data = melted_cosine_cor_tf.idf, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

#heatmap.2(cor.pe.col, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=cor.pe.col, notecol="black", trace='none', key=FALSE)

#dist.pe <- as.dist(1-cor.pe)
#dist.pe.col <- as.dist(1-cor.pe.col)
#dist.cos <- as.dist(1-cosine_cor)
#dist.cos_tf.idf <- as.dist(1-cosine_cor_tf.idf)
#library(proxy)
#dist.cos_tf.idf <- dist(cosine_cor_tf.idf, method = "cosine")

# Hierarchical clustering

#dendrogram.pearson <- hclust(dist.pe,method="average")
#dendrogram.pearson.col <- hclust(dist.pe.col,method="average")
#dendrogram.cos <- hclust(dist.cos, method = "average")
#dendrogram.cos_tf.idf <- hclust(dist.cos_tf.idf, method = "average")

#plot(dendrogram.pearson ,main="Documents grouped by Pearson",hang=-1, cex=0.3)
#plot(dendrogram.pearson.col,main="Documents grouped by Pearson",hang=-1, cex=0.3)
#plot(dendrogram.cos, main = "Documents grouped by cosine similarity", hang = -1, cex = 0.5)
#plot(dendrogram.cos_tf.idf, main = "Documents grouped by cosine similarity using tf-idf matrix", hang = -1, cex = 0.5)
#abline(h=0.5, col="red")

# Cutting trees at a fixed height
#cutree(dendrogram.cos, h = 0.5)
#cutree(dendrogram.cos_tf.idf, h = .95)
cutree(dendrogram_tf.idf_cos.sim, h=0.95)

#clusters <- cutree(dendrogram.cos, h = 0.5)
#clusters_tf.idf <- cutree(dendrogram.cos_tf.idf, h = 0.95)
clusters_tf.idf_cos.sim <- cutree(dendrogram_tf.idf_cos.sim, h=0.90)

# Substite numbered documents by document itself
clusters_variable <- clusters_tf.idf_cos.sim
clusters_doc <- character(max(clusters_variable))
for (i in 1:length(clusters_variable)) {
  cluster_number <- clusters_variable[i]
  if (clusters_doc[cluster_number] == "") {
    clusters_doc[cluster_number] <- as.character(NormalizationResults_FullDocument[i,2])
  }
  else {
    clusters_doc[cluster_number] <- paste(clusters_doc[cluster_number], as.character(NormalizationResults_FullDocument[i,2]), sep=", ")
  }
}
clusters_doc

# Represent every cluster separated
library(dendextend)
library(colorspace)
#dend_color <- color_branches(dendrogram.cos, h = 0.5)
#plot(dend_color)

#labels_dend <- labels(dend_color)
#groups <- cutree(dend_color, h=0.5, order_clusters_as_data = FALSE)
#dends <- list()
#for(i in 1:max(clusters)) {
#  labels_to_keep <- labels_dend[i != groups]
#  dends[[i]] <- prune(dend_color, labels_to_keep)
#}

#par(mfrow = c(2,3))
#for(i in 1:max(clusters)) { 
#  plot(dends[[i]], 
#       main = paste0("Tree number ", i))
#}
#par(mfrow = c(1,1))

# tf.idf
dend_color_tf.idf <- color_branches(dendrogram.cos_tf.idf, h = 0.87)
par(mfrow = c(1,1))
plot(dend_color_tf.idf)

labels_dend_tf.idf <- labels(dend_color_tf.idf)
groups_tf.idf <- cutree(dend_color_tf.idf, h=0.87, order_clusters_as_data = FALSE)
dends_tf.idf <- list()
for(i in 1:max(clusters_tf.idf)) {
  labels_to_keep_tf.idf <- labels_dend_tf.idf[i != groups_tf.idf]
  dends_tf.idf[[i]] <- prune(dend_color_tf.idf, labels_to_keep_tf.idf)
}

par(mfrow = c(2,3))
for(i in 1:max(clusters_tf.idf)) { 
  plot(dends_tf.idf[[i]], 
       main = paste0("Tree number ", i))
}
par(mfrow = c(1,1))

# Hierarchical clustering (alternative method)

library(pvclust)
hc <- pvclust(cosine_cor, method.hclust="average",
              method.dist="correlation")
plot(hc)
pvrect(hc, alpha=.85, pv = "bp")

hc_tf.idf <- pvclust(cosine_cor_tf.idf, method.hclust="average",
              method.dist="correlation")
plot(hc_tf.idf)
pvrect(hc_tf.idf, alpha=.60, pv = "bp")
pvpick(hc_tf.idf, alpha = .60, pv = "bp")

#heatmap(as.matrix(mapping_vector_matrix),
        #Rowv=as.dendrogram(dendrogram.pearson),
        #Colv=as.dendrogram(dendrogram.pearson.col), cexRow = 0.3, cexCol = 0.3)

# Iterative clustering (using kmeans)

k<-c(50)
dat.kmeans <- kmeans(t(mapping_vector_matrix),k,iter.max=1000)
dat.kmeans
clusters_doc_kmeans <- character(max(dat.kmeans$cluster))
for (i in 1:length(dat.kmeans$cluster)) {
  cluster_number <- dat.kmeans$cluster[i]
  if (clusters_doc_kmeans[cluster_number] == "") {
    clusters_doc_kmeans[cluster_number] <- as.character(NormalizationResults_FullDocument[i,2])
  }
  else {
    clusters_doc_kmeans[cluster_number] <- paste(clusters_doc_kmeans[cluster_number], as.character(NormalizationResults_FullDocument[i,2]), sep=", ")
  }
}
clusters_doc_kmeans

#par(mfrow=c(2,3))
#for(i in 1:k){
#  matplot(t(mapping_vector_matrix_tf.idf[dat.kmeans$cluster == i,]), type = "l",
#          main=paste("cluster:",i),ylab="expression log",xlab="Docs")
#}
#par(mfrow=c(1,1))

library(cluster)

#dat.gap.kmeans <- clusGap(mapping_vector_matrix, kmeans, K.max = 20)
#dat.gap.kmeans

library(fpc)
t_mapping_vector_matrix_tf.idf <- t(mapping_vector_matrix_tf.idf)
pamk(t_mapping_vector_matrix_tf.idf, krange = 2:60)
pam_data <- pamk(t_mapping_vector_matrix_tf.idf, krange = 2:60)

# K-Means Cluster Analysis
fit <- kmeans(t_mapping_vector_matrix_tf.idf, pam_data$nc) # 49 cluster solution

clusters_doc_kmeans <- character(max(fit$cluster))
for (i in 1:length(fit$cluster)) {
  cluster_number <- fit$cluster[i]
  if (clusters_doc_kmeans[cluster_number] == "") {
    clusters_doc_kmeans[cluster_number] <- as.character(NormalizationResults_FullDocument[i,2])
  }
  else {
    clusters_doc_kmeans[cluster_number] <- paste(clusters_doc_kmeans[cluster_number], as.character(NormalizationResults_FullDocument[i,2]), sep=", ")
  }
}
clusters_doc_kmeans

# get cluster means 
#aggregate(t_mapping_vector_matrix,by=list(fit$cluster),FUN=mean)
# append cluster assignment
#mapping_vector_matrix_cluster <- rbind(mapping_vector_matrix, cluster = t(fit$cluster))
#rownames(mapping_vector_matrix_cluster)[nrow(mapping_vector_matrix_cluster)] <- "Cluster"

# Number of Docs in every cluster
#H <- table(t(fit$cluster))
#barplot(H)

# Groups of docs
#data.kmedoid <- pam(mapping_vector_matrix, k=49)
#clusplot(data.kmedoid)

#DBScan

library(dbscan)
k=60
kNNdistplot(t(mapping_vector_matrix_tf.idf), k)
abline(h=45, col="red")
dat.dbscan <- dbscan(t(mapping_vector_matrix_tf.idf),eps=42, minPts = 5, borderPoints = FALSE)
dat.dbscan
