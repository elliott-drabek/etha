#! /usr/local/bin/Rscript --vanilla

###############################################################################

args <- commandArgs(TRUE)

pdf(args[2])

distances <- read.table(args[1], header=TRUE, sep='\t', comment.char='', quote='', row.names=1)
distances <- as.dist(as.matrix(distances))

output_file_name <- args[3]

###############################################################################

min_dist <- 0.02

clustering <- hclust(distances)
if (length(distances) > 1) {
  plot(clustering)
  abline(h=min_dist, col='red')
}

clusters <- cutree(clustering, h=min_dist)
#cat('clusters:\n')
#print(clusters)

blue <- function(distances, clusters, i) {
  n <- length(clusters)

  index <- n*(i-1) - i*(i-1)/2 + j-i
}

find_centroid <- function(distances, clusters, cluster) {
  members <- names(clusters[clusters == cluster])
  id <- if (length(members) == 1)
    members
  else
    names(which.min(rowSums(as.matrix(distances)[members, members])))
}

cluster.set <- unique(clusters)
output <- data.frame(
  centroid=sapply(cluster.set, function(cluster) find_centroid(distances, clusters, cluster)),
  members=sapply(cluster.set, function(cluster) paste(names(clusters[clusters == cluster]), collapse=' '))
)

write.table(output, file=output_file_name, quote=F, sep='\t', row.names=F, col.names=F)

###############################################################################

#mds.fit <- cmdscale(distances, eig=TRUE, k=6)

#barplot(mds.fit$eig / sum(mds.fit$eig), ylab='Portion of variance explained', xlab='Component number', main='PCoA component loading (for following figures)')

#for (second.pc in 2:6) {
  #plot(mds.fit$points[,1], mds.fit$points[,second.pc], main='PCoA', xlab='Component 1', ylab=sprintf('Component %d', second.pc), col='red')
#}
