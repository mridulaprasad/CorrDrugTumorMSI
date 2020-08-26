#'Find optimal number of clusters
#'
#'@description function returns optimal number of clusters based on calinhara internal cluster validity index method for k-means correlation distance based clustering method 
#'@param  Inputdata data to calculate number of clusters 
#'@return number of clusters
#'@export 
#'
#'

optimalClust <- function(Inputdata)
{
  wss  <- 0
  mod1 <- sapply(2:10, function(x) amap::Kmeans(Inputdata, x,method="correlation",iter.max = 150)$cluster); 
  wss  <- which.max(apply(mod1,2,function(x) fpc::calinhara(Inputdata,x))); 
  wss  <- wss +1
  return(wss)
}

#'Bisect k-means clustering of the data 
#'
#'@description Perform bisect k-means clustering on a data-matrix. 
#'@details The data given by x is clustered by the k-means algorimthm based correlation distance. The 
#'number of clusters automatically decides according calinhara internal cluster validity index.
#'The split of cluster into another cluster stop if the size of become less than cluster.size vlaue which the
#'fraction of data points into compare compare to total size of input data.  
#'@param  x A numeric matrix of data. 
#'@param cluster.size The fraction of pixels in individual cluster. Particular cluster will not split further if 
#'contains this fracion of data-points with respect total data-points. 
#'@return A vector of integers indicating the cluster to which each point is allocated.
#'@export 
#'@seealso optimalClust
#'@examples 
#'x <- rbind(matrix(rnorm(500,sd=0.3),ncol=2),
#'           matrix(rnorm(500,mean=1,sd=0.3),ncol=2))
#'cl <- BisectKmeansClustering(x)
#'plot(x,col=cl)               
#'


BisectKmeansClustering <- function(Inputdata,cluster.size=0.5)
{
mainCluster <-c();count <-1; structClusternames <- rep(0,dim(Inputdata)[1])
for(m in 1:10)
{
  if(count ==1)
  {
    wss <- optimalClust(Inputdata)
    cl <- amap::Kmeans(Inputdata, wss,method="correlation",iter.max = 150)
    mainCluster <- cl$cluster
  }
  ClusterName <- as.numeric(names(table(mainCluster)))
  ClusterCount <- as.vector(table(mainCluster))/length(mainCluster)
  if(all(ClusterCount <=cluster.size))
  {break}
  else
  {
    id           <- which(ClusterCount >cluster.size)
    rows         <- which(mainCluster %in% ClusterName[id[1]])
    SplitCluster <- Inputdata[rows,]
    wss          <- optimalClust(SplitCluster)
    cl           <- amap::Kmeans(SplitCluster, wss,method="correlation",iter.max = 150)
    cl$cluster   <- cl$cluster+max(mainCluster)
    mainCluster[rows] <- cl$cluster
  }
  count <- count+1
}
  clusternames <- as.numeric(names(table(mainCluster)))
  clusternames <- sort(clusternames)
  for(i in 1:length(clusternames))
   { id <- which(mainCluster %in% clusternames[i])
   structClusternames[id] <- i}

  return(structClusternames)
}
  