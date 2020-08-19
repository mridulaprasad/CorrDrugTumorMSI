#'Spatial correlogram for Moran's I
#'
#'@description creates spatial correlogram of input variable at different lag distances 
#'@param  w a two-dimensional matrix with xy-coordinates of input variable
#'@param x a numeric vector
#'@param distr maximum lag order
#'@param plot logical operator to decide make spatial correlogram plot
#'@return two-dimension matrix of lag distance and corresponding autocorrelation value 
#'@export 
#'
#'

SpatialCorrelogram <- function(w,x,distr=2,Fig=T){
  w <- as.matrix(dist(w))
  aa <- ceiling(max(w)/distr)
  dists <- seq(0,aa*distr,distr)
  dists <- dists[1:31]
  cors <- NULL
  for(i in 1:30){
    w1 <- ifelse(w > dists[i] & w <= dists[i+1], 1, 0) 
    w2 <- w1
    for(j in 1:dim(w1)[1]){
      nu <- sum(w1[j,])
      if(is.na(nu) != TRUE & nu>0){
        w2[j,] <- w1[j,]/nu
      }  
    }
    lag <- w2 %*% x
    cors <- c(cors,cor(x,lag))
    cors[is.na(cors)]<-0
  }
 if(Fig=='T') plot(dists[c(1:30)],cors,type='o',xlab='lagdist',ylab='Autocorrelation value')
  return(cbind(dists[1:30],cors))
} 