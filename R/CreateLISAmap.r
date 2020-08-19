#'Creates LISA map of input image
#'
#'@description Creates LISA image of data matrix 
#'@param  Ionimage two-dimensinal data matrix to create LISA map
#'@param distthrs lag distance value 
#'@param plot return LISA image in output
#'@param mask mask image of tissue to exclude backaground region in the calculation 
#'@return two-dimensional LISA map 
#'@export 

CreateLISAmap <- function(Ionimage,distthrs=5,plot='T',mask=NULL)
{
  ycord    <- dim(Ionimage)[1]; xcord <- dim(Ionimage)[2]
  Ionimage <- as.vector(Ionimage)
  xycoord  <-as.matrix(expand.grid(1:ycord,1:xcord))
  
  if(is.null(mask) == 'FALSE')
  {
    BkgPixelsIndex <- which(as.vector(mask) %in% 0)
    Ionimage <- Ionimage[-BkgPixelsIndex]
    xycoord <- xycoord[-BkgPixelsIndex,]
    
    neighbor <- matrix(0,nrow=dim(xycoord)[1],ncol=dim(xycoord)[1]); 
    nlistk1  <- spdep::dnearneigh(xycoord,0,distthrs);xy_weights <- spdep::nb2listw(nlistk1);
    wzx      <- spdep::lag.listw(xy_weights,scale(Ionimage))
    wzx[which(is.na(wzx))] <- 0
    LISAmap   <- scale(Ionimage)
    LISAmap[LISAmap >=0 & wzx >=0 ] <- 1; LISAmap[LISAmap >=0 & wzx <=0 ] <- 3; LISAmap[LISAmap <=0 & wzx >=0 ] <- 4; LISAmap[LISAmap <=0 & wzx <=0 ] <- 2
    FullLisaMap <- as.vector(mask)
    FullLisaMap[which(FullLisaMap %in% 1)] <- LISAmap
    dim(FullLisaMap) <- c(ycord,xcord)
    if(plot== 'T')image(FullLisaMap,axes=FALSE,col=c('black','blue','green','yellow','red'))
    
  }else{
  
  neighbor <- matrix(0,nrow=dim(xycoord)[1],ncol=dim(xycoord)[1]); 
  nlistk1  <- spdep::dnearneigh(xycoord,0,distthrs);xy_weights <- spdep::nb2listw(nlistk1);
  wzx      <- spdep::lag.listw(xy_weights,scale(Ionimage))
  wzx[which(is.na(wzx))] <- 0
  LISAmap   <- scale(Ionimage)
  LISAmap[LISAmap >=0 & wzx >=0 ] <- 1; LISAmap[LISAmap >=0 & wzx <=0 ] <- 3; LISAmap[LISAmap <=0 & wzx >=0 ] <- 4; LISAmap[LISAmap <=0 & wzx <=0 ] <- 2
  dim(LISAmap) <- c(ycord,xcord)
  if(plot== 'T')image(LISAmap,axes=FALSE,col=c('black','blue','green','yellow','red'))
  return(LISAmap)
}
  }
