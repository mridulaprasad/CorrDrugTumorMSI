#'Quantative analysis of segemented image and LISA map 
#'
#'@description The function calculates the  fraction of pixels from different clusters fall into different zones of the LISA map 
#'@param  ClusteredImg two-dimensinal clustered image
#'@param LISAmap two-dimensional LISA map
#'@return fraction of pixels from different clusters fall into different zones of the LISA map 
#'@export 




ClusteredLISAimgQA <- function( ClusteredImg,LISAmap,plot='T')
{
  img1 <- LISAmap; 
  img2 <- ClusteredImg
  tissuePixeids <- which(img1 %in% 1 |img1 %in% 2|img1 %in% 3|img1 %in% 4)
  Nucluster <- length(table(img2[tissuePixeids]))
  QAmatrix <- as.data.frame(matrix(0,nrow = 4,ncol=Nucluster))
  LISAzones <- c('HH','LL','HL','LH')
  rownames(QAmatrix) <- LISAzones
  colnames(QAmatrix) <- paste0('cluster',c(1:Nucluster))
  id1 <- which(img1 %in% 1); QAmatrix[1,] <- as.vector(table(img2[id1]))
  id1 <- which(img1 %in% 2); QAmatrix[2,] <- as.vector(table(img2[id1]))
  id1 <- which(img1 %in% 3); QAmatrix[3,] <- as.vector(table(img2[id1]))
  id1 <- which(img1 %in% 4); QAmatrix[4,] <- as.vector(table(img2[id1]))
  QAmatrix <- QAmatrix/tissuePixeids
  
  if(plot== 'T'){
    par(mfrow=c(2,2))
    for(i in 1:4){
    barplot(as.matrix(QAmatrix[i,]),ylim = c(0,max(QAmatrix)),main = LISAzones[i])}
  }
  
  return(QAmatrix)
}