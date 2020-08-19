
#'Median filtering of image
#'
#'@description performs median filtering of the image with window size 3 x 3 
#'@param  sampleMat image for median filtering
#'@return returns median filtered image 
#'@export 


medianFilterR <- function(sampleMat)
{
  B = matrix(0, nrow= dim(sampleMat)[1],ncol=dim(sampleMat)[2])
  modifyA = matrix(0,nrow=(dim(sampleMat)[1]+2),ncol=(dim(sampleMat)[2]+2))
  
  for(x in 1:dim(sampleMat)[1])
  {
    for(y in 1:dim(sampleMat)[2])
    {
      modifyA[x+1,y+1] = sampleMat[x,y]
    }
  }
  for(i in 1:(dim(modifyA)[1]-2))
  {
    for(j in 1:(dim(modifyA)[2]-2))
    {
      window = rep(0,9)
      inc = 1
      for(x in 1:3)
      {
        for(y in 1:3)
        {
          window[inc] = modifyA[i+x-1,j+y-1]
          inc = inc +1
        }
      }
      med = sort(window)
      B[i,j] = med[5]
    }
  }
  return(B)
}