#'Preprocessing of MSI data 
#'
#'@description This function does peaks selection, removal of less frequent peaks and high intensity pixels from the tissue edge in MSI data. 
#'
#'For the peaks selection, first the adaptive bins are created based on maximum intenstiy reference spectrum from single or multiple MSI datasets.
#'After bins creation, peaks selection performed from individual MSI dataset using  
#'local maxima search method. Threshold value to exclude the noisy peaks in the peaks selection step is automatically selected from ion-intensity region of the individual mass spectrum.
#'
#'@importFrom ifultools isVectorAtomic 
#'@importFrom  msProcess msSet
#'@importFrom  msProcess msPeakSimple
#'@importFrom msProcess msDenoiseWavelet
#'@param folderpath path of MSI data folder. At the moment our function only supports the MSI data from in Analyze 7.5 format.
#'@param removePeaks remove the bins from MSI data with lower coverage inside tissue below defined threshold 
#'@param edgecorrection remove high intensity pixels from the edge of the tissue. Logical operator with default value is T. If user does not want to perfom this action on their data then set to 'F'.
#'@return two-dimensional matrix where total rows equal to total mass spectra in MSI data and total columns equal to bins created. 
#'@export 
#'@examples 
#'folderpath <- 'MSIdata//AVA-PTX'
#'PreprocessedMSIData <- PreprocessingMSIData(folderpath) # Does all default tasks:peaks selection, remove peaks with less coverage area and does edge correction
#'PreprocessedMSIData <- PreprocessingMSIData(folderpath, removePeaks=0.0) # Does peaks selection and edge correction. 
#'PreprocessedMSIData <- PreprocessingMSIData(folderpath,removePeaks=0.0,edgecorrection=F ) # Only peaks selection 
#'IonInensityMatrix   <- PreprocessedMSIData[[1]][[1]]
#'ycord               <- PreprocessedMSIData[[1]][[3]]; 
#'xcord               <- PreprocessedMSIData[[1]][[2]]; 
#'SingleionImage      <- IonInensityMatrix[,41]
#'dim(SingleionImage) <- c(ycord,xcord)
#'image(SingleionImage)
#'

PreprocessingMSIData <- function(folderpath,removePeaks=0.2,edgecorrection=T)
{
  
  ## Create reference maximum intensity spectra 
  
  dirlist = list.dirs(folderpath)  # read data
  analyfie1 = suppressMessages(MALDIquantForeign::importAnalyze(dirlist[1]))
  maxIntensity = matrix(0,nrow=length(dirlist),ncol=length(analyfie1[[100]]@mass))
  maxMass = matrix(0,nrow=length(dirlist),ncol=length(analyfie1[[100]]@mass))
  MSIdatafiles <- list()
  count =0 
  for(i in 1:length(dirlist))
  {
    print(i)
    MSIdatafiles[[i]] <- suppressMessages(MALDIquantForeign::importAnalyze(dirlist[i]))
    analyfie1         <- MSIdatafiles[[i]]
    IntensityMat      <- matrix(0,nrow=length(analyfie1),ncol=length(analyfie1[[100]]@mass))
    for(j in 1:length(analyfie1))
    {
      if(max(analyfie1[[j]]@intensity) < 100)
      {IntensityMat[j,] <- 0; count = count+1}
      else
      {
        IntensityMat[j,] <- analyfie1[[j]]@intensity
      }
    }
    
    maxIntensity[i,] <- apply(IntensityMat[,c(1:dim(maxIntensity)[2])],2,max) ; maxMass[i,] = analyfie1[[1100]]@mass[1:dim(maxIntensity)[2]]
  }
  
  maxMeanInt     <- apply(maxIntensity,2,max); maxMeanMass = apply(maxMass,2,mean)
  
  
  ### Bin boundary creation 
  
  
  DenoiseSpectra <- msProcess::msDenoiseWavelet(maxMeanInt, wavelet="s8", shrink.fun="hard", thresh.fun="universal",thresh.scale=2, xform="modwt")
  z              <- msProcess::msSet(as.matrix(DenoiseSpectra),mz = as.vector(maxMeanMass))
  z              <- msProcess::msPeakSimple(as.vector(z$mz),as.vector(z$intensity),span = 5,snr.thresh =median(sort(unique(maxMeanInt[2:200]))))
  msPeaks        <- cbind(z$mass.loc,z$intensity)
  
  
  for(i in 1:(dim(msPeaks)[1]-1))
  {
    ## Merge peaks with difference less than 0.05 
    if((msPeaks[i+1,1] - msPeaks[i,1]) < 0.05)
    {
      ## If difference between two identified peaks less than 0.05, then select one with maximum intensity value 
      
      id  <- which(msPeaks[,2] %in% max(msPeaks[i+1,2],msPeaks[i,2]) )
      if(length(id) !=1)
      {
        temp <- which(c(i,i+1) %in% id )
        if(temp == 1){id <- i}
        else{id <- i+1}
      }
      msPeaks[i,1]  <- msPeaks[id,1];msPeaks[i,2] = msPeaks[id,2]
      msPeaks[i+1,1]<- msPeaks[id,1];msPeaks[i+1,2] = msPeaks[id,2] 
    }
  }
  
  Peaklist    <- unique(msPeaks[,1]);extrabin <- c()
  
  binvalue    <- rep(0,(length(Peaklist)[1])+1)
  binvalue[1] <-  Peaklist[1] - 0.1
  c=2
  for(i in 1:length(Peaklist))
  { 
    if(i == length(Peaklist))
    {binvalue[c] <- Peaklist[i] + 0.1 }
    else
    {
      if((Peaklist[i+1] - Peaklist[i]) >0.7)
      {  extrabin <- c(extrabin,seq(Peaklist[i]+0.2,Peaklist[i+1]-0.2,by=0.2))}
      else{
        binvalue[c] <- sum(Peaklist[i],Peaklist[i+1])/2 ;c = c+1}
      
    }
    
  }
  
  
  ## Extra extra bins after last peak identified till the end of the spectrum
  binvalue  <- c(extrabin,binvalue)
  binvalue  <- sort(binvalue)
  binadd    <- seq(binvalue[length(binvalue)],502.5,by =0.5)
  binlength <- c(binvalue,binadd)
  binlength <- binlength[-which(duplicated(binlength))]
  
  cuts          <- cut(analyfie1[[1234]]@mass,binlength)
  duration.freq <- table(cuts)
  duration.freq <- cbind(duration.freq)
  cutdur        <- row.names(duration.freq)
  binsplit      <- function(x){as.numeric(unlist(strsplit(gsub("\\(|\\]", "", x),','))[1])}
  bin_even      <- sapply(cutdur,binsplit)
  
  binsplit     <- function(x){as.numeric(unlist(strsplit(gsub("\\(|\\]", "", x),','))[2])}
  bin_odd      <- sapply(cutdur,binsplit)
  
  ################## Do peak picking inside bins created 
  
  binned  <- cbind(bin_even,bin_odd)
  binname <- apply(binned,1,mean)
  
  ProcessedFiles <- list()
  
  for(s in 1:length(MSIdatafiles))
  {
    analyfie1 = MSIdatafiles[[s]]
    IntenMatrix = matrix(0,nrow=length(analyfie1),ncol = length(bin_odd)) 
    colnames(IntenMatrix) = binname
    
    for(i in 1:length(analyfie1))
    {
      
      if(max(analyfie1[[i]]@intensity) < 100)
      {
        IntenMatrix[i,] =0}
      
      else
      {
       z <- msProcess::msSet(as.matrix(analyfie1[[i]]@intensity),mz = as.vector(analyfie1[[i]]@mass))
        z <- msProcess::msPeakSimple(as.vector(z$mz),as.vector(z$intensity),span = 5,snr.thresh =median(sort(unique(z$intensity[2:200]))))
        mspeaks <- cbind(z$mass.loc,z$mass.right,z$mass.left)
        for(j in 1:dim(mspeaks)[1])
        {
          if(mspeaks[j,1]>range(maxMeanMass)[2]-1)
          {
            break
          }
          
          else
          {
            id <- which(mspeaks[j,1] > bin_even & mspeaks[j,1] < bin_odd)
            ## id length zero if peaks appear on bin boundary, then substract some value from peak to find m/z bin
            if(length(id) ==0)
            {id = which((mspeaks[j]-0.01) > bin_even & (mspeaks[j]-0.01) < bin_odd)}
            id = id[[1]]
            if(IntenMatrix[i,id] !=0){
              IntenMatrix[i,id] = max(IntenMatrix[i,id],z$intensity[j])}   #### If already peak present then select the one with high intensity value 
            else
            {IntenMatrix[i,id] = z$intensity[j]}
          }
        }
      }
    }
    IntenMatrix[is.na(IntenMatrix)] = 0
    
    x = analyfie1[[length(analyfie1)]]@metaData$imaging$pos[[1]]
    y = analyfie1[[length(analyfie1)]]@metaData$imaging$pos[[2]]
    
    IonimageFormasking <- IntenMatrix[,which.max(colSums(IntenMatrix))]
    dim(IonimageFormasking) <- c(y,x)
    mask <- medianFilterR(IonimageFormasking)
    
    ##### Masking based on first histogram block cut-off
    
    minimums <- function(x) which(x - shift(x, 1) < 0  & x - shift(x, 1, type='lead') < 0)
    
    d = density(mask)
    maskvalue = d$x[minimums(d$y)[1]]
    mask = mask > maskvalue
   
    if(edgecorrection == 'T'){
      for(r in 1:dim(IonimageFormasking)[1]) {
        id = which(mask[r,] %in% 1);     mask[r,id[1]] = 0; mask[r,id[length(id)]] =0}
      for(r in 1:dim(IonimageFormasking)[2]){
        id = which(mask[,r] %in% 1);     mask[id[1],r] = 0; mask[id[length(id)],r] =0}
      for(r in 1:dim(IonimageFormasking)[1]) {
        id = which(mask[r,] %in% 1);     mask[r,id[1]] = 0; mask[r,id[length(id)]] =0}
      for(r in 1:dim(IonimageFormasking)[2]){
        id = which(mask[,r] %in% 1);     mask[id[1],r] = 0; mask[id[length(id)],r] =0}
      for(r in 1:dim(IonimageFormasking)[1]) {
        id = which(mask[r,] %in% 1);     mask[r,id[1]] = 0; mask[r,id[length(id)]] =0}
      for(r in 1:dim(IonimageFormasking)[2]){
        id = which(mask[,r] %in% 1);     mask[id[1],r] = 0; mask[id[length(id)],r] =0}
      idImg = which(as.vector(mask) %in% 0) 
      IntenMatrix[idImg,] = 0  ### In intensity matrix make pixel value outside tissue region zero
      
    }
    
    idImg = which(as.vector(mask) %in% 0) 
    countMat = IntenMatrix[-idImg,]
    countMat[countMat != 0] = 1
    Peakfraction = ceiling(dim(countMat)[1] * removePeaks) 
    PeakId = which(colSums(countMat) < Peakfraction)
    IntenMatrix = IntenMatrix[,-c(PeakId)]
   
    ProcessedFiles[[s]] = list(IntenMatrix,x,y) 
  }
  return(ProcessedFiles)
}