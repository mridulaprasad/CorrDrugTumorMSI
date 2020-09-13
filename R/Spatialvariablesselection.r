#'Variables selection 
#'
#'@description This function does variables selection using either spatial or non-spatial regression method. 
#'At the moment, function does variables selection for two class problem only where class of interest label as 2.
#'@param  Inputdata ion intensity matrix
#'@param classids class labels for variables selection
#'@param method method used for variables selection, 'OLS', 'SL', 'SE'.
#'@param lagdist the range of autocorrelation in spatial weight matrix for spatial methods
#'@param xcord,ycord the xy co-ordinates of the MSI ion image
#'@param maskimage exclude background pixels during variables selection
#'@return regression coefficients and p-value for all variables 
#'@export 
#'
#'

VariablesSelection <- function(Inputdata,classids, method='OLS',lagdist=5,xycoord=0)
{
  SpMod_estm<- c();SpMod_pval<- c()
## OLS Model 
if(method=='OLS'){
for(i in 1:dim(Inputdata)[2]){
  results       <- lm(Inputdata[,i] ~ 0+as.factor(classids))
  mt            <- summary(results); 
  SpMod_estm[i] <- mt$coefficients[2];
  SpMod_pval[i] <- mt$coefficients[8];   
} 
metlist <- list(SpMod_estm,SpMod_pval )}
if(method=='SE'){
## SE Model 
  # derive spatial weight matrix 
  xycoord  <- as.matrix(xycoord)
  neighbor <- matrix(0,nrow=max(xycoord[,1]),ncol=max(xycoord[,2])); 
  nlistk1  <- spdep::dnearneigh(xycoord,0,lagdist);xy_weights <- spdep::nb2listw(nlistk1);
  for(i in 1:dim(Inputdata)[2]){
  results      <- spatialreg::errorsarlm(Inputdata[,i] ~ as.factor(classids),listw = xy_weights)
  mt           <- summary(results); 
  SpMod_estm[i]<-  mt$Coef[2];
  SpMod_pval[i] <- mt$Coef[8];  
}
metlist <-list(SpMod_estm,SpMod_pval)}

## SL Model 
  if(method=='SL'){  
    xycoord  <- as.matrix(xycoord)
    neighbor <- matrix(0,nrow=max(xycoord[,1]),ncol=max(xycoord[,2])); 
    nlistk1  <- spdep::dnearneigh(xycoord,0,lagdist);xy_weights <- spdep::nb2listw(nlistk1);
    
   for(i in 1:dim(Inputdata)[2]){
   
  results    <- spatialreg::lagsarlm(Inputdata[,i] ~ as.factor(classids),listw = xy_weights)
  mt         <- summary(results); 
  SpMod_estm[i] <-  mt$Coef[2];
  SpMod_pval[i] <- mt$Coef[8];   
}
metlist <-list(SpMod_estm,SpMod_pval )
  }
  return(metlist)

}