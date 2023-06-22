## getGrowthModel.R by Thomas Veith and Noemi Andor
## This script is a function which returns a nonlinear regression model summary and fitted plot
## The function takes as its input cell line lineage and model type

getGrowthModel<-function(cll,model) {
  # Import passaging table
  out=read.table(paste0("../data/cellCounts/", cll,".txt"))
  rownames(out)=out$id
  
  out$cellCount=out$correctedCount
  
  if(out$cellLine[1] %in% c("NCI-N87")){
    out$cellCount=out$areaOccupied_um2 /out$cellSize_um2
  }else if(out$cellLine[1] %in% c("SNU-668","HGC-27")){
    ## high confluence:
    ii = which(as.POSIXct(out$date)> quantile(as.POSIXct(out$date),0.75))
    out$cellCount[ii]=out$areaOccupied_um2[ii] /out$cellSize_um2[ii]
  }else if(out$cellLine[1] %in% c("SNU-601")){
    ## low confluence:
    ii = which(as.POSIXct(out$date)< quantile(as.POSIXct(out$date),0.5))
    out$cellCount[ii]=out$areaOccupied_um2[ii] /out$cellSize_um2[ii]
  }
  
  # Now we want to pull out just the date and cellCount columns of the out object to use for regression
  nDates<-as.numeric(as.POSIXct(out$date)) # need to convert dates from Objects of Class Date to Numeric
  nDates<- (nDates-min(nDates)) # normalize dates to start at Day 0
  nDates = nDates/(24*60^2)
  
  # set the parameters needed to perform regression and infer growth parameters
  
  y0 = out$cellCount[1]
  K = max(out$cellCount)
  p     <- c(y0 = y0, mumax = 0.2, K = K)
  lower <- c(y0 = y0/2, mumax = 0,   K = K/2)
  upper <- c(y0 = 2*y0, mumax = 5,   K = 2*K)
  
  
  # Here is where you build the regression model, and return the model with a fitted plot
  model1 <- fit_growthmodel(FUN = model, p = p, nDates, out$cellCount, lower = lower, upper = upper)
  adjR2 <- getAdjRSqr(model1)
  # growthrates::plot(model1)
  
  
  newList<-list(modelsum=model1, out=out, adjr2=adjR2)
  
  return(newList)
}

