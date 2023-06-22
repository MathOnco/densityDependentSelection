grpstats <- function(x,g,statscols,q1=0.5){
  allOut=list()
  o=matrix(NA, length(unique(g)),ncol(x) )
  rownames(o)=unique(g); colnames(o)=colnames(x);
  for(col in statscols){
    for(m in unique(g)){
      ii=which(g==m);
      v=NA;
      if(col=='mean'){
        v=apply(x[ii,,drop=F],2,mean,na.rm=T)
      }else if(col=='sum'){
        v=apply(x[ii,,drop=F],2,sum,na.rm=T)
      }else if(col=='var'){
        v=apply(x[ii,,drop=F],2,var)
      }else if(col=='max'){
        v=apply(x[ii,,drop=F],2,max,na.rm=T)
      }else if(col=='min'){
        v=apply(x[ii,,drop=F],2,min,na.rm=T)
      }else if(col=='quantile'){
        v=apply(x[ii,,drop=F],2,quantile,q1,na.rm=T)
      }else if(col=='median'){
        v=apply(x[ii,,drop=F],2,median,na.rm=T)
      }else if(col=='numel+'){##Count elements >0 
        v=apply(x[ii,,drop=F]>0,2,sum,na.rm=T)
      }else if(col=='numel_u'){##Count unique elements 
        v=apply(x[ii,,drop=F],2, function(k) length(unique(k)))
      }else if(col=='fraction+'){##Fraction of elements> 0 out of all finite elements
        v1=apply(!is.na(x[ii,,drop=F]),2,sum,na.rm=T)
        v=apply(x[ii,,drop=F]>0,2,sum,na.rm=T)/v1
      }else if(col=='maxcount'){##Most frequent value
        v1=plyr::count(x[ii])
        v=v1$x[which.max(v1$freq)]
      }else{
        v=get(col)(x[ii,,drop=F])
      }
      o[as.character(m),]=v;
    }
    allOut[[col]]=o;
  }
  return(allOut)
}


########################
# getTranscriptomeData #
########################

### This script returns the transcriptome data for cell line(s) of interest

getTranscriptomeData <- function(cellLines, ccstate="G0G1") {
  ## Load expression data along with clone membership:
  load(paste0("../data/seuratObjects/listOfSeurats_",ccstate,"_e.RObj"))
  genes = sapply(listOfSeurats[cellLines],function(x) rownames(x@assays$RNA@data))
  clones = sapply(listOfSeurats[cellLines],function(x) x@active.ident[colnames(x@assays$RNA@data)])
  cloneMembership = unlist(clones,use.names = F)
  ## Concatenate expression data across cell lines, use only genes present across all cell lines
  fr=plyr::count(unlist(genes))
  ex= sapply(listOfSeurats[cellLines],function(x) as.data.frame(as.matrix(x@assays$RNA@data))[fr$x,], USE.NAMES = F) 
  clMembership=unlist(sapply(names(ex), function(x) rep(x, ncol(ex[[x]])), simplify = F),use.names = F)
  ex= do.call(cbind, sapply(ex, as.matrix))
  ex[is.na(ex)]=0
  rownames(ex) = fr$x
  ## label by cell ID:
  names(clMembership) <- names(cloneMembership) <- unlist(sapply(clones, names))
  
  ## read gene sets

    table1 <- read.csv("../data/pathways/CTD_genes_pathways.CSV")
    table1$X..GeneSymbol = as.character(table1$X..GeneSymbol)
    unique_pathways <- unique(table1$PathwayName)
    gs <- list()
    for (name in unique_pathways) {
      ii=which(table1$PathwayName==name)
      gs[[name]] = unique(table1$X..GeneSymbol[ii])
    }
  
  ## Run AUC
  cells_rankings <- AUCell_buildRankings(ex,keepZeroesAsNA=T) 
  cells_AUC <- AUCell_calcAUC(gs, cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.1, normAUC = T)
  pq <- cells_AUC@assays@data@listData$AUC
  colnames(pq) = clMembership[colnames(pq)]
  ultimateList <- list(cloneMembership=cloneMembership, clMembership = clMembership, pathways=pq)
  
  return(ultimateList)
}



# This function takes as its input an object of class 'lm' and returns the p-value associated with the model

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
}


# This function takes an object of class "nonlinear_fit" from the package "growthrates" and returns the adjusted R-squared
getAdjRSqr <- function (model1) {
  if (class(model1) != "nonlinear_fit") stop("Not an object of class 'nonlinear_fit' ")
  adjR2 <- 1-(((1-model1@rsquared)*(length(model1@obs$y)-1))/(length(model1@obs$y)-(length(model1@par)-1)-1))
  return(adjR2)
}


getPathwayActivityRanksAcrossCellLines<-function(scDat){
  cellLineMedians = as.data.frame(grpstats(as.matrix(scDat$prediction), as.matrix(scDat$cellLine), "median"))
  cellLineMedians = cellLineMedians[order(cellLineMedians$median), ,drop = FALSE]
  cloneMedians = as.data.frame(grpstats(as.matrix(scDat$prediction), as.matrix(scDat$clone), "median"))
  cloneMedians = cloneMedians[order(cloneMedians$median), ,drop = FALSE]
  scDat$cellLineRank = match(scDat$cellLine, rownames(cellLineMedians))
  scDat$cloneRank = match(scDat$clone, rownames(cloneMedians))
  scDat = scDat[order(scDat$cloneRank), ]
  scDat = scDat[order(scDat$cellLineRank), ]
  scDat$PL_rank = 0
  for (i in 2:nrow(scDat)) {
    scDat$PL_rank[1] = 1
    if (scDat$clone[i] == scDat$clone[i - 1]) {
      scDat$PL_rank[i] = scDat$PL_rank[i - 1]
    } else{
      scDat$PL_rank[i] = scDat$PL_rank[i - 1] + 1
    }
  }
  return(scDat)
}

getLargestSignificantPair<-function(scDat, clonesToCellLinesMap, cellLinePairs, maxOnly=T, maxPval=0.05){
  ## find largest significant pair of cell lines and clones
  cellLineComps = list()
  for(i in 1:nrow(cellLinePairs)){
    pval=as.data.frame(t.test(scDat[scDat$cellLine==cellLinePairs[i,1], "prediction"], scDat[scDat$cellLine==cellLinePairs[i,2], "prediction"])$p.value)
    names(pval)="pval"
    pval$CLs = list(c(cellLinePairs[i,1],cellLinePairs[i,2]))
    rownames(pval) = list(c(cellLinePairs[i,1],cellLinePairs[i,2]))
    cellLineComps[[i]]= pval
  }
  cellLineComps = cellLineComps %>% map_df(~rbind(.x), .id = 'src')
  cellLineComps = cellLineComps[cellLineComps$pval<maxPval,]
  if(maxOnly){
    cellLineComps = cellLineComps[cellLineComps$pval==max(cellLineComps$pval),]
  }
  
  subpopComps = list()
  for(CL in unique(clonesToCellLinesMap$cellLine)){
    clonePairs = t(combn(unique(clonesToCellLinesMap[clonesToCellLinesMap$cellLine==CL, "clone"]), 2))
    for(i in 1:nrow(clonePairs)){
      if(numel(scDat[scDat$clone==clonePairs[i,1], "prediction"])>1 && numel(scDat[scDat$clone==clonePairs[i,2], "prediction"])>1){
        pval_=try(t.test(scDat[scDat$clone==clonePairs[i,1], "prediction"], scDat[scDat$clone==clonePairs[i,2], "prediction"]))
        if(class(pval_)=="try-error"){
          pval=as.data.frame(matrix(data=1))
          colnames(pval)="pval"
        }else{
          pval=as.data.frame(list(estimate=pval_$estimate[1]/pval_$estimate[2],pval=pval_$p.value))
        }
      }
      pval$clones = list(c(clonePairs[i,1],clonePairs[i,2]))
      rownames(pval) = list(c(clonePairs[i,1],clonePairs[i,2]))
      subpopComps[[CL]][[i]]= pval
    } 
    subpopComps[[CL]] = subpopComps[[CL]] %>% map_df(~rbind(.x), .id = 'src')
    subpopComps[[CL]] = subpopComps[[CL]][round(subpopComps[[CL]]$pval,3)<=maxPval,,drop=F]
    if(maxOnly && nrow(subpopComps[[CL]])>0){
      ii = order(subpopComps[[CL]]$pval,decreasing = T);
      subpopComps[[CL]] = subpopComps[[CL]][ii[1:min(3,length(ii))],,drop=F]
    }
  }
  subpopComps = subpopComps %>% map_df(~rbind(.x), .id = 'src')
  return(list("subpopComps"=subpopComps, "cellLineComps"=cellLineComps))
}

# Define function to calculate Rsquared and other stats for model fits for a given growth parameter
getStats4ParamFit<-function(param,uList,poi,verbose=F){
  ppCL=grpstats(t(uList$pathways), uList$clMembership, "median")$median[cellLines,]
  depVar = as.numeric(fitSummaryTable[rownames(ppCL),param])
  if(param == "K"){
    depVar = log(depVar)
  } 
  te=sapply(colnames(ppCL), function(x) lm(depVar~indepVar1, data=data.frame(indepVar1 = ppCL[,x])), simplify = F)
  te=te[poi$KEGG.Pathway]
  
  # get p-value for linear models
  pvals=sapply(names(te), function(x) lmp(te[[x]]))
  names(pvals)=gsub("\\.value*", "", names(pvals))
  model = te[[which.min(pvals)]]
  # get pearson's correlation for each pathway
  pearson=sapply(names(te), function(x) cor.test(ppCL[,x], as.numeric(fitSummaryTable[rownames(ppCL),param]))$estimate)
  names(pearson)=gsub("\\.cor*", "", names(pvals))
  # gather it all together
  te = as.data.frame(sapply(te, function(x) summary(x)$adj.r.squared))
  names(te) = "adj.r.squared"
  te$p.value = pvals[rownames(te)]
  te$pearson = pearson[rownames(te)]
  te=te[order(te$p.value),]
  if(verbose){
    print(param)
    print(head(te[1:5,]))
  }
  return(list(te=te, model=model))
}
