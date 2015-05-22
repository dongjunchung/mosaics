
# extract posterior probability

.extractPostProb <- function( postProb, peakRegion, summaryStat="aveLogP", parallel, nCore ) {
  
  # split by chromosome for faster computation
  
  peakRegionChr <- split( peakRegion, peakRegion[,1] )
  locChr <- split( 1:nrow(peakRegion), peakRegion[,1] )
  postProbChr <- split( postProb, postProb[,1] )
  
  chrPeakRegion <- unique(peakRegion[,1])
  chrCommon <- intersect( chrPeakRegion, unique(postProb[,1]) )    
  
  binsize <- postProb[2,2] - postProb[1,2]
  
  # calculate summary statistics of posterior probabilities
  
  if ( parallel == TRUE ) {
    out <- mclapply( chrCommon, function(chrID) {      
      peakSel <- peakRegionChr[[ chrID ]]
      ppSel <- postProbChr[[ chrID ]]
      
      # extract post prob for peak regions
      
      rangePeak <- IRanges( start=peakSel[,2], end=peakSel[,3] )
      rangePP <- IRanges( start=ppSel[,2], end=ppSel[,2]+binsize-1 )
      
      # obtain pp for each peak region
      
      mm <- as.matrix(findOverlaps( rangePeak, rangePP ))
      mmlist <- split( mm[,2], mm[,1] )
      ppList <- lapply( mmlist, function(x) ppSel[x,3] )
      
      # calculate summary statistics of pp
      
      ppVecOrg <- switch( summaryStat,
        aveLogP = sapply( ppList, function(x) mean(-log10(x)) ),
        medianLogP = sapply( ppList, function(x) median(-log10(x)) ),
        sumLogP = sapply( ppList, function(x) sum(-log10(x)) ),
        logMinP = sapply( ppList, function(x) max(-log10(x)) ),
        logAveP = sapply( ppList, function(x) -log10(mean(x)) ),
        logMedianP = sapply( ppList, function(x) -log10(median(x)) )
      )
    
      # final version to return
      
      ppVec <- rep( NA, nrow(peakSel) )
      ppVec[ match( as.numeric(names(mmlist)), c(1:nrow(peakSel)) ) ] <- ppVecOrg
      
      return(ppVec)
    }, mc.cores=nCore )
  } else {
    out <- lapply( chrCommon, function(chrID) {      
      peakSel <- peakRegionChr[[ chrID ]]
      ppSel <- postProbChr[[ chrID ]]
      
      # extract post prob for peak regions
      
      rangePeak <- IRanges( start=peakSel[,2], end=peakSel[,3] )
      rangePP <- IRanges( start=ppSel[,2], end=ppSel[,2]+binsize-1 )
      
      # obtain pp for each peak region
      
      mm <- as.matrix(findOverlaps( rangePeak, rangePP ))
      mmlist <- split( mm[,2], mm[,1] )
      ppList <- lapply( mmlist, function(x) ppSel[x,3] )
      
      # calculate summary statistics of pp
      
      ppVecOrg <- switch( summaryStat,
        aveLogP = sapply( ppList, function(x) mean(-log10(x)) ),
        medianLogP = sapply( ppList, function(x) median(-log10(x)) ),
        sumLogP = sapply( ppList, function(x) sum(-log10(x)) ),
        logMinP = sapply( ppList, function(x) max(-log10(x)) ),
        logAveP = sapply( ppList, function(x) -log10(mean(x)) ),
        logMedianP = sapply( ppList, function(x) -log10(median(x)) )
      )
    
      # final version to return
      
      ppVec <- rep( NA, nrow(peakSel) )
      ppVec[ match( as.numeric(names(mmlist)), c(1:nrow(peakSel)) ) ] <- ppVecOrg
      
      return(ppVec)
    } )
  }
  
  # put the calculated PP back to the original location
  
  ppFinal <- data.frame( peakRegion, rep( NA, nrow(peakRegion) ) )  
  for ( chrID in chrCommon ) {
    ppFinal[ locChr[[chrID]], ncol(ppFinal) ] <- out[[ which( chrCommon == chrID ) ]]
  }
  colnames(ppFinal)[ ncol(ppFinal) ] <- summaryStat
  
  return(ppFinal)
}
