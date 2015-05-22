.annotatePeak <- function( peakStart_chr, peakStop_chr, coord_chr, analysisType,
  Y_chr, X_chr, M_chr, GC_chr, pp_chr, nRatio=1 ) {
  
  # calculate additional info
      
  peaksize <- peakStop_chr - peakStart_chr + 1
  binsize <- coord_chr[2] - coord_chr[1]
  
  peak_range <- IRanges( start=peakStart_chr, end=peakStop_chr )
  info_range <- IRanges( start=coord_chr, end=coord_chr+binsize-1 )
  mm <- as.matrix( findOverlaps( peak_range, info_range ) )
  matchlist <- split( mm[,2], mm[,1] )
  
  switch( analysisType,
      OS = {        
          # extract info
          
          peak_info <- sapply( matchlist, 
              function(x) {
                  betapH_i <- pp_chr[x]
                  logAveP <- -log10(mean(betapH_i))
                  logMinP <- -log10(min(betapH_i))
                  aveLogP <- mean(-log10(betapH_i))
                  aveChipCount <- mean(Y_chr[x])
                  maxChipCount <- max(Y_chr[x])
                  aveM <- mean(M_chr[x])
                  aveGC <- mean(GC_chr[x]) 
                  return( c( logAveP, logMinP, aveLogP, aveChipCount, maxChipCount, aveM, aveGC ) )
              }
          )
              
          peak_info <- t(peak_info)
          
          # combine all
                  
          final_peakset_chr <- data.frame(
              peakStart_chr, peakStop_chr, peaksize, peak_info, 
              stringsAsFactors=FALSE )
          colnames(final_peakset_chr) <-
              c('peakStart','peakStop','peakSize','logAveP','logMinP','aveLogP',
              'aveChipCount','maxChipCount','map','GC')
      },
      TS = {        
          # extract info
      
          peak_info <- sapply( matchlist, 
              function(x) {
                  betapH_i <- pp_chr[x]
                  logAveP <- -log10(mean(betapH_i))
                  logMinP <- -log10(min(betapH_i))
                  aveLogP <- mean(-log10(betapH_i))
                  aveChipCount <- mean(Y_chr[x])
                  maxChipCount <- max(Y_chr[x])
                  aveInputCount <- mean(X_chr[x])
                  aveInputCountScaled <- aveInputCount * nRatio
                  aveLog2Ratio <- mean( log2( (Y_chr[x]+1) / (X_chr[x]*nRatio+1) ) )
                  aveM <- mean(M_chr[x])
                  aveGC <- mean(GC_chr[x]) 
                  return( c( logAveP, logMinP, aveLogP, aveChipCount, maxChipCount, 
                      aveInputCount, aveInputCountScaled, aveLog2Ratio, aveM, aveGC ) )
              }
          )
              
          peak_info <- t(peak_info)
          
          # combine all
                  
          final_peakset_chr <- data.frame( 
              peakStart_chr, peakStop_chr, peaksize, peak_info, 
              stringsAsFactors=FALSE )
          colnames(final_peakset_chr) <-
              c('peakStart','peakStop','peakSize','logAveP','logMinP','aveLogP',
              'aveChipCount','maxChipCount',
              'aveInputCount','aveInputCountScaled','aveLog2Ratio','map','GC')
      },
      IO = {        
          # extract info
          
          peak_info <- sapply( matchlist, 
              function(x) {
                  betapH_i <- pp_chr[x]
                  logAveP <- -log10(mean(betapH_i))
                  logMinP <- -log10(min(betapH_i))
                  aveLogP <- mean(-log10(betapH_i))
                  aveChipCount <- mean(Y_chr[x])
                  maxChipCount <- max(Y_chr[x])
                  aveInputCount <- mean(X_chr[x])
                  aveInputCountScaled <- aveInputCount * nRatio
                  aveLog2Ratio <- mean( log2( (Y_chr[x]+1) / (X_chr[x]*nRatio+1) ) )
                  return( c( logAveP, logMinP, aveLogP, aveChipCount, maxChipCount, 
                      aveInputCount, aveInputCountScaled, aveLog2Ratio ) )
              }
          )
              
          peak_info <- t(peak_info)
      
          # combine all
                  
          final_peakset_chr <- data.frame( 
              peakStart_chr, peakStop_chr, peaksize, peak_info, 
              stringsAsFactors=FALSE )
          colnames(final_peakset_chr) <-
              c('peakStart','peakStop','peakSize','logAveP','logMinP','aveLogP',
              'aveChipCount','maxChipCount',
              'aveInputCount','aveInputCountScaled','aveLog2Ratio')
      }
  )
	
	return(final_peakset_chr)
}
