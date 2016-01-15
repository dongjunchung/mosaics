
.trimExtend <- function( chipExist, inputProvided, inputExist,
  chrID, peakStart, peakEnd, summitSignal, summit, stackedFragment, normC, 
  extendFromSummit, minRead, trimMinRead1, trimFC1, extendMinRead1, extendFC1,
  trimMinRead2, trimFC2, extendMinRead2, extendFC2 ) {
  
  if ( chipExist ) { 
      
    # Initialize with the originals
    
    peakStartNew <- peakStart
    peakEndNew <- peakEnd
    
    # load data
    
    xvar <- (stackedFragment$ChIP[[1]]):(stackedFragment$ChIP[[2]])
    yvar <- inverse.rle(stackedFragment$ChIP[[3]])    
    profileChip <- cbind( xvar, yvar )
    
    if ( inputProvided && inputExist && !is.na(stackedFragment$Input[[1]]) ) {
      xvar <- (stackedFragment$Input[[1]]):(stackedFragment$Input[[2]])
      yvar <- inverse.rle(stackedFragment$Input[[3]])    
      profileInput <- cbind( xvar, yvar )
    }
    
    # check summit, not restricted to peak region
    
    #summitOrg <- summit
    locSummit <- which( profileChip[,2] == max( profileChip[,2] ) )
    if ( length(locSummit) > 1 ) {
      firstPeak <- which.max(profileChip[,2])
      firstBreak <- locSummit[ which( diff(locSummit) != 1 )[1] ]
      if ( is.na(firstBreak) ) {
        firstBreak <- max(locSummit)
      }
      firstBlock <- firstPeak:firstBreak
      locSummit <- floor(mean(firstBlock))
    }
    summitNew <- profileChip[ locSummit, 1 ]
    if ( inputProvided && inputExist && !is.na(stackedFragment$Input[[1]]) ) {
      summitCoordInput <- match( summitNew, profileInput[,1] )
    }
    
    # check whether summit is located outside the boundary
    # if so, adjust start and end
    
    if ( summitNew < peakStart ) {
      peakStartNew <- summitNew
    }
    if ( peakEnd < summitNew ) {
      peakEndNew <- summitNew
    }
    
    # Extract ChIP & input profiles for the original and new peak region
    
    locStartNew <- match( peakStartNew, profileChip[,1] )
    if( is.na(locStartNew) ) {
      locStartNew <- 1
    }
    
    locEndNew <- match( peakEndNew, profileChip[,1] )
    if( is.na(locEndNew) ){
      locEndNew <- nrow(profileChip)
    }
    
    # check summit, while restricted to peak region
    
    profileChipRegion <- profileChip[ locStartNew:locEndNew, ]
    
    locSummitRegion <- which( profileChipRegion[,2] == max( profileChipRegion[,2] ) )
    if ( length(locSummitRegion) > 1 ) {
      firstPeak <- which.max(profileChipRegion[,2])
      firstBreak <- locSummitRegion[ which( diff(locSummitRegion) != 1 )[1] ]
      if ( is.na(firstBreak) ) {
        firstBreak <- max(locSummitRegion)
      }
      firstBlock <- firstPeak:firstBreak
      locSummitRegion <- floor(mean(firstBlock))
    }
    
    # calculate improvement of ChIP over input
    
    if ( inputProvided && inputExist && !is.na(stackedFragment$Input[[1]]) ) {
            
      profileCoordInput <- match( profileChip[,1], profileInput[,1] )
      profileInput <- cbind( profileChip[,1], profileInput[ profileCoordInput, 2 ] )
      profileInput[ is.na(profileInput[,2]), 2 ] <- 0
      chipInputFC <- 100 * ( profileChip[,2] - profileInput[,2] * normC ) / profileChip[,2]
      
      profileCoordInputRegion <- match( profileChipRegion[,1], profileInput[,1] )
      profileInputRegion <- cbind( profileChipRegion[,1], profileInput[ profileCoordInputRegion, 2 ] )
      profileInputRegion[ is.na(profileInputRegion[,2]), 2 ] <- 0
      chipInputFCRegion <- 100 * ( profileChipRegion[,2] - profileInputRegion[,2] * normC ) / profileChipRegion[,2]
      
    } else {
      
      profileInput <- profileChip
      profileInput[,2] <- 0
      chipInputFC <- 100
      
      profileInputRegion <- profileChipRegion
      profileInputRegion[,2] <- 0
      chipInputFCRegion <- 100
      
    }
    
    # Decide whether to trim
    # Trim: Trimming and extension are mostly minRead driven.
    
    trimStart <- FALSE
    trimEnd <- FALSE
    
    if ( profileChipRegion[1,2] < minRead * trimMinRead1 | chipInputFCRegion[1] < trimFC1 ) {
      trimStart <- TRUE
    }
    if ( profileChipRegion[ nrow(profileChipRegion), 2 ] < minRead * trimMinRead1 | chipInputFCRegion[ length(chipInputFCRegion) ] < trimFC1 ) {
      trimEnd <- TRUE
    }
    
    # Decide whether to extend
    
    extendStart <- FALSE
    extendEnd <- FALSE
    
    if ( profileChipRegion[1,2] >= minRead * extendMinRead1 & chipInputFCRegion[1] >= extendFC1 ) {
      extendStart <- TRUE
    }
    if ( profileChipRegion[ nrow(profileChipRegion), 2 ] >= minRead * extendMinRead1 & chipInputFCRegion[ length(chipInputFCRegion) ] >= extendFC1 ) { 
      extendEnd <- TRUE
    }
    
    # if needed, trim peak region
    
    if ( trimStart ) {
      i1 <- which( profileChipRegion[,2] >= minRead * trimMinRead2 & chipInputFCRegion >= trimFC2 )[1]
      #i2 <- which( chipInputFCRegion >= trimFC2 )[1]
      i2 <- NA
      
      if ( !is.na(i1) ) { 
        if ( i1 > locSummitRegion ) { i1 <- NA }
      }
      if ( !is.na(i2) ) {
        if ( i2 > locSummitRegion ) { i2 <- NA }
      }
      if ( !is.na(i1) | !is.na(i2) ) {
        peakStartNew <- profileChipRegion[ min( i1, i2, na.rm = TRUE ), 1 ] 
      }
    }
    
    if ( trimEnd ) {
      i1 <- rev(which( profileChipRegion[,2] >= minRead * trimMinRead2 & chipInputFCRegion >= trimFC2 ))[1]
      #i2<- rev(which( chipInputFCRegion >= trimFC2 ))[1]
      i2 <- NA
      
      if ( !is.na(i1) ) { 
        if ( i1 < locSummitRegion ) { i1 <- NA }
      }
      if ( !is.na(i2) ) {
        if ( i2 < locSummitRegion ) { i2 <- NA }
      }
      if ( !is.na(i1) | !is.na(i2) ) {
        peakEndNew <- profileChipRegion[ max( i1, i2, na.rm = TRUE ), 1 ]
      }
    }
    
    # if needed, extend peak region
    
    if ( extendStart ) {
      i1 <- which( profileChip[,2] >= minRead * extendMinRead2 & chipInputFC >= extendFC2 )[1]
      #i2 <- which( chipInputFC < FC & chipInputFC >=5 )[1]
      i2 <- NA
      
      if ( !is.na(i1) ) { 
        if ( i1 > locSummit ) { i1 <- NA }
      }
      if ( !is.na(i2) ) {
        if ( i2 > locSummit ) { i2 <- NA }
      }
      if ( !is.na(i1) | !is.na(i2) ) {
        peakStartNew <- profileChip[ max( i1, i2, na.rm = TRUE ), 1 ]
      }
    }
    
    if ( extendEnd ) {
      i1 <- rev(which( profileChip[,2] >= minRead * extendMinRead2 & chipInputFC >= extendFC2 ))[1]
      #i2 <- rev(which( chipInputFC < FC & chip.input.region.FC >= 5 ))[1]
      i2 <- NA 
      
      if ( !is.na(i1) ) { 
        if ( i1 < locSummit ) { i1 <- NA }
      }
      if ( !is.na(i2) ) {
        if ( i2 < locSummit ) { i2 <- NA }
      }
      if ( !is.na(i1) | !is.na(i2) ) {
        peakEndNew <- profileChip[ min( i1, i2, na.rm = TRUE ), 1 ]
      }
    }
    
    # If peakEndNew < peakStartNew, use original peak region
    
    if ( peakEndNew < peakStartNew ) {
      #message( paste( chrID,":",peakStart,"-",peakEnd, sep="" ) )
      #message( "Modified peakEnd is located on the left hand side of modified peakStart." )
      peakStartNew <- peakStart
      peakEndNew <- peakEnd
    }
    
    # If the start and stop still do not cover the summit, extend them
    
    if ( peakStartNew >= summitNew ) {
      #peakStartNew <- max( summitNew - extendFromSummit, 0, peakStartNew )
      peakStartNew <- max( min( summitNew - extendFromSummit, peakStartNew ), 0 )
    }
    if ( peakEndNew <= summitNew ) { 
      #peakEndNew <- min( summitNew + extendFromSummit, peakEndNew )
      peakEndNew <- max( summitNew + extendFromSummit, peakEndNew )
    }
    
    # update peak summit, while restricted to "new" peak region
    
    locStartNew <- match( peakStartNew, profileChip[,1] )
    if( is.na(locStartNew) ) {
      locStartNew <- 1
    }
    
    locEndNew <- match( peakEndNew, profileChip[,1] )
    if( is.na(locEndNew) ){
      locEndNew <- nrow(profileChip)
    }
    
    profileChipRegion <- profileChip[ locStartNew:locEndNew, , drop=FALSE ]    
    locSummitRegion <- which( profileChipRegion[,2] == max( profileChipRegion[,2] ) )
    if ( length(locSummitRegion) > 1 ) {
      firstPeak <- which.max(profileChipRegion[,2])
      firstBreak <- locSummitRegion[ which( diff(locSummitRegion) != 1 )[1] ]
      if ( is.na(firstBreak) ) {
        firstBreak <- max(locSummitRegion)
      }
      firstBlock <- firstPeak:firstBreak
      locSummitRegion <- floor(mean(firstBlock))
    }
    
    summitSignalNew <- profileChipRegion[ locSummitRegion, 2 ]
    summitNew <- profileChipRegion[ locSummitRegion, 1 ]
    
  } else {
    
    # if we don't have reads in peak region, no need to change peak boundaries or summits
    
    peakStartNew <- peakStart
    peakEndNew <- peakEnd
    summitSignalNew <- summitSignal
    summitNew <- summit
    
  }
  
  return( c( peakStartNew, peakEndNew, summitSignalNew, summitNew ) )
}

