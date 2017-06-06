
## Extract reads corresponding to each peak region

.loadReadData <- function( object, readfile=NULL,
  fileFormat=NULL, PET=FALSE, fragLen=200, keepReads=FALSE,
  parallel=FALSE, nCore=8, tempDir=NULL, perl="perl" )
{ 
  
  # process aligned read file
  
  if ( PET ) {
    message( "Info: Paired-end tag (PET) is assumed (PET=TRUE)." )
  } else {
    message( "Info: Single-end tag (SET) is assumed (PET=FALSE)." )
    message( "Info: Average fragment length is set as ",fragLen," (fragLen=",fragLen,")." )
  }
  #print( Sys.time() )
  
  if ( PET == TRUE ) {
      allFormat <- c( "eland_result", "bam" )
      #allFormat <- c( "eland_result", "sam", "bam" )
      allFormatName <- c( "Eland result", "BAM" )
      #allFormatName <- c( "Eland result", "SAM", "BAM" )
  } else {
      allFormat <- c( "eland_result", "eland_extended", "eland_export", 
          "bowtie", "sam", "bam", "bed" )
      #    "bowtie", "sam", "bam", "bed", "csem" )
      allFormatName <- c( "Eland result", "Eland extended", "Eland export", 
          "Bowtie default", "SAM", "BAM", "BED" )
      #    "Bowtie default", "SAM", "BAM", "BED", "CSEM" )
  }
        
  # check whether minimal options are missing
  
  if ( length(readfile) != 1 || is.null(readfile) )
  {
      stop( "Please specify the aligned read file!" )
  }   
  
  if ( length(fileFormat) != 1 || is.null(fileFormat) )
  {
      stop( "Please specify aligned read file format! Read '?constructExtRead' for supported file formats" )
  }   
  
  # check file format specification
  
  if ( length(which(!is.na(match( fileFormat, allFormat )))) == 0 )
  {
      stop( "Unsupported aligned read file format! Read '?constructExtRead' for supported file formats" )
  } 
  
  fileFormatName <- allFormatName[ match( fileFormat, allFormat ) ]

  # summarize peak info
  
  peakList <- print(object)
  nPeak <- nrow(peakList)
  
  locByChr <- split( 1:nrow(peakList), peakList[,1] )
  peakByChr <- split( peakList[ , 2:3, drop=FALSE ], peakList[,1] )
  peakChr <- names(peakByChr)
  
  #chrCommon <- Reduce( intersect, list( peakChr, readChr ) )
  chrCommon <- peakChr
  chrCommon <- sort(chrCommon)
  
  nameVecAll <- 
    paste( peakList[,1], ":", peakList[,2], "-", peakList[,3], sep="" )
  nameVecChr <- split( nameVecAll, peakList[,1] )
  
  # construct GRanges for peak list, while extending boundaries little bit for safety

  peakgr <- GRanges( 
    seqnames = Rle(peakList[,1]), 
    ranges = IRanges( start=peakList[,2], end=peakList[,3], names=nameVecAll ),
    strand = Rle(strand(rep("*",nrow(peakList))))
  )
  peakgr <- split( peakgr, seqnames(peakgr) )
  
  # construct GRanges for extended peak list, for safer mapping

  peakgrExt <- GRanges( 
    seqnames = Rle(peakList[,1]), 
    ranges = IRanges( start=peakList[,2] - fragLen, end=peakList[,3] + fragLen, names=nameVecAll ),
    strand = Rle(strand(rep("*",nrow(peakList))))
  )
  peakgrExt <- trim(peakgrExt)
  peakgrExt <- split( peakgrExt, seqnames(peakgrExt) )
  
  # process read files
  
  if ( fileFormat != "bam" ) {
    
    #####################################################
    #
    # If not BAM, use perl script to process data
    #
    #####################################################
    
    message( "Info: Reading and processing aligned read file..." )
    
    if ( is.null(tempDir) ) {
      tempfileName <- tempfile( c("output","summary") )
    } else {
      tempfileName <- c( paste(tempDir,"output.txt",sep=""),
        paste(tempDir,"summary.txt",sep="") )
    }
    
    # intermediate file name
    
    .constructExtRead( infile=readfile, outfile=tempfileName[1],
      summaryfile=tempfileName[2], fileFormat=fileFormat, 
      PET=PET, fragLen=fragLen, perl=perl )            
    
    # read summary file (chrID, # lines)
    
    summaryInfo <- read.table( tempfileName[2], header=FALSE,
      stringsAsFactors=FALSE, comment.char="", check.names=FALSE )
    colnames(summaryInfo) <- c("chrID","nline")   
    readChr <- summaryInfo[,1] 
    #print( Sys.time() )
    
    # match reads to each peak region for each chromosome and construct GRangesList
    # (using parallel computing, if parallel exists)
    
    if ( parallel == TRUE ) {        
      greads <- mclapply( chrCommon, 
        function(x) .extread2GRanges( chr=x, 
          outfileName=paste(tempfileName[1],"_",x,sep=""), 
          nRow=summaryInfo[ summaryInfo[,1]==x, 2 ], PET=PET ), 
        mc.cores = nCore )
    } else {        
      greads <- lapply( chrCommon, 
        function(x) .extread2GRanges( chr=x, 
          outfileName=paste(tempfileName[1],"_",x,sep=""), 
          nRow=summaryInfo[ summaryInfo[,1]==x, 2 ], PET=PET )
      )
    }
  
	  greads <- GRangesList(greads)
    
    # name of list
    
    namevec <- rep( "", length(greads) )
    for ( i in 1:length(greads) ) {
      namevec[i] <- as.character(seqnames(greads[[i]])[1])
    }
    names(greads) <- namevec
    
    # sequencing depth
    
    seqDepth <- sum(summaryInfo[,2])
    
    # remove temporary files after use
    
    unlink( tempfileName[1] )
    unlink( tempfileName[2] )
    
    gc()    
    
  } else {
    
    #####################################################
    #
    # In the case of BAM, use BioC packages to process data
    #
    #####################################################
    
    # check whether BAM index exists. Otherwise, generate BAM index
    
    bamName <- list.files( path=dirname(readfile), pattern=basename(readfile) )
    if ( length(grep( "bai", bamName )) > 0 ) {
      message( "Info: Use the provided BAM index file." )
    } else {
      message( "Info: BAM index file does not exist. Generating BAM index file..." )
      indexBam( readfile )
    }
    
    # load reads corresponding to peak regions
    
    message( "Info: Reading and processing aligned read file..." )
  
	  param <- ScanBamParam( which=peakgrExt )
    
	  if ( PET == FALSE ) {
		  suppressWarnings( greads <- readGAlignments( readfile, param = param, use.names = FALSE ) )
		  suppressWarnings( greads <- as( greads, "GRanges" ) )
		  suppressWarnings( greads <- resize( greads, fragLen ) )
	  } else {
		  #suppressWarnings( greads <- readAlignmentsPairsFromBam( readfile, param = param ) )
		  #suppressWarnings( greads <- GRanges( seqnames = seqnames(greads),
		  #	ranges = IRanges( start=start(left(greads)), end=end(right(greads)) ),
			#  strand = Rle( "*", length(greads) ) )
			#)
      
      suppressWarnings( greads <- readGAlignmentPairs( readfile, param = param ) )

      snms = seqnames(greads)
      starts = start(greads@first)
      ends = end(greads@last)
      
      # remove reads with negative widths         
      idx = (starts >= ends)
      if(any(idx)){
        warning("Removing ",sum(idx)," reads, due to negative read lengths")
        snms = snms[!idx]
        starts = starts[!idx]
        ends = ends[!idx]
      }
      suppressWarnings( greads <- GRanges( seqnames = snms,ranges = IRanges( start=starts, end=ends),strand = "*"))
      rm(snms,starts,ends)
	  }
    
    # split by chromosome
    
    greads <- split( greads, seqnames(greads) )
    
    # sequencing depth
    
    seqDepth <- countBam(readfile)$records
  }
    
  # match reads with peaks & calculate coverage
  
  message( "Info: Processing and combining peak list and reads..." )    
  
  if ( parallel == TRUE ) {        
    fragSetOrg <- mclapply( chrCommon, function(chr) {
      suppressWarnings( greads[[chr]] <- trim(greads[[chr]]) )
      #suppressWarnings( midp <- resize( shift( greads[[chr]], width(greads[[chr]])/2 ), 1 ) )
      midp <- shift( greads[[chr]], width(greads[[chr]])/2 )
      strand(midp) <- "*"
      midp <- resize( midp, 1 )
      overlaps <- findOverlaps( peakgr[[chr]], midp )
      suppressWarnings( fragSetEach <- split( 
        greads[[chr]][subjectHits(overlaps)], queryHits(overlaps) ) )
      idFound <- unique(queryHits(overlaps))
      names(fragSetEach) <- paste( seqnames(peakgr[[chr]]), ":", 
        start(peakgr[[chr]]), "-", end(peakgr[[chr]]), sep="" )[ idFound ]    
      return(fragSetEach)      
    }, mc.cores = nCore )
  } else {        
    fragSetOrg <- lapply( chrCommon, function(chr) {
      suppressWarnings( greads[[chr]] <- trim(greads[[chr]]) )
      #suppressWarnings( midp <- resize( shift( greads[[chr]], width(greads[[chr]])/2 ), 1 ) )
      midp <- shift( greads[[chr]], width(greads[[chr]])/2 )
      strand(midp) <- "*"
      midp <- resize( midp, 1 )
      overlaps <- findOverlaps( peakgr[[chr]], midp )
      suppressWarnings( fragSetEach <- split( 
        greads[[chr]][subjectHits(overlaps)], queryHits(overlaps) ) )
      idFound <- unique(queryHits(overlaps))
      names(fragSetEach) <- paste( seqnames(peakgr[[chr]]), ":", 
        start(peakgr[[chr]]), "-", end(peakgr[[chr]]), sep="" )[ idFound ]        
      return(fragSetEach)
    } )
  }
  
  names(fragSetOrg) <- chrCommon

  rm( greads )
  gc()
  
  # flatten reads to match them to peak list for easier processing in later steps

  if ( parallel == TRUE ) {
    fragSet <- mclapply( as.list(1:nPeak), function(i) {
      return(fragSetOrg[[ peakList[i,1] ]][[ nameVecAll[i] ]])
    }, mc.cores = nCore )
  } else {
    fragSet <- lapply( as.list(1:nPeak), function(i) {
      return(fragSetOrg[[ peakList[i,1] ]][[ nameVecAll[i] ]])
    } )
  }
  names(fragSet) <- nameVecAll

  rm( fragSetOrg )
  gc()
  
  # calculate number of reads
  
  numReads <- sapply( fragSet, length )
  names(numReads) <- names(fragSet)
  
  # calculate coverage
  
  message( "Info: Calculating coverage..." )    
  
  if ( parallel == TRUE ) {
    stackedFragment <- mclapply( fragSet, function(fragEach) {
      outmat <- vector( "list", 3 )
      
      if ( length(fragEach) > 0 ) {
        fragmat <- cbind( start(fragEach), end(fragEach) )
        xvarq <- c( min(fragmat[,1]), max(fragmat[,2]) )
        xvar <- c(xvarq[1]:xvarq[2])
        yvar <- .ff_stack( fragmat[,1], fragmat[,2], xvarq[1], xvarq[2] )  
        #outmat <- cbind( xvar, yvar )
        outmat[[1]] <- min(xvar)
        outmat[[2]] <- max(xvar)
        outmat[[3]] <- rle(yvar)
      } else {
        #outmat <- matrix( NA )
        outmat[[1]] <- outmat[[2]] <- outmat[[3]] <- NA
      }
      return( outmat )
    }, mc.cores = nCore )
  } else {
    stackedFragment <- lapply( fragSet, function(fragEach) {
      outmat <- vector( "list", 3 )
      
      if ( length(fragEach) > 0 ) {
        fragmat <- cbind( start(fragEach), end(fragEach) )
        xvarq <- c( min(fragmat[,1]), max(fragmat[,2]) )
        xvar <- c(xvarq[1]:xvarq[2])
        yvar <- .ff_stack( fragmat[,1], fragmat[,2], xvarq[1], xvarq[2] )  
        #outmat <- cbind( xvar, yvar )
        outmat[[1]] <- min(xvar)
        outmat[[2]] <- max(xvar)
        outmat[[3]] <- rle(yvar)
      } else {
        #outmat <- matrix( NA )
        outmat[[1]] <- outmat[[2]] <- outmat[[3]] <- NA
      }
      return( outmat )
    } )
  }
  names(stackedFragment) <- names(fragSet)
  
  gc()
  
  message( "Info: Done!\n" )
  
  # return output
  
  if ( keepReads == TRUE ) {
    return(list( 
      fragSet=fragSet, numReads=numReads, stackedFragment=stackedFragment, seqDepth=seqDepth ))
  } else {
    return(list( 
      fragSet=GRanges(), numReads=numReads, stackedFragment=stackedFragment, seqDepth=seqDepth ))
  }
}

# match peak & reads (other file formats)

.extread2GRanges <- function( chr, outfileName, nRow, PET ) {
  
  # read processed read file
  # - PET read: chr, start, end
  # - SET read: chr, position, strand, read length
  
  if ( PET == TRUE ) {
    # if PET, (chrID, start, end)
    
    readCur <- read.table( outfileName, sep='\t', nrows=nRow,
      header=FALSE, stringsAsFactors=FALSE, comment.char="",
      colClasses=c("character","numeric","numeric"), check.names=FALSE )
    colnames(readCur) <- c("chrID","start","end")
	
    greads <- GRanges( seqnames = Rle(readCur[,1]),
    	ranges = IRanges( start=readCur[,2], end=readCur[,3] ),
    	strand = Rle( "*", nrow(readCur) )
	  )
  } else {    
    # if SET, (chrID, start, end, strand)
    
    readCur <- read.table( outfileName, sep='\t', nrows=nRow,
      header=FALSE, stringsAsFactors=FALSE, comment.char="",
      colClasses=c("character","numeric","numeric","character"), 
      check.names=FALSE )
    colnames(readCur) <- c("chrID","start","end","str")
    
    readCur[ readCur[,4] == "F", 4 ] <- "+"
    readCur[ readCur[,4] == "R", 4 ] <- "-"
	
  	greads <- GRanges( seqnames = Rle(readCur[,1]),
  		ranges = IRanges( start=readCur[,2], end=readCur[,3] ),
  		strand = Rle(readCur[,4])
  	)
	
  }
	
  return(greads)
}
