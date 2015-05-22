
# read bin-level data & process it

readBins <- function( type=c("chip","input"), fileName=NULL, 
    dataType='unique', rounding=100, parallel=FALSE, nCore=8 )
{
    # [Note] Assumption: chip, input, M, GC, & N have corresponding coordinates 
    
    # error treatment
    
    .validType( type )
    .validLocation( type=type, fileName=fileName )
    .validDataType( dataType )
    
    # check options: parallel computing (optional)
    
    if ( parallel == TRUE ) {
        message( "Use 'parallel' package for parallel computing." )        
        if ( length(find.package('parallel',quiet=TRUE)) == 0 ) {
            stop( "Please install 'parallel' package!" )
        }
    }
    
    # existence
    
    existInput <- ( length(which(type=="input")) > 0 )
    existM <- ( length(which(type=="M")) > 0 )
    existGC <- ( length(which(type=="GC")) > 0 )
    existN <- ( length(which(type=="N")) > 0 )
    
    # read bin-level data
    
    message( "Info: reading and preprocessing bin-level data..." )
    
    chipFileName <- fileName[ type=="chip" ]
    chip <- read.table( chipFileName, header=FALSE, sep='\t', 
        colClasses=c("character","numeric","numeric"), comment.char="#",
        stringsAsFactors=FALSE, check.names=FALSE )
    
    if ( existM )
    {
        mapScoreFileName <- fileName[ type=="M" ]
        mapScore <- read.table( mapScoreFileName, header=FALSE, 
            stringsAsFactors=FALSE, check.names=FALSE, comment.char="#" )
    }
    if ( existGC ) {    
        gcScoreFileName <- fileName[ type=="GC" ]
        gcScore <- read.table( gcScoreFileName, header=FALSE, 
            stringsAsFactors=FALSE, check.names=FALSE, comment.char="#" )
    }   
    if ( existN ) {
        nNucFileName <- fileName[ type=="N" ]
        nNuc <- read.table( nNucFileName, header=FALSE, 
            stringsAsFactors=FALSE, check.names=FALSE, comment.char="#" )
    }    
    if ( existInput )
    {
        inputFileName <- fileName[ type=="input" ]
        input <- read.table( inputFileName, header=FALSE, sep='\t', 
            colClasses=c("character","numeric","numeric"), comment.char="#",
            stringsAsFactors=FALSE, check.names=FALSE )
    }
    
    # error treatment
    
    if ( existInput )
    {
        if ( chipFileName == inputFileName ) {
            warning( "the same data was used for both ChIP & Input." )
            warning( "parameters cannot be properly estimated in this case." )
            warning( "remember that you will get errors in 'mosaicsFit' method!" )
        }
    } 
  
    # extract sequencing depth
  
    seqDepth <- rep( NA, 2 )
    seqDepth[1] <- read.table( file=chipFileName, nrows=1, comment.char="" )[1,4]
    if ( existInput ) {
      seqDepth[2] <- read.table( file=inputFileName, nrows=1, comment.char="" )[1,4]
    } else {
      seqDepth[2] <- NA
    }
    
    
    # process bin-level data
    
    chrChIP <- unique(chip[,1])
    chrCommon <- chrChIP
    
    if ( existInput )
    {        
        if ( existM & existGC & existN )
        {            
            ####################################################################
            #                                                                  #
            #     two-sample analysis (if M, GC, N are available)              #
            #                                                                  #
            ####################################################################
                        
            if ( length(chrChIP) > 1 ) {
                # genome-wide analysis
                
                message( "Info: assume that data contains more than one chromosome" )
                
                if ( ncol(mapScore) == 3 & ncol(gcScore) == 3 & ncol(nNuc) == 3 ) {                
                    # split data by chromosome
                    
                    chipByChr <- split( chip[,2:3], chip[,1] )
                    inputByChr <- split( input[,2:3], input[,1] )
                    mapByChr <- split( mapScore[,2:3], mapScore[,1] )
                    gcByChr <- split( gcScore[,2:3], gcScore[,1] )
                    nByChr <- split( nNuc[,2:3], nNuc[,1] )
                    
                    # extract list of shared chromosomes
                    
                    chrInput <- unique(input[,1])
                    chrM <- unique(mapScore[,1])
                    chrGC <- unique(gcScore[,1])
                    chrN <- unique(nNuc[,1])
                    
                    chrCommon <- Reduce( intersect,
                        list( chrChIP, chrInput, chrM, chrGC, chrN ) )
                    chrCommon <- sort(chrCommon)
                    
                    # provide error & stop if there is no common chromosome
                    
                    if ( length(chrCommon) == 0 ) {
                        stop( "No chromosome is common among ChIP, control, mappability, GC, and N files." )
                    }
                    
                    # construct data to process

                    coordMat <- matrix( NA, length(chrCommon), 4 )
                    
                    dataList <- list()
                    for ( i in 1:length(chrCommon) ) {
                        dataList[[i]] <- list()
                        dataList[[i]]$chip <- chipByChr[[ chrCommon[i] ]]
                        dataList[[i]]$input <- inputByChr[[ chrCommon[i] ]]
                        dataList[[i]]$mapScore <- mapByChr[[ chrCommon[i] ]]
                        dataList[[i]]$gcScore <- gcByChr[[ chrCommon[i] ]]
                        dataList[[i]]$nNuc <- nByChr[[ chrCommon[i] ]]
                        
                        coordMat[i,1] <- min(chipByChr[[ chrCommon[i] ]][,1])
                        coordMat[i,2] <- max(chipByChr[[ chrCommon[i] ]][,1])
                    }
                    
                    rm( chipByChr, inputByChr, mapByChr, gcByChr, nByChr )
                    rm( chrInput, chrM, chrGC, chrN )
                    gc()
                    
                    # processing (using parallel computing, if multicore exists)
                    
                    if ( parallel == TRUE ) {
                        # if "multicore" package exists, utilize parallel computing with "mclapply"
                        #require(multicore)
                        
                        out <- mclapply( dataList, function(x) { 
                            .processBin_MGCX( chip=x$chip, input=x$input,
                                mapScore=x$mapScore, gcScore=x$gcScore, nNuc=x$nNuc, 
                                dataType=dataType, rounding=rounding )
                            }, mc.cores = nCore )
                    } else {
                        # otherwise, use usual "lapply"
                        
                        out <- lapply( dataList, function(x) { 
                            .processBin_MGCX( chip=x$chip, input=x$input,
                                mapScore=x$mapScore, gcScore=x$gcScore, nNuc=x$nNuc, 
                                dataType=dataType, rounding=rounding )
                            } )
                    }    
                    
                    rm( dataList )
                    gc()
                    
                    # stack processed data
                    
                    chrID <- coord <- Y <- X <- M <- GC <- c()
                    
                    for ( i in 1:length(chrCommon) ) {
                        chrID <- c( chrID, rep( chrCommon[i], length(out[[i]]$coord) ) )
                        coord <- c( coord, out[[i]]$coord )
                        Y <- c( Y, out[[i]]$Y )
                        X <- c( X, out[[i]]$X )
                        M <- c( M, out[[i]]$M )
                        GC <- c( GC, out[[i]]$GC )
                        
                        coordMat[i,3] <- min(out[[i]]$coord)
                        coordMat[i,4] <- max(out[[i]]$coord)
                    }               
                    
                    rm( out )
                    gc()   
    
                    outBin <- new( "BinData", chrID=chrID, coord=coord, tagCount=Y, input=X, 
                        mappability=M, gcContent=GC, dataType=dataType, seqDepth=seqDepth )
                        
                    rm( chrID, coord, Y, X, M, GC )
                    gc()    
                } else {
                    stop( "Numbers of columns do not match!" )
                }
            } else {
                message( "Info: assume that data contains only one chromosome" )
                
                if ( ncol(mapScore) == 2 & ncol(gcScore) == 2 & ncol(nNuc) == 2 ) {
                    # processing
                    
                    out <- .processBin_MGCX( chip=chip[,2:3], input=input[,2:3],
                        mapScore=mapScore, gcScore=gcScore, nNuc=nNuc, 
                        dataType=dataType, rounding=rounding )                        
                    
                    chrID <- rep( chip[1,1], length(out$coord) )
                    outBin <- new( "BinData", chrID=chrID, 
                        coord=out$coord, tagCount=out$Y, input=out$X, 
                        mappability=out$M, gcContent=out$GC, dataType=dataType, seqDepth=seqDepth )  
                        
                    rm( chrID, out )
                    gc()
                } else {
                    stop( "Numbers of columns do not match!" )
                }   
            }         
        } else if ( existN ) {
            ####################################################################
            #                                                                  #
            #          two-sample analysis (input only, with N)                #
            #                                                                  #
            ####################################################################
                        
            if ( length(chrChIP) > 1 ) {
                # genome-wide analysis
                
                message( "Info: data contains more than one chromosome" )
                
                if ( ncol(nNuc) == 3 ) {                
                    # split data by chromosome
                    
                    chipByChr <- split( chip[,2:3], chip[,1] )
                    inputByChr <- split( input[,2:3], input[,1] )
                    nByChr <- split( nNuc[,2:3], nNuc[,1] )
                    
                    # extract list of shared chromosomes
                    
                    chrInput <- unique(input[,1])
                    chrN <- unique(nNuc[,1])                    
                    chrCommon <- intersect( chrChIP, chrInput, chrN )
                    chrCommon <- sort(chrCommon)
                    
                    # provide error & stop if there is no common chromosome
                    
                    if ( length(chrCommon) == 0 ) {
                        stop( "No chromosome is common among ChIP, control, and N files." )
                    }
                    
                    # construct data to process
                    
                    coordMat <- matrix( NA, length(chrCommon), 4 )
                    
                    dataList <- list()
                    for ( i in 1:length(chrCommon) ) {
                        dataList[[i]] <- list()
                        dataList[[i]]$chip <- chipByChr[[ chrCommon[i] ]]
                        dataList[[i]]$input <- inputByChr[[ chrCommon[i] ]]
                        dataList[[i]]$nNuc <- nByChr[[ chrCommon[i] ]]
                        
                        coordMat[i,1] <- min(chipByChr[[ chrCommon[i] ]][,1])
                        coordMat[i,2] <- max(chipByChr[[ chrCommon[i] ]][,1])
                    }
                    
                    rm( chipByChr, inputByChr, nByChr )
                    rm( chrInput, chrN )
                    gc()
                    
                    # processing (using parallel computing, if multicore exists)
                    
                    if ( parallel == TRUE ) {
                        # if "multicore" package exists, utilize parallel computing with "mclapply"
                        #require(multicore)
                        
                        out <- mclapply( dataList, function(x) { 
                            .processBin_XN( chip=x$chip, input=x$input, nNuc=x$nNuc, 
                                dataType=dataType )
                            }, mc.cores = nCore )
                    } else {
                        # otherwise, use usual "lapply"
                        
                        out <- lapply( dataList, function(x) { 
                            .processBin_XN( chip=x$chip, input=x$input, nNuc=x$nNuc, 
                                dataType=dataType )
                            } )
                    }    
                    
                    rm( dataList )
                    gc()
                    
                    # stack processed data
                    
                    chrID <- coord <- Y <- X <- c()
                    
                    for ( i in 1:length(chrCommon) ) {
                        chrID <- c( chrID, rep( chrCommon[i], length(out[[i]]$coord) ) )
                        coord <- c( coord, out[[i]]$coord )
                        Y <- c( Y, out[[i]]$Y )
                        X <- c( X, out[[i]]$X )
                        
                        coordMat[i,3] <- min(out[[i]]$coord)
                        coordMat[i,4] <- max(out[[i]]$coord)
                    }           
                    
                    rm( out )
                    gc()       
    
                    outBin <- new( "BinData", chrID=chrID, coord=coord, tagCount=Y, input=X, 
                        dataType=dataType, seqDepth=seqDepth )
                    
                    rm( chrID, coord, Y, X )
                    gc()    
                } else {
                    stop( "Numbers of columns do not match!" )
                }
            } else {
                message( "Info: data contains only one chromosome" )
                
                if ( ncol(nNuc) == 2 ) {
                    # processing
                    
                    out <- .processBin_XN( chip=chip[,2:3], input=input[,2:3], nNuc=nNuc, 
                        dataType=dataType )
    
                    chrID <- rep( chip[1,1], length(out$coord) ) 
                    outBin <- new( "BinData", chrID=chrID, coord=out$coord, tagCount=out$Y, 
                        input=out$X, dataType=dataType, seqDepth=seqDepth ) 
                    
                    rm( chrID, out )
                    gc() 
                } else {
                    stop( "Numbers of columns do not match!" )
                }   
            }       
        } else
        {
            ####################################################################
            #                                                                  #
            #            two-sample analysis (input only, without N)           #
            #                                                                  #
            ####################################################################
                        
            if ( length(chrChIP) > 1 ) {
                # genome-wide analysis
                
                message( "Info: data contains more than one chromosome." )
                           
                # split data by chromosome
                
                chipByChr <- split( chip[,2:3], chip[,1] )
                inputByChr <- split( input[,2:3], input[,1] )
                
                # extract list of shared chromosomes
                
                chrInput <- unique(input[,1])                 
                chrCommon <- intersect( chrChIP, chrInput )
                chrCommon <- sort(chrCommon)
                    
                # provide error & stop if there is no common chromosome
                
                if ( length(chrCommon) == 0 ) {
                    stop( "No chromosome is common between ChIP and control files." )
                }
                
                # construct data to process
                
                coordMat <- matrix( NA, length(chrCommon), 4 )
                
                dataList <- list()
                for ( i in 1:length(chrCommon) ) {
                    dataList[[i]] <- list()
                    dataList[[i]]$chip <- chipByChr[[ chrCommon[i] ]]
                    dataList[[i]]$input <- inputByChr[[ chrCommon[i] ]]
                        
                    coordMat[i,1] <- min(chipByChr[[ chrCommon[i] ]][,1])
                    coordMat[i,2] <- max(chipByChr[[ chrCommon[i] ]][,1])
                }
                    
                rm( chipByChr, inputByChr )
                rm( chrInput )
                gc()
                
                # processing (using parallel computing, if multicore exists)
                
                if ( parallel == TRUE ) {
                    # if "multicore" package exists, utilize parallel computing with "mclapply"
                    #require(multicore)
                    
                    out <- mclapply( dataList, function(x) { 
                        .processBin_X( chip=x$chip, input=x$input, dataType=dataType ) 
                        }, mc.cores = nCore )
                } else {
                    # otherwise, use usual "lapply"
                    
                    out <- lapply( dataList, function(x) { 
                        .processBin_X( chip=x$chip, input=x$input, dataType=dataType ) 
                        } )
                }    
                
                rm( dataList )
                gc()
                    
                # stack processed data
                
                chrID <- coord <- Y <- X <- c()
                
                for ( i in 1:length(chrCommon) ) {
                    chrID <- c( chrID, rep( chrCommon[i], length(out[[i]]$coord) ) )
                    coord <- c( coord, out[[i]]$coord )
                    Y <- c( Y, out[[i]]$Y )
                    X <- c( X, out[[i]]$X )
                        
                    coordMat[i,3] <- min(out[[i]]$coord)
                    coordMat[i,4] <- max(out[[i]]$coord)
                }   
                
                rm( out )
                gc()               

                outBin <- new( "BinData", chrID=chrID, coord=coord, tagCount=Y, input=X, 
                    dataType=dataType, seqDepth=seqDepth )    
                
                rm( chrID, coord, Y, X )
                gc()
            } else {
                message( "Info: data contains only one chromosome." )
                
                out <- .processBin_X( chip=chip[,2:3], input=input[,2:3], dataType=dataType ) 
                
                chrID <- rep( chip[1,1], length(out$coord) )      
                outBin <- new( "BinData", chrID=chrID, coord=out$coord, tagCount=out$Y, 
                    input=out$X, dataType=dataType, seqDepth=seqDepth )  
                
                rm( chrID, out )
                gc() 
            }         
        }
    } else
    {
        ####################################################################
        #                                                                  #
        #                      one-sample analysis                         #
        #                                                                  #
        ####################################################################
        
        if ( existM & existGC & existN ) {  
                        
            if ( length(chrChIP) > 1 ) {
                # genome-wide analysis
                
                message( "Info: data contains more than one chromosome." )
                
                if ( ncol(mapScore) == 3 & ncol(gcScore) == 3 & ncol(nNuc) == 3 ) {                
                    # split data by chromosome
                    
                    chipByChr <- split( chip[,2:3], chip[,1] )
                    mapByChr <- split( mapScore[,2:3], mapScore[,1] )
                    gcByChr <- split( gcScore[,2:3], gcScore[,1] )
                    nByChr <- split( nNuc[,2:3], nNuc[,1] )
                    
                    # extract list of shared chromosomes
                    
                    chrM <- unique(mapScore[,1])
                    chrGC <- unique(gcScore[,1])
                    chrN <- unique(nNuc[,1])
                    
                    chrCommon <- Reduce( intersect,
                        list( chrChIP, chrM, chrGC, chrN ) )
                    chrCommon <- sort(chrCommon)
                    
                    # provide error & stop if there is no common chromosome
                    
                    if ( length(chrCommon) == 0 ) {
                        stop( "No chromosome is common among ChIP, mappability, GC, and N files." )
                    }
                    
                    # construct data to process
                    
                    coordMat <- matrix( NA, length(chrCommon), 4 )
                    
                    dataList <- list()
                    for ( i in 1:length(chrCommon) ) {
                        dataList[[i]] <- list()
                        dataList[[i]]$chip <- chipByChr[[ chrCommon[i] ]]
                        dataList[[i]]$mapScore <- mapByChr[[ chrCommon[i] ]]
                        dataList[[i]]$gcScore <- gcByChr[[ chrCommon[i] ]]
                        dataList[[i]]$nNuc <- nByChr[[ chrCommon[i] ]]
                        
                        coordMat[i,1] <- min(chipByChr[[ chrCommon[i] ]][,1])
                        coordMat[i,2] <- max(chipByChr[[ chrCommon[i] ]][,1])
                    }
                    
                    rm( chipByChr, mapByChr, gcByChr, nByChr )
                    rm( chrM, chrGC, chrN )
                    gc()
                    
                    # processing (using parallel computing, if multicore exists)
                    
                    if ( parallel == TRUE ) {
                        # if "multicore" package exists, utilize parallel computing with "mclapply"
                        #require(multicore)
                        
                        out <- mclapply( dataList, function(x) { 
                            .processBin_MGC( chip=x$chip,
                                mapScore=x$mapScore, gcScore=x$gcScore, nNuc=x$nNuc, 
                                dataType=dataType, rounding=rounding )
                            }, mc.cores = nCore )
                    } else {
                        # otherwise, use usual "lapply"
                        
                        out <- lapply( dataList, function(x) { 
                            .processBin_MGC( chip=x$chip,
                                mapScore=x$mapScore, gcScore=x$gcScore, nNuc=x$nNuc, 
                                dataType=dataType, rounding=rounding )
                            } )
                    }    
                    
                    rm( dataList )
                    gc()
                    
                    # stack processed data
                    
                    chrID <- coord <- Y <- M <- GC <- c()
                    
                    for ( i in 1:length(chrCommon) ) {
                        chrID <- c( chrID, rep( chrCommon[i], length(out[[i]]$coord) ) )
                        coord <- c( coord, out[[i]]$coord )
                        Y <- c( Y, out[[i]]$Y )
                        M <- c( M, out[[i]]$M )
                        GC <- c( GC, out[[i]]$GC )
                        
                        coordMat[i,3] <- min(out[[i]]$coord)
                        coordMat[i,4] <- max(out[[i]]$coord)
                    }         
                    
                    rm( out )
                    gc()
    
                    outBin <- new( "BinData", chrID=chrID, coord=coord, tagCount=Y, 
                        mappability=M, gcContent=GC, dataType=dataType, seqDepth=seqDepth )              
                    
                    rm( chrID, coord, Y, M, GC )
                    gc()
                } else {
                    stop( "Numbers of columns do not match!" )
                }
            } else {
                message( "Info: data contains only one chromosome." )
                
                if ( ncol(mapScore) == 2 & ncol(gcScore) == 2 & ncol(nNuc) == 2 ) {
                    # processing
                    
                    out <- .processBin_MGC( chip=chip[,2:3],
                        mapScore=mapScore, gcScore=gcScore, nNuc=nNuc, 
                        dataType=dataType, rounding=rounding )
                        
                    chrID <- rep( chip[1,1], length(out$coord) )    
                    outBin <- new( "BinData", chrID=chrID, coord=out$coord, tagCount=out$Y, 
                        mappability=out$M, gcContent=out$GC, dataType=dataType, seqDepth=seqDepth )   
                    
                    rm( chrID, out )
                    gc()
                } else {
                    stop( "Numbers of columns do not match!" )
                }   
            }   
        } else {
                stop( "All of mappability, GC content, and sequence ambiguity should be provided for one-sample analysis!" )     
        }
    }
    
    message( "Info: done!\n" )
    
    # info about preprocessing (for single chromosome data)
    
    if( length(chrChIP) == 1 & existN ) {
        nRatio <- length(which(nNuc[,2]==1)) / length(nNuc[,2])
        cat( "------------------------------------------------------------\n" )
        cat( "Info: preprocessing summary\n" )
        cat( "------------------------------------------------------------\n" )
        cat( "- percentage of bins with ambiguous sequences: ",round(nRatio*100),"%\n", sep="" )
        cat( "  (these bins will be excluded from the analysis)\n" )
        cat( "- before preprocessing:\n" )
        cat( "\tfirst coordinates = ",min(chip[,2]),
            ", last coordinates = ",max(chip[,2]), "\n", sep="" )
        cat( "- after preprocessing:\n" )
        cat( "\tfirst coordinates = ",min(outBin@coord),
            ", last coordinates = ",max(outBin@coord), "\n", sep="" )
        cat( "------------------------------------------------------------\n" )
    }
    
    # info about preprocessing (for multiple chromosome data)
    
    if( length(chrCommon) > 1 & existN ) {
        cat( "------------------------------------------------------------\n" )
        cat( "Info: preprocessing summary\n" )
        cat( "------------------------------------------------------------\n" )
        cat( "[Note] Bins with ambiguous sequences will be excluded from the analysis.\n" )
        cat( "Coordinates before & after preprocessing:\n" )
        for ( i in 1:nrow(coordMat) ) {
            cat( chrCommon[i],": ",sep="" )
        cat( "\t",coordMat[i,1]," - ",coordMat[i,2], "\t->", sep="" )
        cat( "\t",coordMat[i,3]," - ",coordMat[i,4], "\n", sep="" )
        }
        cat( "------------------------------------------------------------\n" )
    }
    
    return(outBin)
}

.validType <- function( type )
{
    # error treatment: check invalid type
    
    allType <- c("chip","M","GC","N","input")
    invalidType <- TRUE
    for ( i in 1:length(type) )
    {
        if ( length(which(!is.na(match(type[i],allType))))==0 )
        {
            invalidType <- FALSE
        }
    }
    if ( !invalidType )
    {
        stop( "Invalid 'type'! Choose among 'chip','M','GC','N','input'!" )
    }
    
    # error treatment: check whether type contains at least chip+M+GC+N or chip+input
    
    minType1 <- c("chip","M","GC","N")
    minType2 <- c("chip","input")
    
    notMin1 <- ( length(which(!is.na(match(minType1,type)))) < length(minType1) )
        # minimum requirement for one-sample analysis is not satisfied
    notMin2 <- ( length(which(!is.na(match(minType2,type)))) < length(minType2) )
        # minimum requirement for two-sample analysis (Input only) is not satisfied
    
    if ( notMin1 && notMin2 )
    {
        stop( "Minimum requirements of 'type':\nplease provide at least either c('chip','M','GC','N') or c('chip','input')!" )
    }
}

.validLocation <- function( type, fileName )
{    
    # error treatment: check whether 'fileName' info is provided
    
    if ( is.null(fileName) )
    {
        stop( "Please provide 'fileName' information!" )
    }
    
    # error treatment: check whether 'fileName' info matches 'type'    
    
    if ( length(type)!=length(fileName) )
    {
        #print(type)
        #print(fileName)
        stop( "Length of 'type' & length of 'fileName' do not match!" )
    }    
}
 
.validDataType <- function( dataType )
{   
    # error treatment: check invalid dataType
    
    if ( dataType!='unique' & dataType!='multi' )
    {
        stop( "Invalid 'dataType'! Choose either 'unique' or 'multi'!" )
    }
}

.processBin_MGC <- function( chip, mapScore, gcScore, nNuc, dataType, rounding )
{        
    # choose data with min # row
    
    minID <- which.min( c( nrow(chip), nrow(mapScore), nrow(gcScore) ) )
    dataList <- list( chip, mapScore, gcScore )
    dataMinID <- dataList[[ minID ]]
    rm( dataList, minID )
    
    # match coordinates
    
    chip <- chip[ !is.na(match(chip[,1],dataMinID[,1])), ]
    mapScore <- mapScore[ !is.na(match(mapScore[,1],dataMinID[,1])), ]
    gcScore <- gcScore[ !is.na(match(gcScore[,1],dataMinID[,1])), ]
    nNuc <- nNuc[ !is.na(match(nNuc[,1],dataMinID[,1])), ]
    
    # exclude bins with ambiguous sequences
    
    nRegID <- which(nNuc[,2]==1)
    if ( length(nRegID)>0 ) {
        chip <- chip[-nRegID,]
        mapScore <- mapScore[-nRegID,]
        gcScore <- gcScore[-nRegID,]
        nNuc <- nNuc[-nRegID,]
    }
    coord <- chip[,1]
    
    # round values
    
    Y <- round(chip[,2])    # round tag count in the case of multi match
    M <- round(mapScore[,2]*rounding)/rounding
    denom <- 1 - nNuc[,2]
    GC <- round(gcScore[,2]*rounding/denom)/rounding    
    
    if( dataType == 'unique' )
    {
        # in case of unique match
        if( length(which(Y>0&M==0)) > 0 )
        {
            Y[ Y>0 & M==0 ] <- 0    
        }
    }
    
    return( list( coord=coord, Y=Y, M=M, GC=GC ) )
}

.processBin_X <- function( chip, input, dataType )
{            
    # choose data with min # row
    
    minID <- which.min( c( nrow(chip), nrow(input) ) )
    dataList <- list( chip, input )
    dataMinID <- dataList[[ minID ]]
    rm( dataList, minID )
    
    # match coordinates
    
    chip <- chip[ !is.na(match(chip[,1],dataMinID[,1])), ]
    input <- input[ !is.na(match(input[,1],dataMinID[,1])), ]    
    coord <- chip[,1]
    
    # round values
    
    Y <- round(chip[,2])    # round tag count in the case of multi match
    X <- round(input[,2])
    
    return( list( coord=coord, Y=Y, X=X ) )
}

.processBin_XN <- function( chip, input, nNuc, dataType )
{            
    # choose data with min # row
    
    minID <- which.min( c( nrow(chip), nrow(input) ) )
    dataList <- list( chip, input )
    dataMinID <- dataList[[ minID ]]
    rm( dataList, minID )
    
    # match coordinates
    
    chip <- chip[ !is.na(match(chip[,1],dataMinID[,1])), ]
    input <- input[ !is.na(match(input[,1],dataMinID[,1])), ]
    nNuc <- nNuc[ !is.na(match(nNuc[,1],dataMinID[,1])), ]
    
    # exclude bins with ambiguous sequences
    
    nRegID <- which(nNuc[,2]==1)
    if ( length(nRegID)>0 ) {
        chip <- chip[-nRegID,]
        input <- input[-nRegID,]
        nNuc <- nNuc[-nRegID,]
    }
    coord <- chip[,1]
    
    # round values
    
    Y <- round(chip[,2])    # round tag count in the case of multi match
    X <- round(input[,2])     
    
    return( list( coord=coord, Y=Y, X=X ) )
}

.processBin_MGCX <- function( chip, input, mapScore, gcScore, nNuc, dataType, rounding )
{
    # choose data with min # row
    
    minID <- which.min( c( nrow(chip), nrow(input), nrow(mapScore), nrow(gcScore) ) )
    dataList <- list( chip, input, mapScore, gcScore )
    dataMinID <- dataList[[ minID ]]
    rm( dataList, minID )
    
    # match coordinates
    
    chip <- chip[ !is.na(match(chip[,1],dataMinID[,1])), ]
    input <- input[ !is.na(match(input[,1],dataMinID[,1])), ]
    mapScore <- mapScore[ !is.na(match(mapScore[,1],dataMinID[,1])), ]
    gcScore <- gcScore[ !is.na(match(gcScore[,1],dataMinID[,1])), ]
    nNuc <- nNuc[ !is.na(match(nNuc[,1],dataMinID[,1])), ]
    
    # exclude bins with ambiguous sequences
    
    nRegID <- which(nNuc[,2]==1)
    if ( length(nRegID)>0 ) {
        chip <- chip[-nRegID,]
        mapScore <- mapScore[-nRegID,]
        gcScore <- gcScore[-nRegID,]
        nNuc <- nNuc[-nRegID,]
        input <- input[-nRegID,]
    }
    coord <- chip[,1]
    
    # round values
    
    Y <- round(chip[,2])    # round tag count in the case of multi match
    M <- round(mapScore[,2]*rounding)/rounding
    denom <- 1 - nNuc[,2]
    GC <- round(gcScore[,2]*rounding/denom)/rounding
    X <- round(input[,2])    
    
    if( dataType == 'unique' )
    {
        # in case of unique match
        if( length(which(Y>0&M==0)) > 0 )
        {
            Y[ Y>0 & M==0 ] <- 0    
        }
    }    
    
    return( list( coord=coord, Y=Y, X=X, M=M, GC=GC ) )
}
