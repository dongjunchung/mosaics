
# calculate posterior probabilities (1 signal model)

.getPH_1S <- function( margDensity, pi0 )
{
    # take parameters
    
    MDZ0 <- margDensity$MDZ0
    MDZ1 <- margDensity$MDZ1
    
    
    # calculate posterior probabilities
    
    denom <- MDZ0*pi0 + MDZ1*(1-pi0)
    pH0 = MDZ0 * pi0 / denom  
    pH1 = MDZ1 * (1-pi0) / denom

    if( length(which(is.na(pH0))) > 0 )
    {
        pH0[ is.na(pH0) ] = 1
        pH1[ is.na(pH0) ] = 0
    }
    post_prob = list( pH0=pH0, pH1=pH1 )
    
    return( post_prob )
}


# calculate posterior probabilities (2 signal model)

.getPH_2S <- function( margDensity, pi0, p1 )
{
    # take parameters
    
    MDZ0 <- margDensity$MDZ0
    MDZ1 <- margDensity$MDZ1
    MDZ2 <- margDensity$MDZ2
    
    
    # calculate posterior probabilities
    
    denom <- MDZ0*pi0 + (1-pi0) * ( MDZ1*p1 + MDZ2*(1-p1) )
    pH0 = MDZ0 * pi0 / denom  
    pH1 = MDZ1 * (1-pi0) * p1 / denom
    pH2 = MDZ2 * (1-pi0) * (1-p1) / denom
    

    if( length(which(is.na(pH0))) > 0 )
    {
        pH0[ is.na(pH0) ] = 1
        pH1[ is.na(pH0) ] = 0
        pH2[ is.na(pH0) ] = 0
    }
    post_prob = list( pH0=pH0, pH1=pH1, pH2=pH2 )
    
    return( post_prob )
}


# claim peaks

.peakCall <- function( postProb, dataSet, FDR, binsize=NA, 
    maxgap=200, minsize=50, thres=10, analysisType, nRatio=1 )
{   
    #library(IRanges)
    
    # FDR = direct posterior probability approach (Newton et al., 2004, Biostatistics)
    
    # determine peaks

    betapH <- postProb$pH0
    betapH_s <- sort(betapH)
    sbetapH <- cumsum(betapH_s) / c(1:length(betapH))       # expected rate of false discoveries
    id <- which( sbetapH <= FDR )
    
    if(length(id)>0)
    {
        #################### if peaks exist
        
        chrID <- dataSet$chrID
        coord <- dataSet$coord
        Y <- dataSet$Y
        switch( analysisType,
            OS = {
                M <- dataSet$M
                GC <- dataSet$GC
            },
            TS = {
                X <- dataSet$X
                M <- dataSet$M
                GC <- dataSet$GC
            },
            IO = {
                X <- dataSet$X
            }
        ) 
                
        
        # threshold peaks by min tag count & determine peaks
        
        cutoff <- betapH_s[max(id)]    
        bd_bin <- rep( 0, length(Y) )
        bd_bin[ betapH<=cutoff & Y>thres ] <- 1
        
        if ( length(which( betapH<=cutoff & Y>thres )) > 0 ) {
            # if we still have peaks after tag count thresholding
             
            # empirical FDR
            
            empFDR <- sum(betapH[which(betapH<=cutoff)]) / length(which(betapH<=cutoff))
            empFDR_thres <- sum(betapH[which( betapH<=cutoff & Y>thres )]) / 
                length(which( betapH<=cutoff & Y>thres ))
            #cat( "Info: empirical FDR (before thresholding) = ", 
            #    round(1000*empFDR)/1000, "\n", sep="" )        
            #cat( "Info: empirical FDR (after thresholding) = ", 
            #    round(1000*empFDR_thres)/1000, "\n", sep="" )
            
            # binsize calculation
            
            if ( is.na(binsize) ) {
                binsize <- min(abs( diff(coord) ))
            }
            
            # process peaks for each chromosome
            
            chrList <- sort(unique(chrID))
            final_peakset <- c()
                
            bd_bin_list <- split( bd_bin, chrID )
            betapH_list <- split( betapH, chrID )
            coord_list <- split( coord, chrID )
            Y_list <- split( Y, chrID )
            switch( analysisType,
                OS = {
                    M_list <- split( M, chrID )
                    GC_list <- split( GC, chrID )
                },
                TS = {
                    X_list <- split( X, chrID )
                    M_list <- split( M, chrID )
                    GC_list <- split( GC, chrID )
                },
                IO = {
                    X_list <- split( X, chrID )
                }
            ) 
            
            for ( chr in 1:length(chrList) ) {
                # extract data for given chromosome
                
                bd_bin_chr <- bd_bin_list[[ chrList[chr] ]]
                betapH_chr <- betapH_list[[ chrList[chr] ]]
                coord_chr <- coord_list[[ chrList[chr] ]]
                Y_chr <- Y_list[[ chrList[chr] ]]
                switch( analysisType,
                    OS = {
                        X_chr <- NA
                        M_chr <- M_list[[ chrList[chr] ]]
                        GC_chr <- GC_list[[ chrList[chr] ]]
                    },
                    TS = {
                        X_chr <- X_list[[ chrList[chr] ]]
                        M_chr <- M_list[[ chrList[chr] ]]
                        GC_chr <- GC_list[[ chrList[chr] ]]
                    },
                    IO = {
                        X_chr <- X_list[[ chrList[chr] ]]
                        M_chr <- NA
                        GC_chr <- NA
                    }
                ) 
                                
                # generate initial peak list
          
                bd_ID <- which(bd_bin_chr==1)           # initial peak (bin-level)
                
                if ( length(bd_ID) > 0 ) {
                    # to take care of the case that there is no peak in this chromosome
                
                    #indRanges <- IRanges( start=bd_ID, end=bd_ID+1 )
                    #reducedRanges <- reduce(indRanges)  # merge nearby peaks
                    
                    #binsize <- coord_chr[2] - coord_chr[1]                
                    #peak_start <- coord_chr[ start(reducedRanges) ]
                    #peak_stop <- coord_chr[ end(reducedRanges)-1 ] + binsize - 1
                    #peak_start <- coord_chr[ start(indRanges) ]
                    #peak_stop <- coord_chr[ end(indRanges)-1 ] + binsize - 1
                    
                    peak_range <- IRanges( start=coord_chr[ bd_ID ], 
                        end=coord_chr[ bd_ID ] + binsize - 1 )
                    peak_reduced <- reduce(peak_range)
                    peak_start <- start(peak_reduced)
                    peak_stop <- end(peak_reduced)
                    
                    # merge close peaks if distance<=maxgap
                    
                    coord_org <- IRanges( start=peak_start, end=peak_stop+maxgap )
                    coord_merged <- reduce(coord_org)
                    peak_start <- start(coord_merged)
                    peak_stop <- end(coord_merged) - maxgap                    
                    
                    # filter peaks smaller than minsize & order by coordinates
                    
                    peaksize <- peak_stop - peak_start + 1
                    filterID = which( peaksize <= minsize )
                    if ( length(filterID) > 0 )
                    {
                        peak_start <- peak_start[ -filterID ]
                        peak_stop <- peak_stop[ -filterID ]
                    }
                    peak_start <- peak_start[ order(peak_start) ]        
                    peak_stop <- peak_stop[ order(peak_start) ]
                    peaksize <- peak_stop - peak_start + 1
                    
                    #print(quantile(peaksize))
                    
                    # calculate additional info
      
                    final_peakset_chr <- .annotatePeak( 
                      peakStart_chr=peak_start, peakStop_chr=peak_stop, 
                      coord_chr=coord_chr, analysisType=analysisType,
                      Y_chr=Y_chr, X_chr=X_chr, M_chr=M_chr, GC_chr=GC_chr, 
                      pp_chr=betapH_chr, nRatio=nRatio )
                    
                    final_peakset_chr <- data.frame( chrList[chr], final_peakset_chr )
                    colnames(final_peakset_chr)[1] <- "chrID"
                  
                    final_peakset <- rbind( final_peakset, final_peakset_chr )  
                }              
            }
            
            
            #################### if peaks exist
            
            bdBin <- data.frame( chrID, coord, bd_bin, 
            	stringsAsFactors=FALSE )
			colnames(bdBin) <- c( "chrID", "coord", "peak" )
            
            return( list( peakSet=final_peakset, bdBin=bdBin, empFDR=empFDR ) )
        } else {
            # if we lose all peaks after tag count thresholding
             
            # empirical FDR
            
            #empFDR <- sum(betapH[which(betapH<=cutoff)]) / length(which(betapH<=cutoff))
            #empFDR_thres <- 0
            #message( "Info: no peaks remain after thresholding")
            #message( "Info: empirical FDR (before thresholding) = ", round(1000*empFDR)/1000 )        
            #message( "Info: empirical FDR (after thresholding) = ", round(1000*empFDR_thres)/1000 )
            
            chrID <- dataSet$chrID
            coord <- dataSet$coord
            Y <- dataSet$Y
            
            bdBin <- data.frame( chrID, coord, rep(0,length(Y)), 
            	stringsAsFactors=FALSE )
        colnames(bdBin) <- c( "chrID", "coord", "peak" )
            
            return( list( peakSet=NULL, bdBin=bdBin, empFDR=0 ) )
        
        }
    } else
    {
        #################### if there is no peak
        
        chrID <- dataSet$chrID
            coord <- dataSet$coord
        Y <- dataSet$Y
            
        bdBin <- data.frame( chrID, coord, rep(0,length(Y)), 
        	stringsAsFactors=FALSE )
        colnames(bdBin) <- c( "chrID", "coord", "peak" )
        
        return( list( peakSet=NULL, bdBin=bdBin, empFDR=0 ) )  
    }
}
