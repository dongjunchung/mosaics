
#########################################################
# Compute mean/variance for sequenced samples for 
# exploratory analysis
#########################################################

.computeStat <- function( Y, S )
{    
    # construct unique M/GC pair
    
    tmp <- cbind( Y, S )
    tempdf <- within( as.data.frame(tmp), {S <- factor(S)} )    
    ySubList <- with( tempdf, split(Y, factor(S)) )  
    
    # calculate strata-specific parameters
    
    yStatList <- lapply( ySubList, function(Y) {
            Y <- as.numeric( as.vector(Y) )
            
            nitem <- length(Y)
            nitemgeq0 <- length( Y[Y!=0] )
            p0 <- length(which( Y==0 )) / length(Y)
            p1 <- length(which( Y==1 )) / length(Y)
            meanYgeq0 <- mean( Y[Y!=0] )
            varYgeq0 <- var( Y[Y!=0] )
            medYgeq0 <- median( Y[Y!=0] )
            meanYall <- mean(Y)
            varYall <- var(Y)
            
            result <- c( nitem, nitemgeq0, p0, p1,
                meanYgeq0, varYgeq0, medYgeq0, meanYall, varYall )
            
            return( result )
        }
    )
    yStatMat <- matrix( unlist(yStatList), ncol = 9, byrow = TRUE )
    
    # construct summary
    
    yStats <- list( uS = as.numeric(as.vector(names(ySubList))),
        nitem = as.numeric(as.vector(yStatMat[, 1])),
        nitemgeq0 = as.numeric(as.vector(yStatMat[, 2])),
        p0 = as.numeric(as.vector(yStatMat[, 3])),
        p1 = as.numeric(as.vector(yStatMat[, 4])), 
        meanYgeq0 = as.numeric(as.vector(yStatMat[, 5])), 
        varYgeq0 = as.numeric(as.vector(yStatMat[, 6])), 
        medYgeq0 = as.numeric(as.vector(yStatMat[, 7])),
        meanYall = as.numeric(as.vector(yStatMat[, 8])),
        varYall = as.numeric(as.vector(yStatMat[, 9]))
    )
        
    return(yStats)
}
