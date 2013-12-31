
#.mosaicsZ0 <- function( Y, analysisType=NA, X=NA, M=NA, GC=NA, truncProb=0.9999, Y_freq )
.mosaicsZ0 <- function( Y, bgEst="automatic", analysisType=NA, 
    X=NA, M=NA, GC=NA, inputTrunc=NA, Y_freq,
    parallel=parallel, nCore=nCore )
{    
    #library(MASS)    
    
    # truncate X for input only analysis
    
    if ( analysisType=="IO" ) {
        #trunc <- quantile( X, prob=truncProb )
        #X[which(X>trunc)] <- trunc
        X[which(X>inputTrunc)] <- inputTrunc
    }
    
    # construct unique M/GC pair
    
    switch( analysisType,
        OS = { S1 <- paste( M, GC, sep = ':' ) },
        TS = { S1 <- paste( X, M, GC, sep = ':' ) },
        IO = { S1 <- X }
    )
    #tmp <- cbind( Y, S1 )
    #tempdf <- within( as.data.frame(tmp), {S1 <- factor(S1)} )    
    tmp <- data.frame( Y, as.factor(S1) )
    ySubList <- with( tmp, split(Y, S1) )    
    #ySubList <- split( Y, S1 )
    
    # construct strata index
    
    if ( analysisType=="IO" ) {
        indStrata <- as.numeric(as.vector(names(ySubList)))
    } else {
        indStrata <- apply( as.matrix(names(ySubList)), 1, 
            function(x){as.numeric(strsplit(x, ":")[[1]])} )
        indStrata <- t(indStrata)    
    }
    
    # calculate strata-specific parameters
    
    if ( parallel ) {
        yParamList <- mclapply( ySubList, 
            function(x) .getParamZ0( x, bgEst=bgEst, Y_freq=Y_freq ),
            mc.cores=nCore )    
    } else {
        yParamList <- lapply( ySubList, 
            function(x) .getParamZ0( x, bgEst=bgEst, Y_freq=Y_freq ) )
    }
    yParamMat <- matrix( unlist(yParamList), ncol = 9, byrow = TRUE )
    
    # decode ty
    
    ty <- rep( "", nrow(yParamMat) )
    ty[ yParamMat[, 9]==1 ] <- "Q"
    ty[ yParamMat[, 9]==2 ] <- "MM"
    ty[ yParamMat[, 9]==3 ] <- NA
    
    # construct common summary
    
    parEstZ0 <- list( 
        a_u = yParamMat[, 1], b_u = yParamMat[, 2],
        mean0_u = yParamMat[, 3], var0_u = yParamMat[, 4], 
        u0_u = yParamMat[, 5], u1_u = yParamMat[, 6], u2_u = yParamMat[, 7],
        n_u = yParamMat[, 8], ty_u = ty,
        Y_val = as.numeric(names(table(Y))), Y_freq = table(Y), Y_sub_list=ySubList )
    
    # constuct analysis type specific summary
    
    switch( analysisType,
        OS = {
            parEstZ0$M_u <- indStrata[,1,drop=FALSE]
            parEstZ0$GC_u <- indStrata[,2,drop=FALSE]
        },
        TS = {
            parEstZ0$X_u <- indStrata[,1,drop=FALSE]
            parEstZ0$M_u <- indStrata[,2,drop=FALSE]
            parEstZ0$GC_u <- indStrata[,3,drop=FALSE]       
        },
        IO = {
            parEstZ0$X_u <- indStrata
            #parEstZ0$trunc <- trunc      
        }
    )
        
    return(parEstZ0)
}
