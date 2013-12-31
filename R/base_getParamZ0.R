
.getParamZ0 <- function( Y, bgEst=NA, Y_freq ) {      
    
    # encode ty to return numeric vector
    # Q -> 1, MM -> 2, NA -> 3
    
    u0 <- length(which(Y==0))
    u1 <- length(which(Y==1))
    u2 <- length(which(Y==2))

    #print( (u0+u1+u2)/length(Y) )
    
    # determine background estimation method
    
    if ( u0>0 & u1>0 & u2>0 ) {                
        # if we have all of 0, 1, 2, then consider "bgEst"
        
        if ( bgEst=="matchLow" ) {
            estMethod <- "matchLow"
        } else if ( bgEst=="rMOM" ) {
            estMethod <- "rMOM"
        } else {
            # the case that "bgEst" is not specified:
            # if "bgEst" not specified & bins w/ 0, 1, 2 is more than 50%, then use "low"
            # otherwise, use "high"
            
            #if ( sum(Y_freq[ as.numeric(names(Y_freq))<=2 ]) / sum(Y_freq) > 0.5 ) {
            #    estMethod <- "matchLow"
            #} else {
            #    estMethod <- "rMOM"
            #}
        }
    } else {
        # if more than one of 0, 1, 2 are missing, then use only "high"
        
        estMethod <- "rMOM"
    }
    
    # estimate parameters
    
    if ( estMethod == "matchLow" ) {  
          
        # matching 0, 1, 2 counts

        #print("matching 0, 1, 2 counts")
        
        r1 <- u1/u0
        r2 <- u2/u1
        a <- r1/(2*r2-r1)
        b <- 1/(2*r2-r1)-1
        #ty <- 'Q'
        ty <- 1
    
        mean0 <- mean(Y)
        var0 <- var(Y)
        
        if( a<=0 || b<=0 || a==Inf || b==Inf ) {
            # MOM, if a & b estimates based non 0, 1, 2 counts are improper 
            # - use all observations, assuming that pi0>=1.
            
            if( !is.na(var0) && !is.na(mean0) && var0 > mean0 ) {
                # use MOM only when we have proper mean and variance
                
                a <- mean0^2/(var0-mean0)
                b <- mean0/(var0-mean0)
                #ty <- 'MM'
                ty <- 2
            } else {
                a <- NA
                b <- NA
                #ty <- NA
                ty <- 3
            }
        }
    } else if ( estMethod == "rMOM" ) {
        # MOM (mosaics, ver 1.0.2): using robust statistics 
        # (mosaics, ver 1.0.6): median and MAD for estimates of mean and variance
    
        #print( "rMOM" )
        
        mean0 <- median(Y)
        var0 <- ( mad(Y) )^2
        
        if( !is.na(var0) && !is.na(mean0) &&  var0 > mean0 ) {
            # use MOM only when we have proper mean and variance
            
            a <- mean0^2 / (var0-mean0)
            b <- mean0 / (var0-mean0)
            #ty <- 'MM'
            ty <- 2
        } else {
            a <- NA
            b <- NA
            #ty <- NA
            ty <- 3
        }
    }
    
    return( c( a, b, mean0, var0, u0, u1, u2, length(Y), ty ) )
}
