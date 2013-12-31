
# calculate log likelihood for MOSAiCS

.logLik <- function( mosaicsEst, tagCount, pNfit, k=3, signalModel="2S" )
{              
    switch( signalModel,
        "1S" = {
            # one signal model
            
            pi0 <- mosaicsEst@pi0            
            fitMD <- .margDist_1S( mosaicsEst=mosaicsEst, 
                tagCount=tagCount, pNfit=pNfit, k=k )
            
            loglik <- sum( log( pi0 * fitMD$MDZ0 + ( 1 - pi0 ) * fitMD$MDZ1 ) )
        },
        "2S" = {
            # two signal model
            
            pi0 <- mosaicsEst@pi0
            p1 <- mosaicsEst@p1            
            fitMD <- .margDist_2S( mosaicsEst=mosaicsEst, 
                tagCount=tagCount, pNfit=pNfit, k=k )
                
            loglik <- sum( log( pi0 * fitMD$MDZ0 + 
            	( 1 - pi0 ) * ( fitMD$MDZ1 * p1 + fitMD$MDZ2 * (1-p1) ) ) )
        }
    )
    
    return(loglik)
}
