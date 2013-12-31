#########################################################
# Compute marginal density for Y=N+S1+S2
#########################################################

.margDist_1S <- function( mosaicsEst, tagCount, pNfit, k=3 )
{            
    # extract parameters
    
    a <- mosaicsEst@a
    mu_est <- mosaicsEst@muEst
    b_est <- a / mu_est
        
    b <- mosaicsEst@b
    c <- mosaicsEst@c
    
    Yori <- tagCount
    
    
    # use only Y >= k
    
    Y <- Yori - k 
    Y[which(Y<0)] = -1  
    
    
    # round mu

    mu_round <- round(mu_est,2)
    if( length(which(Y<0))>0 ) mu_round[which(Y<0)] <- 0
    mu_round_U <- unique(mu_round)
    #n_mu_U <- length(mu_round_U)
    
    
    # prob of N using rounding mu for prob of S
     
    Ymax <- max(Y)
    pN <- pNfit$pN
    
    
    # prob of N (Yori & b_est have same length)
    
    MDZ0 = dnbinom( Yori, size=a, b_est/(b_est+1) )
    
    
    # prob of S (positive only when Y>=k)
    
    pS <- dnbinom( 0:Ymax, b, c/(c+1) )
   
    id_Y_ge_k <- which( Y>=0 )
    #n_Y_ge_k <- length( id_Y_ge_k )
    MDZ1 <- rep( 0, length(Y) )
    MDZ1[ id_Y_ge_k ] <- conv_1S( y=Y[id_Y_ge_k], mu_round=mu_round[id_Y_ge_k],
        mu_round_U=mu_round_U, pN=pN, pS=pS )
            
    
    # return object

    MD=data.frame( cbind( MDZ0, MDZ1 ) )
    colnames(MD)=c('MDZ0','MDZ1')
    return(MD)   
}
