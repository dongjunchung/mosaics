#########################################################
# Compute marginal density for Y=N+S1+S2
#########################################################

.margDist_2S <- function( mosaicsEst, tagCount, pNfit, k=3 )
{
    # extract parameters
    
    a <- mosaicsEst@a    
    mu_est <- mosaicsEst@muEst
    b_est <- a / mu_est
        
    b1 <- mosaicsEst@b1
    c1 <- mosaicsEst@c1
    b2 <- mosaicsEst@b2
    c2 <- mosaicsEst@c2
    
    Yori <- tagCount
    
    
    # use only Y >= k
    
    Y = Yori - k        
    Y[which(Y<0)] = -1  
     
    
    # round mu to the nearest hundredth

    mu_round <- round(mu_est,2)
    if( length(which(Y<0))>0 ) mu_round[which(Y<0)] <- 0
    mu_round_U <- unique(mu_round)
    #n_mu_U <- length(mu_round_U)
    
    
    # prob of N using rounding mu for prob of S1 & S2
     
    Ymax <- max(Y)
    pN <- pNfit$pN
    
    
    # prob of N (Yori & b_est have same length)
    
    MDZ0 = dnbinom( Yori, size=a, b_est/(b_est+1) )
    
    
    # prob of S1 & S2 (positive only when Y>=k)
        
    pS1 = dnbinom( 0:Ymax, b1, c1/(c1+1) )
    pS2 = dnbinom( 0:Ymax, b2, c2/(c2+1) )
    
    id_Y_ge_k <- which(Y>=0)    
    MDG <- matrix( 0, length(Y), 2 )
    MDGfit <- conv_2S( y=Y[id_Y_ge_k], mu_round=mu_round[id_Y_ge_k],
        mu_round_U=mu_round_U, pN=pN, pS1=pS1, pS2=pS2 )
    MDGfit <- matrix( MDGfit, nrow=2 )  
    MDG[ id_Y_ge_k, ] <- t(MDGfit)
    MDZ1 <- MDG[,1]
    MDZ2 <- MDG[,2]
    
    
    # return object

    MD=data.frame( cbind( MDZ0, MDZ1, MDZ2 ) )
    colnames(MD)=c('MDZ0','MDZ1','MDZ2')
    return(MD)
}
