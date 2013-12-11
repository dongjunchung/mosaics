
.calcModelBIC <- function( fitZ1, Y, pNfit, k=3, model="2S", type="BIC", npar )
{       
    # parameter estimates (common)
    
    mu_est <- fitZ1$muEst
    a <- fitZ1$a
    pi0 <- fitZ1$pi0
    b_est <- a / mu_est        
    PYZ0 <- pNfit$PYZ0
    
    # choose Y>=k
    
    id_geqk <- which(Y>=k)      
    Y_ori <- Y[id_geqk]
    b_est <- b_est[id_geqk]
    mu_est <- mu_est[id_geqk]
    
    Yk <- Y_ori - k        # use only Y >= k
    #if(length(which(Y<0))>0 ) Y[which(Y<0)] <- -1 
    Ykmax <- max(Yk)
    #ind_ge_k <- which(Y>=0)
    
    if ( model=="2S" )
    {
        # parameter estimates
        
        p1 <- fitZ1$p1
        b1 <- fitZ1$b1
        c1 <- fitZ1$c1
        b2 <- fitZ1$b2
        c2 <- fitZ1$c2
        
        # calculate log likelihood
        
        #PYZ1 <- .margDistZ1_2S( Y_ori, pNfit, b1, c1, b2, c2 )
        PYZ1 <- .margDistZ1_2S( Yk, Ykmax, pNfit, b1, c1, b2, c2 )
        PYZ1G1 <- PYZ1$MDG1
        PYZ1G2 <- PYZ1$MDG2    
        
        logLik0 <- log( pi0*PYZ0 + (1-pi0)*( p1*PYZ1G1 + (1-p1)*PYZ1G2) )
        logLik1 <- sum( logLik0[!is.na(logLik0)] )
    } else if ( model=="1S" )
    {
        # parameter estimates
        
        b <- fitZ1$b
        c <- fitZ1$c
        
        # calculate log likelihood
        
        #PYZ1 <- .margDistZ1_1S( Y_ori, pNfit, b, c )
        PYZ1 <- .margDistZ1_1S( Yk, Ykmax, pNfit, b, c )
            
        logLik0 <- log( pi0*PYZ0 + (1-pi0)*PYZ1 )
        logLik1 <- sum( logLik0[!is.na(logLik0)] )    
    }
    
    # calculate BIC
    
    switch( type,
        "AIC" = {
            penalty <- 2            
        },
        "BIC" = {
            n <- length(Y_ori)
            penalty <- log(n)
        }
    )
    val <- -2 * logLik1 + penalty * npar
    
    return(val)
}
