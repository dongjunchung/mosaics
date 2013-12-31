.calcYbdAll <- function( MOSAiCS_Z0, k ) {
    
    a <- MOSAiCS_Z0$a
    mu_est <- MOSAiCS_Z0$muEst
    b_est <- a / mu_est
    pi0 <- MOSAiCS_Z0$pi0
    tab_Y <- MOSAiCS_Z0$Y_freq
    Y_u <- MOSAiCS_Z0$Y_val
        
    b_estt <- table(b_est)
    b_estv <- as.numeric(names(b_estt))
    b_estn <- as.numeric(b_estt)
    
    Eyz0 <- rep(0,length(Y_u))                  # N_Z0(y)
    #Eyz0[1] <- round(sum(dnbinom(Y_u[1],a,b_est/(b_est+1))))
    Eyz0[1] <- round(sum( b_estn * dnbinom(Y_u[1],a,b_estv/(b_estv+1)) ))

    i <- 2
    while( Eyz0[i-1]>0 & i<=length(Y_u) )
    {
        Eyz0[i] <- round(sum( b_estn * dnbinom(Y_u[i],a,b_estv/(b_estv+1)) ))
        i <- i + 1      
    }

    Y_Z1_tmp <- tab_Y - pi0*Eyz0                
        # (1-pi0) * N_Z1(y) = N(y) - pi0 * N_Z0(y)
    Y_Z1_tmp[which(Y_Z1_tmp<0)] <- 0
    Y_bd <- cbind( Y_u[-c(1:3)]-k, Y_Z1_tmp[-c(1:3)] ) 
        ### already adjusted for k
    
    if ( sum(Y_bd[,2])==0 ) {
        # nothing left after background are excluded from the data
        # risk level = extremely high
        stop( "over-estimation of background detected. Please tune the parameters!" )
    } else if ( Y_bd[which.max(Y_bd[,2]),2] / sum(Y_bd[,2]) >= 0.90 ) {
        # almost nothing left after background are excluded from the data
        # risk level = very high
        stop( "over-estimation of background detected. Please tune the parameters!" )
    } else if ( Y_bd[which.max(Y_bd[,2]),2] / sum(Y_bd[,2]) >= 0.80 ) {
        # still not much left after background are excluded from the data
        # rick level = medium high
        warning( "there is some possibility of background over-estimation.\n" )
        warning( "parameter tuning might provide better fits.\n" )
        Y_bd_all <- rep(Y_bd[,1],Y_bd[,2])
    } else {
        Y_bd_all <- rep(Y_bd[,1],Y_bd[,2])
    }    
    
    return( Y_bd_all )
}
