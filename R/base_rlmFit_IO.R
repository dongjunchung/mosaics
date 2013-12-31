
.rlmFit_IO <- function( parEst, d=0.25, Y, X, inputTrunc )
{    
    #library(MASS)
    
    
    # exclude NA & MM points
    
    if ( length(which(is.na(parEst$a_u)==TRUE))==length(parEst$a_u) ) {
        stop( "insufficient # of proper strata! Cannot proceed!" )
    }
    
    #idNA <- unique( c( which(is.na(parEst$a_u)==TRUE), which(as.character(parEst$ty_u)=='MM') ) )
    idNA <- which( is.na(parEst$a_u) )
    if ( length(idNA)>0 ) {
        a_u <- parEst$a_u[-idNA]
        mu_u <- parEst$a_u[-idNA] / parEst$b_u[-idNA]
        X_u <- parEst$X_u[-idNA]
        n_u <- parEst$n_u[-idNA]
        mean0 <- parEst$mean0_u[-idNA]
        var0 <- parEst$var0_u[-idNA]
    } else {
        a_u <- parEst$a_u
        mu_u <- parEst$a_u / parEst$b_u
        X_u <- parEst$X_u
        n_u <- parEst$n_u
        mean0 <- parEst$mean0_u
        var0 <- parEst$var0_u    
    }
    #trunc <- parEst$trunc
    
    # estimators of a

    a_w_strata <- sum(n_u*a_u)/sum(n_u)         # [Note] No truncation for input only analysis

    
    # sample size weighted rlm fit    
    
    #X[ which(X>trunc) ] <- trunc                  # [Note] we need truncation of X for input only analysis
    X[ which(X>inputTrunc) ] <- inputTrunc
    
    Xtrans <- X_u^d    
    fit_trunc_adj <- rlm( log(mu_u) ~ Xtrans, weights=n_u/sum(n_u) )
        
        
    # return estimates
    
    coef_rlm <- coef(fit_trunc_adj)
    muEst <- exp( coef_rlm[1] + coef_rlm[2]*(X^d) )

    pi0 <- .getPi0( muEst, a_w_strata, parEst$Y_freq )
    
    betaEst <- coef(fit_trunc_adj)
    names(betaEst) <- c( "(intercept)", "input" )

    return( list( pi0 = pi0, a = a_w_strata, muEst=muEst, betaEst=betaEst,
        Y_val=parEst$Y_val, Y_freq=parEst$Y_freq ))
    
}
