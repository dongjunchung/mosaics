
.rlmFit_TS  <- function( parEst, mean_thres=1, s=2, d=0.25, bgEst, Y, M, GC, X )
{
    #library(MASS)
    #library(splines)
    
    
    # exclude NA & MM points
    
    if ( length(which(is.na(parEst$a_u)==TRUE))==length(parEst$a_u) ) {
        stop( "insufficient # of proper strata! Cannot proceed!" )
    }
    
    # do not excluded estimtes from MOM (ver 1.0.2)
    
    #idNA <- unique(c(which(is.na(parEst$a_u)==TRUE),which(as.character(parEst$ty_u)=='MM')))
    #idNA <- which( is.na(parEst$a_u) | as.character(parEst$ty_u)=='MM' )
    idNA <- which( is.na(parEst$a_u) )
    if ( length(idNA)>0 ) {
        a_u <- parEst$a_u[-idNA]
        mu_u <- parEst$a_u[-idNA] / parEst$b_u[-idNA]
        X_u <- parEst$X_u[-idNA]
        M_u <- parEst$M_u[-idNA]
        GC_u <- parEst$GC_u[-idNA]
        n_u <- parEst$n_u[-idNA]
        mean0 <- parEst$mean0_u[-idNA]
        var0 <- parEst$var0_u[-idNA]
    } else {
        a_u <- parEst$a_u
        mu_u <- parEst$a_u / parEst$b_u
        X_u <- parEst$X_u
        M_u <- parEst$M_u
        GC_u <- parEst$GC_u
        n_u <- parEst$n_u
        mean0 <- parEst$mean0_u
        var0 <- parEst$var0_u    
    }
    
    
    # exclude mu>Yhat points
    
    id_bigmean <- which( ( mu_u - mean0 ) > mean_thres )  # [Note_1] mean_thres=1
    if( length(id_bigmean)>0 ) {
        a_u <- a_u[-id_bigmean]
        mu_u <- mu_u[-id_bigmean]
        X_u <- X_u[-id_bigmean]
        M_u <- M_u[-id_bigmean]
        GC_u <- GC_u[-id_bigmean]
        n_u <- n_u[-id_bigmean]
        mean0 <- mean0[-id_bigmean]
        var0 <- var0[-id_bigmean]
    }
    
    
    # estimators of a

    a_wq_strata <- sum(n_u*a_u)/sum(n_u)

    #idNA <- unique(c(which(a_u>quantile(a_u,prob=0.95)),which(a_u<quantile(a_u,prob=0.05))))
    idNA <- which( a_u>quantile(a_u,prob=0.95) | a_u<quantile(a_u,prob=0.05) )
    if( length(idNA)>0 ) {
        a_u_wq <- a_u[-idNA]
        n_u_wq <- n_u[-idNA]
        a_wq_strata <- sum(n_u_wq*a_u_wq) / sum(n_u_wq) ### Note using a_wq_strata
    }
    
    
    # sample size weighted rlm fit
    
    GC_knots <- c(quantile(GC_u,prob=0.25),quantile(GC_u,prob=0.75))
    GC_bs <- bs(GC_u,degree=1,knots=GC_knots)           # [Note_2] different GC_knots
    
    IX <- rep(0,length(X_u))
    IX[ X_u <= s ] <- 1
    
    # exception handling: no 0, 1, 2 counts in input
    if ( length(which( IX==0 )) == length(IX) ) {
        stop( paste("insufficient # of bins with counts <=", s, "in control sample. Please try input-only analysis!") )
    }
    
    MX_u <- log2(M_u+1) * IX
    GC1X_u <- GC_bs[,1] * IX
    GC2X_u <- GC_bs[,2] * IX
    GC3X_u <- GC_bs[,3] * IX
    SmallX_u <- (X_u^d) * IX
    BigX_u <- (X_u^d) * (1-IX)
    
    fit_trunc_adj <- rlm( log(mu_u) ~ SmallX_u + BigX_u + MX_u +
        GC1X_u + GC2X_u + GC3X_u, weights=n_u/sum(n_u) )
        
        
    # return estimates

    GC_bs <- bs( GC, knots=GC_knots, degree=1 ) 

    IX_all <- rep(0,length(X))
    IX_all[ X <= s ] <- 1
    
    MX_u_all <- log2(M+1) * IX_all
    GC1X_u_all <- GC_bs[,1] * IX_all
    GC2X_u_all <- GC_bs[,2] * IX_all
    GC3X_u_all <- GC_bs[,3] * IX_all
    SmallX_u_all <- (X^d) * IX_all
    BigX_u_all <- (X^d) * (1-IX_all)
    
    
    coef_rlm <- coef(fit_trunc_adj)
    muEst <- exp( coef_rlm[1] + coef_rlm[2]*SmallX_u_all + coef_rlm[3]*BigX_u_all +
        coef_rlm[4]*MX_u_all + 
        coef_rlm[5]*GC1X_u_all + coef_rlm[6]*GC2X_u_all + coef_rlm[7]*GC3X_u_all )

    pi0 <- .getPi0( muEst, a_wq_strata, parEst$Y_freq, bgEst=bgEst )
    
    betaEst <- coef(fit_trunc_adj)
    names(betaEst) <- c( "(intercept)", "control(small)", "control(large)", "log2(M+1)",
        "spline(GC)_1", "spline(GC)_2", "spline(GC)_3" )
    
    return( list( pi0 = pi0, a = a_wq_strata, muEst=muEst, betaEst=betaEst,
            Y_val = parEst$Y_val, Y_freq = parEst$Y_freq ) )
}
