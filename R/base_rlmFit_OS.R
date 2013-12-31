
.rlmFit_OS <- function( parEst, mean_thres=0, Y, M, GC )
{    
    #library(MASS)
    #library(splines)
    
    
    # exclude NA & MM points
    
    if ( length(which(is.na(parEst$a_u)==TRUE))==length(parEst$a_u) ) {
        stop( "insufficient # of proper strata! Cannot proceed!" )
    }
    
    # do not excluded estimtes from MOM (ver 1.0.2) 
    
    #idNA <- unique( c( which(is.na(parEst$a_u)==TRUE), which(as.character(parEst$ty_u)=='MM') ) )
    #idNA <- which( is.na(parEst$a_u) | as.character(parEst$ty_u)=='MM' )
    idNA <- which( is.na(parEst$a_u) )
    if ( length(idNA)>0 ) {
        a_u <- parEst$a_u[-idNA]
        mu_u <- parEst$a_u[-idNA] / parEst$b_u[-idNA]
        M_u <- parEst$M_u[-idNA]
        GC_u <- parEst$GC_u[-idNA]
        n_u <- parEst$n_u[-idNA]
        mean0 <- parEst$mean0_u[-idNA]
        var0 <- parEst$var0_u[-idNA]
    } else {
        a_u <- parEst$a_u
        mu_u <- parEst$a_u / parEst$b_u
        M_u <- parEst$M_u
        GC_u <- parEst$GC_u
        n_u <- parEst$n_u
        mean0 <- parEst$mean0_u
        var0 <- parEst$var0_u 
    }
    
    
    # exclude mu>Yhat points
    
    Y_freq <- table(Y)
    
    if( sum(Y_freq[ as.numeric(names(Y_freq))<=2 ])/sum(Y_freq) > 0.5 ) {    
	id_bigmean <- which( ( mu_u - mean0 ) > mean_thres )
	if( length(id_bigmean)>0 ) {
		a_u <- a_u[-id_bigmean]
		mu_u <- mu_u[-id_bigmean]
		M_u <- M_u[-id_bigmean]
		GC_u <- GC_u[-id_bigmean]
		n_u <- n_u[-id_bigmean]
		mean0 <- mean0[-id_bigmean]
		var0 <- var0[-id_bigmean]
	}
    }
    
    
    # estimators of a

    if( sum(Y_freq[ as.numeric(names(Y_freq))<=2 ])/sum(Y_freq) > 0.5 ) {
	t_n_sq <- length(a_u)
	
	idNA <- which( a_u>quantile(a_u,prob=0.95) | a_u<quantile(a_u,prob=0.05) )
	if(length(idNA)>0){
		a_u_wq <- a_u[-idNA]
		n_u_wq <- n_u[-idNA]
		a_wq_strata <- sum(n_u_wq*a_u_wq) / sum(n_u_wq)
	}    
    } else {
    	a_wq_strata <- sum(n_u*a_u)/sum(n_u)
    }
    
    # sample size weighted rlm fit    
    
    GC_knots <- c(quantile(GC_u,prob=0.25),quantile(GC_u,prob=0.75))
    GC_u_bs <- bs(GC_u,degree=1,knots=GC_knots)           # [Note_2] different GC_knots
    GC_bs <- bs( GC, knots=GC_knots, degree=1 ) 
    
    fit_trunc_adj <- rlm( log(mu_u) ~ log2(M_u+1) + GC_u_bs, weights=n_u/sum(n_u) )
    #print(fit_trunc_adj)    
        
        
    # return estimates
    
    coef_rlm <- coef(fit_trunc_adj)
    muEst <- exp( coef_rlm[1] + coef_rlm[2]*log2(M+1) +
        coef_rlm[3]*GC_bs[,1] + coef_rlm[4]*GC_bs[,2] + coef_rlm[5]*GC_bs[,3] )

    pi0 <- .getPi0( muEst, a_wq_strata, parEst$Y_freq )
    
    betaEst <- coef(fit_trunc_adj)
    names(betaEst) <- c( "(intercept)", "log2(M+1)", "spline(GC)_1", "spline(GC)_2", "spline(GC)_3" )

    return( list( pi0 = pi0, a = a_wq_strata, muEst=muEst, betaEst=betaEst,
        Y_val=parEst$Y_val, Y_freq=parEst$Y_freq ))
    
}
