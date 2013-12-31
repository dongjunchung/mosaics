
#########################################################
# Z1 with mixture of 2 signal components
#########################################################

.mosaicsZ1_2S <- function( MOSAiCS_Z0, Y, pNfit, Y_bd_all, k=3 )
{
    ##################################################################
    ### initialization of the main EM
    ##################################################################
    
    a <- MOSAiCS_Z0$a    
    #mu_est <- get_mu_est(M, GC, MOSAiCS_Z0, trunc_GC=0.55)
    mu_est <- MOSAiCS_Z0$muEst
    b_est <- a / mu_est
    pi0 <- MOSAiCS_Z0$pi0    
    tab_Y <- MOSAiCS_Z0$Y_freq
    Y_u <- MOSAiCS_Z0$Y_val
    
    
    # initial EM calculation

    par_2mixNB <- .emZ1_2S( Y_bd_all, epsilon=0.0001, 
        b1_old=NULL, c1_old=NULL, b2_old=NULL, c2_old=NULL )

    mu1_init <- par_2mixNB$b1/par_2mixNB$c1
    mu2_init <- par_2mixNB$b2/par_2mixNB$c2
    
    var1_init <- mu1_init*(1+1/par_2mixNB$c1)
    var2_init <- mu2_init*(1+1/par_2mixNB$c2)
    
    
    # calculate initial EN and varN, based on initial EM calculation

    if( min(mu1_init,mu2_init) > mean(mu_est) )
    {
        EN <- mean(mu_est)
    } else
    {
        if( min(mu1_init,mu2_init) > median(mu_est) )
        {
            EN <- median(mu_est)
        } else
        {
            EN <- min(mu1_init,mu2_init) - 0.1 
        }
    }

    varN <- EN*(1 + EN/a)
 
    if( min(var1_init,var2_init) < varN )
    {
        varN <- min(var1_init,var2_init) - 0.15
    }
    
    
    # initialize b1, c1, b2, and c2
    
    if( var1_init - varN - mu1_init + EN <= 0 )
    {
        b1_init <- (mu1_init - EN)^2 / 0.1
        c1_init <- (mu1_init - EN) / 0.1  
    } else
    {
        b1_init <- (mu1_init - EN)^2 / (var1_init - varN - mu1_init + EN)
        c1_init <- (mu1_init - EN) / (var1_init - varN - mu1_init + EN)    
    }
       
    if( var2_init - varN - mu2_init + EN <= 0 )
    {
        b2_init <- (mu2_init - EN)^2 / 0.1
        c2_init <- (mu2_init - EN) / 0.1  
    } else
    {
        b2_init <- (mu2_init - EN)^2 / (var2_init - varN - mu2_init + EN)
        c2_init <- (mu2_init - EN) / (var2_init - varN - mu2_init + EN)    
    }   
    
    
    # take care of identifiability problem

    mS1 <- b1_init/c1_init
    mS2 <- b2_init/c2_init
    
    if(mS1 <= mS2)
    {
        # OK. Good to go.
        
        p1_init <- par_2mixNB$p1
        
        b1_init <- b1_init
        c1_init <- c1_init
        b2_init <- b2_init
        c2_init <- c2_init
    } else
    {
        # Let's swap them.
        
        p1_init <- par_2mixNB$p2    # p2 = 1 - p1
        
        b1_init_temp <- b2_init
        c1_init_temp <- c2_init
        b2_init_temp <- b1_init
        c2_init_temp <- c1_init
        
        b1_init <- b1_init_temp
        c1_init <- c1_init_temp
        b2_init <- b2_init_temp
        c2_init <- c2_init_temp
    }


    ##################################################################
    ### The main EM calculation: iterate till convergence
    ##################################################################
    
    # calculate EN and vanN
    # (redefine EN after heuristic EN adjustment for initialization)

    EN <- mean(mu_est) 
    varN <- EN*(1 + EN/a) 
    
    
    # parameters for Y>=k
    
    id_geqk <- which(Y>=k)      
    Y_ori <- Y[id_geqk]
    Y_Z1 <- Y - k               # adjusted Y
    Y_Z1 <- Y_Z1[id_geqk]
    #M_Z1 <- M[id_geqk]
    #GC_Z1 <- GC[id_geqk]
    b_est_Z1 <- b_est[id_geqk]
    mu_est_Z1 <- mu_est[id_geqk]
    
    Yk <- Y_ori - k        # use only Y >= k
    #if(length(which(Y<0))>0 ) Y[which(Y<0)] <- -1 
    Ykmax <- max(Yk)
    #ind_ge_k <- which(Y>=0)
    
    
    # Initialization of the main EM calculation
    
    #PYZ0 <- dnbinom( Y_ori, a, b_est_Z1/(b_est_Z1+1) )
    PYZ0 <- pNfit$PYZ0
    #pNfit <- .calcPN( Y_ori, k, a, mu_est_Z1 ) 
    #PYZ1 <- .margDistZ1_2S( Y_ori, pNfit, b1_init, c1_init, b2_init, c2_init )
    PYZ1 <- .margDistZ1_2S( Yk, Ykmax, pNfit, b1_init, c1_init, b2_init, c2_init )

    PYZ1G1 <- PYZ1$MDG1
    PYZ1G2 <- PYZ1$MDG2
        
    b1_iter <- b1_init
    c1_iter <- c1_init
    b2_iter <- b2_init
    c2_iter <- c2_init
    p1 <- p1_init
    p1_iter <- p1_init
    
    
    # main EM iteration
    
    eps <- 1e-6
    logLik1 <- sum( log( pi0*PYZ0 + (1-pi0)*( p1*PYZ1G1 + (1-p1)*PYZ1G2) ) )
    logLik <- c( -Inf, logLik1 )
    iter <- 2
    
    while( abs(logLik[iter]-logLik[iter-1])>eps & iter < 10 )
    {
        # E-step
        
        denom_G_latent <- p1*PYZ1G1 + (1-p1)*PYZ1G2

        G_latent <- p1*PYZ1G1 / denom_G_latent
        Z_latent <- (1 - pi0)*denom_G_latent / ( pi0*PYZ0 + (1 - pi0)*denom_G_latent )
        
        # M-step
        
        p1 <- sum(G_latent*Z_latent) / sum(Z_latent)

        mu1 <- sum(Z_latent*G_latent*Y_Z1) / sum(Z_latent*G_latent)
        var1 <- sum(Z_latent*G_latent*(Y_Z1 - mu1)^2) / sum(Z_latent*G_latent)

        mu2 <- sum(Z_latent*(1-G_latent)*Y_Z1) / sum(Z_latent*(1-G_latent))
        var2 <- sum(Z_latent*(1-G_latent)*(Y_Z1 - mu2)^2) / sum(Z_latent*(1-G_latent))

        b1 <- (mu1 - EN)^2 / (var1 - varN - mu1 + EN)
        c1 <- (mu1 - EN) / (var1 - varN - mu1 + EN)
        
        b2 <- (mu2 - EN)^2 / (var2 - varN - mu2 + EN)
        c2 <- (mu2 - EN) / (var2 - varN - mu2 + EN)
        
        # stop iteration if assumptions are not satisfied
        
        if ( p1<0.01 | b1<0 | c1<0 | b2<0 | c2<0 )
        {
            p1 <- p1_iter[(iter-1)]
            b1 <- b1_iter[(iter-1)]
            c1 <- c1_iter[(iter-1)]
            b2 <- b2_iter[(iter-1)]
            c2 <- c2_iter[(iter-1)]
            break
        }
        
        # calculate P(Y|Z=1,G=1) and P(Y|Z=1,G=2)
        
        #print( "calculate P(Y|Z=1,G=1) and P(Y|Z=1,G=2)" )

        #PYZ1 <- .margDistZ1_2S( Y_ori, pNfit, b1, c1, b2, c2)
        PYZ1 <- .margDistZ1_2S( Yk, Ykmax, pNfit, b1, c1, b2, c2)

        PYZ1G1 <- PYZ1$MDG1
        PYZ1G2 <- PYZ1$MDG2
        
        # update iteration

        logLik_t <- sum( log( pi0*PYZ0 + (1-pi0)*( p1*PYZ1G1 + (1-p1)*PYZ1G2) ) )
        if ( is.na(logLik_t) | is.nan(logLik_t) ) {
            p1 <- p1_iter[(iter-1)]
            b1 <- b1_iter[(iter-1)]
            c1 <- c1_iter[(iter-1)]
            b2 <- b2_iter[(iter-1)]
            c2 <- c2_iter[(iter-1)]
            logLik_t <- logLik[(iter-1)]
            break
        }
        logLik <- c( logLik, logLik_t )
        p1_iter <- c(p1_iter, p1)
        b1_iter <- c(b1_iter, b1)
        c1_iter <- c(c1_iter, c1)
        b2_iter <- c(b2_iter, b2)
        c2_iter <- c(c2_iter, c2)

        iter <- iter + 1
    }
    
    return( list(
        pi0 = MOSAiCS_Z0$pi0, a = MOSAiCS_Z0$a, muEst = MOSAiCS_Z0$muEst,
        p1 = p1, b1 = b1, c1 = c1, b2 = b2, c2 = c2 ) )
}


# calculate P(Y|Z=1,G=1) and P(Y|Z=1,G=2) for the current parameters

#.margDistZ1_2S <- function( Yori, pNfit, b1, c1, b2, c2 )
.margDistZ1_2S <- function( Y, Ymax, pNfit, b1, c1, b2, c2 )
# Y <- Yori - k
{   
    k <- pNfit$k
    pN <- pNfit$pN
    mu_round <- pNfit$mu_round
    mu_round_U <- pNfit$mu_round_U
    
    # process Y
    
    #Y <- Yori - k        # use only Y >= k
    #if(length(which(Y<0))>0 ) Y[which(Y<0)] <- -1 
    #Ymax <- max(Y)
    #ind_ge_k <- which(Y>=0)
                 
    # prob of S1 & S2
    
    pS1 <- dnbinom( 0:Ymax, b1, c1/(c1+1) )
    pS2 <- dnbinom( 0:Ymax, b2, c2/(c2+1) )
    MDGfit <- conv_2S( y=Y, mu_round=mu_round,
        mu_round_U=mu_round_U, pN=pN, pS1=pS1, pS2=pS2 )
    MDG <- matrix( MDGfit, ncol=2, byrow=TRUE )  
    
    return( list( MDG1 = MDG[,1], MDG2 = MDG[,2] ) )
   
}


#########################################################
# EM algorithm to initialize MOSAiCS_Z1_2S
#########################################################

.calMDZ1_2S <- function(b1,c1,b2,c2,Y_freq,Y_val)
{
    # NB distribution

    MD <- matrix(0,nrow=length(Y_val),ncol=2)
    MD[,1] <- dnbinom(Y_val,size=b1,prob=c1/(c1+1))    
    MD[,2] <- dnbinom(Y_val,size=b2,prob=c2/(c2+1))
   
    id <- which( apply(MD,1,sum)==0 )    
    if(length(id)>0){
        replaceID <- min(id)-1
        MD[id,1] <- MD[replaceID,1]
        MD[id,2] <- MD[replaceID,2]
    }
    return(MD)
}

.eStepZ1_2S <- function(MD,pZ,Y_freq,Y_val)
{

    p1 <- sum(pZ[,1]*Y_freq)/sum(Y_freq)
    p2 <- 1-p1
    
    MD1 <- MD[,1]
    MD2 <- MD[,2]
    
    pZ[,1] <- (MD1*p1)/(MD1*p1+MD2*p2)
    pZ[,2] <- (MD2*p2)/(MD1*p1+MD2*p2)
    
    return(pZ)
}

.mStepZ1_2S <- function(pZ,Y_freq,Y_val,b1_old=NULL,c1_old=NULL,b2_old=NULL,c2_old=NULL)
{
    pZ1 <- pZ[,1]
    pZ2 <- pZ[,2]    
    p1 <- sum(pZ1*Y_freq) / sum(Y_freq)
    p2 <- 1-p1
    
    mu1 <- sum(pZ1*Y_val*Y_freq) / sum(pZ1*Y_freq)
    var1 <- sum(pZ1*Y_freq*(Y_val-mu1)^2) / sum(pZ1*Y_freq)    
    mu2 <- sum(pZ2*Y_val*Y_freq) / sum(pZ2*Y_freq)
    var2 <- sum(pZ2*Y_freq*(Y_val-mu2)^2) / sum(pZ2*Y_freq)   
    
    if ( is.na(mu1) || is.na(var1) || is.nan(mu1) || is.nan(var1) ) {
        # inappropriate values occur
        stop( "over-estimation of background detected. Please tune the parameters!" )
    }
    
    if ( is.na(mu2) || is.na(var2) || is.nan(mu2) || is.nan(var2) ) {
        # inappropriate values occur
        stop( "over-estimation of background detected. Please tune the parameters!" )
    }
    
    if ( var1 < mu1 ) { var1 = mu1 + 0.1 }
    if ( var2 < mu2 ) { var2 = mu2 + 0.1 }    
    
    if ( length(b1_old)==0 ) { b1 <- mu1^2/(var1-mu1) } else { b1 <- b1_old }
    if ( length(c1_old)==0 ) { c1 <- mu1/(var1-mu1) } else { c1 <- c1_old }    
    if ( length(b2_old)==0 ) { b2 <- mu2^2/(var2-mu2) } else { b2 <- b2_old }
    if ( length(c2_old)==0 ) { c2 <- mu2/(var2-mu2) } else { c2 <- c2_old }
        
    MD <- .calMDZ1_2S(b1,c1,b2,c2,Y_freq,Y_val)
    logLik <- sum( log(p1*MD[,1]+p2*MD[,2]) * Y_freq )
    
    return( list( b1=b1, c1=c1, b2=b2, c2=c2, p1=p1, p2=p2, MD=MD, logLik=logLik ) )
}


.emZ1_2S <- function( Y, epsilon, b1_old=NULL, c1_old=NULL, b2_old=NULL, c2_old=NULL )
{    
    Y_freq <- table(Y)
    Y_val <- as.numeric(names(table(Y)))       
    
    # initialize b1, c1, b2, and c2
    
    ind_0.5 <- which( Y <= quantile(Y,0.5) )
    mean1 <- mean(Y[ind_0.5])
    var1 <- var(Y[ind_0.5])
 
    ind_0.8 <- which( Y > quantile(Y,0.8) )
    mean2 <- mean(Y[ind_0.8])
    var2 <- var(Y[ind_0.8])
    
    if(mean1 == 0||is.na(mean1)==TRUE) 
    {
        ind_0.6 <- which( Y <= quantile(Y,0.6) )
        mean1 <- mean(Y[ind_0.6])
        var1 <- var(Y[ind_0.6])
    }

    if(mean1 == 0||is.na(mean1)==TRUE) 
    {
        ind_0.7 <- which( Y <= quantile(Y,0.7) )
        mean1 <- mean(Y[ind_0.7])
        var1 <- var(Y[ind_0.7])
    }

    if(mean1 == 0||is.na(mean1)==TRUE) 
    {
        ind_0.8 <- which( Y <= quantile(Y,0.8) )
        mean1 <- mean(Y[ind_0.8])
        var1 <- var(Y[ind_0.8])
    }

    if(mean1 == 0||is.na(mean1)==TRUE) 
    {
        ind_0.9 <- which( Y <= quantile(Y,0.9) )
        
        mean1 <- mean(Y[ind_0.9])
        var1 <- var(Y[ind_0.9])
        
        ind_0.9 <- which( Y > quantile(Y,0.9) )
        
        mean2 <- mean(Y[ind_0.9])
        var2 <- var(Y[ind_0.9])
    }
    
    if ( is.na(mean1) || is.na(var1) || is.nan(mean1) || is.nan(var1) ) {
        # inappropriate values occur
        stop( "over-estimation of background detected. Please tune the parameters!" )
    }
    
    if ( is.na(mean2) || is.na(var2) || is.nan(mean2) || is.nan(var2) ) {
        # inappropriate values occur
        stop( "over-estimation of background detected. Please tune the parameters!" )
    }
    
    if (var1 < mean1) { var1 <- mean1+0.1 }
    if (length(b1_old)==0) { b1 <- mean1^2/(var1-mean1) } else { b1 <- b1_old }
    if (length(c1_old)==0) { c1 <- mean1/(var1-mean1) } else { c1 <- c1_old }
    if (var2 < mean2) { var2 <- mean2+0.1 }
    if (length(b2_old)==0) { b2 <- mean2^2/(var2-mean2) } else { b2 <- b2_old }
    if (length(c2_old)==0) { c2 <- mean2/(var2-mean2) } else { c2 <- c2_old }
    
    
    # initialize EM

    MD_2mixNB  <- .calMDZ1_2S(b1,c1,b2,c2,Y_freq,Y_val)
    pZ_old <- matrix(0.5,nrow=length(Y_val),ncol=2)
    
    pZ <- .eStepZ1_2S(MD_2mixNB,pZ_old,Y_freq,Y_val)
    pZ_old <- pZ
    par_2mixNB <- .mStepZ1_2S(pZ,Y_freq,Y_val,b1_old,c1_old,b2_old,c2_old)
    
    logLik <- par_2mixNB$logLik
    absdif <- Inf
    i <- 2
    
        
    # EM iterations

    while( absdif > epsilon & i<=10000 )
    {
        # E-step
        
        par_2mixNB_old <- par_2mixNB
        pZ <- .eStepZ1_2S(par_2mixNB$MD,pZ_old,Y_freq,Y_val)
        
        # M-step
        
        pZ_old <- pZ
        par_2mixNB <- .mStepZ1_2S(pZ,Y_freq,Y_val,b1_old,c1_old,b2_old,c2_old)
        
        # update log lik    
        
        logLik_old <- logLik
        if(i==10000)
        {            
            logLik <- logLik_old
        } else
        {
            logLik <- par_2mixNB$logLik
        }
        absdif <- abs( logLik - logLik_old )
        
        i=i+1   
    }
    par_2mixNB$logLik <- logLik
    
    return(par_2mixNB)
}
