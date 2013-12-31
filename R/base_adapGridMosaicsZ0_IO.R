
#.adapGridMosaicsZ0_IO <- function( Y, X, min_n_X=50 )
.adapGridMosaicsZ0_IO <- function( Y, X, bgEst=NA, inputTrunc, min_n_X=50,
    parallel=parallel, nCore=nCore )
{        
    X_u <- a_u <- b_u <- mean0_u <- var0_u <-
        u0_u <- u1_u <- u2_u <- n_u <- ty_u <- c()
    
    Y_freq <- table(Y)
    
    
    # adaptive griding for X (Input)
    
    X_set <- sort( unique(X), decreasing=TRUE )
    ind_X_set <- rep( 0, length(X_set) )
    
    ind_now <- 1
    N_now <- 0
    
    for ( i in 1:length(X_set) )
    {
        N_i <- length( which( X==X_set[i] ) )
        if ( N_now <= min_n_X )
        {
            ind_X_set[i] <- ind_now
            N_now <- N_now + N_i
        } else
        {
            ind_now <- ind_now + 1
            ind_X_set[i] <- ind_now
            N_now <- N_i
        }
    }
    
    X_set_new <- rep( 0, length(X_set) )    
    for ( i in 1:length(unique(ind_X_set)) )
    {
        X_set_new[ind_X_set==i] <- median( X_set[ind_X_set==i] )
    }
    
    X_new <- rep( 0, length(X) )    
    for ( i in 1:length(X_set) )
    {
        X_new[ X==X_set[i] ] <- X_set_new[i]
    }
    
    
    # background fit
        
    par_est2 <- .mosaicsZ0( Y=Y, bgEst=bgEst, analysisType="IO", 
        X=X_new, inputTrunc=inputTrunc, Y_freq=Y_freq,
        parallel=parallel, nCore=nCore )
    
    
    # return object
    
    par_est_final <- list( X_u = par_est2$X_u, a_u = par_est2$a_u, b_u = par_est2$b_u,
        mean0_u = par_est2$mean0_u, var0_u = par_est2$var0_u,
        u0_u = par_est2$u0_u, u1_u = par_est2$u1_u, u2_u = par_est2$u2_u, 
        n_u = par_est2$n_u, ty_u = par_est2$ty_u,
        Y_val = as.numeric(names(table(Y))), Y_freq = table(Y) )
    return( par_est_final )
}
