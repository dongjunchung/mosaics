
.adapGridMosaicsZ0_OS <- function( Y, M, GC, bgEst=NA,
    min_n_MGC=50, grids_MGC=c(0.01,0.02,0.04,0.10,0.20,0.50),
    parallel=parallel, nCore=nCore )
{        
    M_u <- GC_u <- a_u <- b_u <- mean0_u <- var0_u <-
        u0_u <- u1_u <- u2_u <- n_u <- ty_u <- unitM_u <- c()
        
    Y_freq <- table(Y)
    
    
    # iteration of adaptive griding
    
    for ( i in 1:length(grids_MGC) )
    {
        message( "Info: grid = ", grids_MGC[i] )   
        
        if ( i==1 )
        {
            # initialization
    
            Y2 <- Y
            M2 <- M
            GC2 <- GC
        } else
        {
            # check (M,GC) pair satisfying the condition
            # if condition satisfied for (M,GC) pair,
            # then keep results and remove (M,GC) pair
            # otherwise, retain (M,GC) pair
            
            Y2 <- M2 <- GC2 <- c()
            N <- length(par_est2$M_u)
            
            for ( j in 1:N )
            {
                if ( par_est2$n_u[j] < min_n_MGC )
                {
                    Y_sub_j <- par_est2$Y_sub_list[[j]]
                    nj <- length(Y_sub_j)
                    
                    Y2 <- c( Y2, Y_sub_j )
                    M2 <- c( M2, rep( par_est2$M_u[j], nj ) )
                    GC2 <- c( GC2, rep( par_est2$GC_u[j], nj ) )
                } else
                {
                    M_u <- c( M_u, par_est2$M_u[j] )
                    GC_u <- c( GC_u, par_est2$GC_u[j] )
                    a_u <- c( a_u, par_est2$a_u[j] )
                    b_u <- c( b_u, par_est2$b_u[j] )
                    mean0_u <- c( mean0_u, par_est2$mean0_u[j] )
                    var0_u <- c( var0_u, par_est2$var0_u[j] )
                    n_u <- c( n_u, par_est2$n_u[j] )
                    u0_u <- c( u0_u, par_est2$u0_u[j] )
                    u1_u <- c( u1_u, par_est2$u1_u[j] )
                    u2_u <- c( u2_u, par_est2$u2_u[j] )
                    ty_u <- c( ty_u, par_est2$ty_u[j] )
                    unitM_u <- c( unitM_u, grids_MGC[(i-1)] )
                }
            }
            
            
            # redefine grids
            
            if ( length(M2) > 0 )
            {   
                unitM <- grids_MGC[i]
                M3 <- M2
                GC3 <- GC2
                
                for ( j in 1:((1/unitM)-1) )
                {
                    M2.j <- M2[ M2>=unitM*(j-1) & M2<unitM*j ] 
                    GC2.j <- GC2[ GC2>=unitM*(j-1) & GC2<unitM*j ]
                    
                    M3[ M2>=unitM*(j-1) & M2<unitM*j ] <- median(M2.j)
                    GC3[ GC2>=unitM*(j-1) & GC2<unitM*j ] <- median(GC2.j)
                }
                
                M2.j <- M2[ M2>=unitM*((1/unitM)-1) & M2<=unitM*(1/unitM) ] 
                GC2.j <- GC2[ GC2>=unitM*((1/unitM)-1) & GC2<=unitM*(1/unitM) ]
                
                M3[ M2>=unitM*((1/unitM)-1) & M2<=unitM*(1/unitM) ] <- median(M2.j)
                GC3[ GC2>=unitM*((1/unitM)-1) & GC2<=unitM*(1/unitM) ] <- median(GC2.j)
                
                M2 <- M3
                GC2 <- GC3
                rm( M3, GC3 )
            }
        }
        
        
        # if there is no more (M,GC) pair, then stop iterations
        
        if ( length(Y2)==0 )
        {
            break
        }
        
        
        # background fit
        
        par_est2 <- .mosaicsZ0( Y=Y2, bgEst=bgEst, analysisType="OS", 
            M=M2, GC=GC2, Y_freq=Y_freq,
            parallel=parallel, nCore=nCore )
        
        smallN <- which( par_est2$n_u < min_n_MGC )
        n.uns <- length(smallN) / length(par_est2$n_u)
        #print( paste("percentage of (M,GC) pair with small n:",n.uns) )
        #print( "distribution of n:" )
        #print( quantile( par_est2$n_u, seq(0,1,0.05) ) )
        
        
        # if it is the last iteration, keep the remaining
        
        if ( i==length(grids_MGC) )
        {
            M_u <- c( M_u, par_est2$M_u )
            GC_u <- c( GC_u, par_est2$GC_u )
            a_u <- c( a_u, par_est2$a_u )
            b_u <- c( b_u, par_est2$b_u )
            mean0_u <- c( mean0_u, par_est2$mean0_u )
            var0_u <- c( var0_u, par_est2$var0_u )
            n_u <- c( n_u, par_est2$n_u )
            u0_u <- c( u0_u, par_est2$u0_u )
            u1_u <- c( u1_u, par_est2$u1_u )
            u2_u <- c( u2_u, par_est2$u2_u )
            ty_u <- c( ty_u, par_est2$ty_u )
            unitM_u <- c( unitM_u, rep( grids_MGC[i], length(par_est2$M_u) ) )
        }
    }
    
    
    # return object
    
    par_est_final <- list( M_u = M_u, GC_u = GC_u, a_u = a_u, b_u = b_u,
        mean0_u = mean0_u, var0_u = var0_u,
        u0_u = u0_u, u1_u = u1_u, u2_u = u2_u, n_u = n_u, ty_u = ty_u,
        Y_val = as.numeric(names(table(Y))), Y_freq = table(Y), unitM_u = unitM_u )
    return( par_est_final )
}
