
#########################################################
# Plot GOF for MOSAiCS-HMM
#########################################################

.plotGOFHMM <- function( mosaicsHMMEst, mosaicsEst, tagCount, input, signalModel="2S", 
	k=3, seed=12345, parallel=FALSE, nCore=8 )
{    		
    # extract MOSAiCS parameter estimates
	# generate histogram of ChIP (& input) data
	
    analysisType <- mosaicsEst@analysisType
    
    YFreq <- table(tagCount)
    YVal <- as.numeric(names(YFreq))
    
    if ( analysisType=="TS" | analysisType=="IO" ) {
        XFreq <- table( input )
        XVal <- as.numeric(names(XFreq))
    }
        
    N <- sum(YFreq)

    a <- mosaicsEst@a
    pi0 <- mosaicsEst@pi0
    muEst <- mosaicsEst@muEst
    bEst <- a / muEst
    
    b <- mosaicsEst@b
    c <- mosaicsEst@c
    
    p1 <- mosaicsEst@p1
    b1 <- mosaicsEst@b1
    c1 <- mosaicsEst@c1
    b2 <- mosaicsEst@b2
    c2 <- mosaicsEst@c2
    
    # simulate data: null model
    
    message( "Info: simulating background tags..." )
    
    probZ_2S <- c( pi0, (1-pi0)*p1, (1-pi0)*(1-p1) )
    set.seed(12345)
    Zsim_2S <- sample( c(0,1,2), size=N, prob=probZ_2S, replace=TRUE )
    set.seed(12345)
    Ysim0 <- rnbinom( N, a, bEst/(bEst+1) )   
    
    YsimZ0Freq <- table( Ysim0[ Zsim_2S==0 ] )
    YsimZ0Val <- as.numeric(names(YsimZ0Freq))  
    
    rm( Zsim_2S )
    gc()  
    
    # simulate data: null + signal model
    
    if ( signalModel == "1S" ) {
	    # simulate one sample signal
	    
	    # MOSAiCS model
	    
	    message( "Info: simulating ChIP tags for MOSAiCS..." )
	    
	    probZ_1S <- c( pi0, (1-pi0) )
	    set.seed(12345)
	    Zsim_1S <- sample( c(0,1), size=N, prob=probZ_1S, replace=TRUE )
	    set.seed(12345)
	    Ysim <- rnbinom( N, a, bEst/(bEst+1) )
	    nS_1S <- length(which( Zsim_1S==1 ))
	    Ysim[ Zsim_1S==1 ] <- Ysim[ Zsim_1S==1 ] + rnbinom( nS_1S, b, c/(c+1) ) + k
	    
	    YsimFreq <- table(Ysim)
	    YsimVal <- as.numeric(names(YsimFreq))
	    
	    rm( Zsim_1S, nS_1S )
	    gc()
	    
	    # MOSAiCS-HMM model
	    
	    message( "Info: simulating ChIP tags for MOSAiCS-HMM..." )
	    
		if ( parallel ) {
	        Ysim_HMM <- mclapply( mosaicsHMMEst, 
				function(x) .simulateHMM_1S( 
					mosaicsHMMEst_chr=x, a=a, bEst=bEst, 
					b=b, c=c, k=k, seed=seed ),
				mc.cores=nCore )
	    } else {
	        Ysim_HMM <- lapply( mosaicsHMMEst, 
				function(x) .simulateHMM_1S( 
					mosaicsHMMEst_chr=x, a=a, bEst=bEst, 
					b=b, c=c, k=k, seed=seed ) )		
		}
	    
	    YsimFreq_HMM <- table(unlist(Ysim_HMM))
	    YsimVal_HMM <- as.numeric(names(YsimFreq_HMM))
		
	} else if ( signalModel == "2S" ) {
	    # simulation two sample signal
	    
	    # MOSAiCS model
	    
	    message( "Info: simulating ChIP tags for MOSAiCS..." )
	    
	    probZ_2S <- c( pi0, (1-pi0)*p1, (1-pi0)*(1-p1) )
	    set.seed(12345)
	    Zsim_2S <- sample( c(0,1,2), size=N, prob=probZ_2S, replace=TRUE )
	    set.seed(12345)
	    Ysim <- rnbinom( N, a, bEst/(bEst+1) )   
	    nS_2S_1 <- length(which(Zsim_2S==1))
	    nS_2S_2 <- length(which(Zsim_2S==2))
	    Ysim[ Zsim_2S==1 ] <- Ysim[ Zsim_2S==1 ] + rnbinom( nS_2S_1, b1, c1/(c1+1) ) + k
	    Ysim[ Zsim_2S==2 ] <- Ysim[ Zsim_2S==2 ] + rnbinom( nS_2S_2, b2, c2/(c2+1) ) + k 
	    
	    YsimFreq <- table(Ysim)
	    YsimVal <- as.numeric(names(YsimFreq))
	    
	    rm( Zsim_2S, nS_2S_1, nS_2S_2 )
	    gc()
	    
	    # MOSAiCS-HMM model
	    
	    message( "Info: simulating ChIP tags for MOSAiCS-HMM..." )
	    
		if ( parallel ) {			
			Ysim_HMM <- mclapply( mosaicsHMMEst, 
				function(x) .simulateHMM_2S( 
					mosaicsHMMEst_chr=x, a=a, bEst=bEst, 
					p1=p1, b1=b1, c1=c1, b2=b2, c2=c2, k=k, seed=seed ),
				mc.cores=nCore )
	    } else {			
			Ysim_HMM <- lapply( mosaicsHMMEst, 
				function(x) .simulateHMM_2S( 
					mosaicsHMMEst_chr=x, a=a, bEst=bEst, 
					p1=p1, b1=b1, c1=c1, b2=b2, c2=c2, k=k, seed=seed )
				)		
		}
	    
	    YsimFreq_HMM <- table(unlist(Ysim_HMM))
	    YsimVal_HMM <- as.numeric(names(YsimFreq_HMM))
	}	
    
    # draw GOF plot
    
    plot( log10(YVal+1),log10(YFreq), type='l', ylab='Frequency', xlab='Tag count', axes=FALSE )  
    if ( analysisType=="TS" | analysisType=="IO" ) {
        points( log10(XVal+1),log10(XFreq), type='l', col='darkgray' )  
    }
    points( log10(YsimZ0Val+1), log10(YsimZ0Freq), type='l', col='green' )
    points( log10(YsimVal+1), log10(YsimFreq), type='l', col='red' )  
    points( log10(YsimVal_HMM+1), log10(YsimFreq_HMM), type='l', col='blue' )

    if ( analysisType=="TS" | analysisType=="IO" ) {
        legend( 1, log10(max(YFreq)+1), 
            c('Actual data (ChIP)','Actual data (Control)','Sim:Null','Sim:MOSAiCS','Sim:MOSAiCS-HMM'),
            col=c('black','darkgray','green','red','blue'),lty=c(1,1,1,1,1),bty='n')
    } else {
        legend( 1, log10(max(YFreq)+1), 
            c('Actual data','Sim:Null','Sim:MOSAiCS','Sim:MOSAiCS-HMM'),
            col=c('black','green','red','blue'),lty=c(1,1,1,1),bty='n')
    }
    abline(v=log10(0+1),lty=2,col='gray')
    abline(v=log10(1+1),lty=2,col='gray')
    abline(v=log10(2+1),lty=2,col='gray')
    axis(1,0:6,10^c(0:6)-1)
    axis(2,0:10,10^c(0:10))
}

.simulateHMM_1S <- function( mosaicsHMMEst_chr, a, bEst, b, c, k, seed ) {

	# extract MOSAiCS-HMM parameter estimates
	
	piMat <- mosaicsHMMEst_chr$piMat
	pi0vec <- mosaicsHMMEst_chr$pi0vec
	gMat <- mosaicsHMMEst_chr$gMat_chr
	chrlen <- ncol(gMat)
	
	# simulate hidden states
	
	Zsim_HMM_1S_chr <- rep( NA, chrlen )
	set.seed(seed)
	Zsim_HMM_1S_chr[1] <- sample( c(0,1), 1, prob = pi0vec )
	
	set.seed(seed)
	for ( j in 2:chrlen ) {
		if ( Zsim_HMM_1S_chr[(j-1)] == 0 ) {
			Zsim_HMM_1S_chr[j] <- sample( c(0,1), 1, prob = piMat[1,] )
		} else {
			Zsim_HMM_1S_chr[j] <- sample( c(0,1), 1, prob = piMat[2,] )
		}
	}
	
	# simulate ChIP tag counts
	
	set.seed(seed)
	Ysim_HMM_1S_chr <- rnbinom( chrlen, a, bEst/(bEst+1) )
	nS_HMM_1S_chr <- length(which( Zsim_HMM_1S_chr == 1 ))
	Ysim_HMM_1S_chr[ Zsim_HMM_1S_chr == 1 ] <- 
		Ysim_HMM_1S_chr[ Zsim_HMM_1S_chr == 1 ] + rnbinom( nS_HMM_1S_chr, b, c/(c+1) ) + k
	
	return( Ysim_HMM_1S_chr )
}

.simulateHMM_2S <- function( mosaicsHMMEst_chr, a, bEst, p1, b1, c1, b2, c2, k, seed ) {
			
	# extract MOSAiCS-HMM parameter estimates
	
	piMat <- mosaicsHMMEst_chr$piMat
	pi0vec <- mosaicsHMMEst_chr$pi0vec
	gMat <- mosaicsHMMEst_chr$gMat_chr
	Nchr <- ncol(gMat)
	
	# simulate hidden states
	
	Zsim_HMM_2S_chr <- rep( NA, Nchr )
	set.seed(12345)
	Zsim_HMM_2S_chr[1] <- sample( c(0,1), 1, prob = pi0vec )
	
	set.seed(12345)
	for ( j in 2:Nchr ) {
		if ( Zsim_HMM_2S_chr[(j-1)] == 0 ) {
			Zsim_HMM_2S_chr[j] <- sample( c(0,1), 1, prob = piMat[1,] )
		} else {
			Zsim_HMM_2S_chr[j] <- sample( c(0,1), 1, prob = piMat[2,] )
		}
	}
	
	# simulate ChIP tag counts
	
	Nsig <- length(which( Zsim_HMM_2S_chr != 0 ))				
	probZ_2S <- c( p1, (1-p1) )
	set.seed(12345)
	Zsim_HMM_2S_chr[ Zsim_HMM_2S_chr != 0 ] <-
		sample( c(1,2), size=Nsig, prob=probZ_2S, replace=TRUE )
	set.seed(12345)
	Ysim_HMM_2S_chr <- rnbinom( Nchr, a, bEst/(bEst+1) )   
	nS_HMM_2S_1_chr <- length(which( Zsim_HMM_2S_chr == 1 ))
	nS_HMM_2S_2_chr <- length(which( Zsim_HMM_2S_chr == 2 ))
	Ysim_HMM_2S_chr[ Zsim_HMM_2S_chr == 1 ] <- 
		Ysim_HMM_2S_chr[ Zsim_HMM_2S_chr == 1 ] + rnbinom( nS_HMM_2S_1_chr, b1, c1/(c1+1) ) + k
	Ysim_HMM_2S_chr[ Zsim_HMM_2S_chr == 2 ] <-
		Ysim_HMM_2S_chr[ Zsim_HMM_2S_chr == 2 ] + rnbinom( nS_HMM_2S_2_chr, b2, c2/(c2+1) ) + k 

	return( Ysim_HMM_2S_chr )
}
