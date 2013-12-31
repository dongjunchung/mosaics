
#########################################################
# Plot GOF
#########################################################

.plotGOF <- function( mosaicsEst, tagCount, input=0, k=3 )
{    
    analysisType <- mosaicsEst@analysisType
    
    # calculate basic parameter estimates
    
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
    
    # simulate one sample signal
    
    probZ_1S <- c( pi0, (1-pi0) )
    set.seed(12345)
    Zsim_1S <- sample( c(0,1), size=N, prob=probZ_1S, replace=TRUE )
    set.seed(12345)
    Ysim_1S <- rnbinom( N, a, bEst/(bEst+1) )
    nS_1S <- length(which( Zsim_1S==1 ))
    Ysim_1S[ Zsim_1S==1 ] <- Ysim_1S[ Zsim_1S==1 ] + rnbinom( nS_1S, b, c/(c+1) ) + k
    
    #YsimZ0Freq <- table( Ysim_1S[ Zsim_1S==0 ] )
    #YsimZ0Val <- as.numeric(names(YsimZ0Freq))
    
    YsimFreq_1S <- table(Ysim_1S)
    YsimVal_1S <- as.numeric(names(YsimFreq_1S))
    
    # simulation two sample signal
    
    probZ_2S <- c( pi0, (1-pi0)*p1, (1-pi0)*(1-p1) )
    set.seed(12345)
    Zsim_2S <- sample( c(0,1,2), size=N, prob=probZ_2S, replace=TRUE )
    set.seed(12345)
    Ysim_2S <- rnbinom( N, a, bEst/(bEst+1) )   
    nS_2S_1 <- length(which(Zsim_2S==1))
    nS_2S_2 <- length(which(Zsim_2S==2))
    Ysim_2S[ Zsim_2S==1 ] <- Ysim_2S[ Zsim_2S==1 ] + rnbinom( nS_2S_1, b1, c1/(c1+1) ) + k
    Ysim_2S[ Zsim_2S==2 ] <- Ysim_2S[ Zsim_2S==2 ] + rnbinom( nS_2S_2, b2, c2/(c2+1) ) + k 
    
    YsimZ0Freq <- table( Ysim_2S[ Zsim_2S==0 ] )
    YsimZ0Val <- as.numeric(names(YsimZ0Freq))    

    #Y_sim_Z0_freq <- table(Y_sim[which(Z_sim == 0)])
    #Y_sim_Z0_val <- as.numeric(names(Y_sim_Z0_freq))
    
    YsimFreq_2S <- table(Ysim_2S)
    YsimVal_2S <- as.numeric(names(YsimFreq_2S))
    
    # draw GOF plot
    
    plot( log10(YVal+1),log10(YFreq), type='l', ylab='Frequency', xlab='Tag count', axes=FALSE )  
    if ( analysisType=="TS" | analysisType=="IO" ) {
        points( log10(XVal+1),log10(XFreq), type='l', col='darkgray' )  
    }
    points( log10(YsimZ0Val+1), log10(YsimZ0Freq), type='l', col='green' )
    points( log10(YsimVal_1S+1), log10(YsimFreq_1S), type='l', col='red' )  
    points( log10(YsimVal_2S+1), log10(YsimFreq_2S), type='l', col='blue' )

    if ( analysisType=="TS" | analysisType=="IO" ) {
        legend( 1, log10(max(YFreq)+1), 
            c('Actual data (ChIP)','Actual data (Control)','Sim:N','Sim:N+S1','Sim:N+S1+S2'),
            col=c('black','darkgray','green','red','blue'),lty=c(1,1,1,1,1),bty='n')
    } else {
        legend( 1, log10(max(YFreq)+1), 
            c('Actual data','Sim:N','Sim:N+S1','Sim:N+S1+S2'),
            col=c('black','green','red','blue'),lty=c(1,1,1,1),bty='n')
    }
    abline(v=log10(0+1),lty=2,col='gray')
    abline(v=log10(1+1),lty=2,col='gray')
    abline(v=log10(2+1),lty=2,col='gray')
    axis(1,0:6,10^c(0:6)-1)
    axis(2,0:10,10^c(0:10))

    #return(list(Y_sim_val = Y_sim_val, Y_sim_freq = Y_sim_freq))
}
