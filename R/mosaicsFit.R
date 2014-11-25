
# fit MOSAiCS

setMethod(
    f="mosaicsFit",
    signature="BinData",
    definition=function( object, analysisType="automatic", bgEst="rMOM",
        k=3, meanThres=NA, s=2, d=0.25, trans="power", truncProb=0.999, 
        parallel=FALSE, nCore=8 )
    {
        # Note: users can tune parameters only regarding MOSAiCS model fitting.
        # Note: tuning adaptive griding parameters is not supported yet.
    
        # check options: parallel computing (optional)
        
        if ( parallel == TRUE ) {
            message( "Use 'parallel' package for parallel computing." )        
            if ( length(find.package('parallel',quiet=TRUE)) == 0 ) {
                stop( "Please install 'parallel' package!" )
            }
        }
        
        if ( analysisType == "automatic" ) {
                # if "analysisType" is not specified, "analysisType" is determined by dataset
                
                if ( length(object@input)==0 ) {
                        # we don't have input: we should have both of M and GC
                        
                        if ( length(object@mappability)>0 & length(object@gcContent)>0 ) {
                            #cat( "Info: one-sample analysis.\n" )
                                analysisType <- "OS"
                        } else {
                            stop( "any of control data, mappability, or GC content does not exist. Cannot proceed!\n" )
                        } 
                } else {
                        # we have input: TS if we have both M & GC; IO otherwise. 
                        
                        if ( length(object@mappability)>0 & length(object@gcContent)>0 ) {
                            #cat( "Info: two-sample analysis (with mappability & GC content).\n" )
                                analysisType <- "TS"
                        } else {
                            #cat( "Info: two-sample analysis (Input only).\n" )
                                analysisType <- "IO"
                        }
                }
        } else {
                # if "analysisType" is specified, check its validity
                                
                # error treatment: Input-only analysis is impossible if input data does not exist
                
                if ( analysisType=="IO" & length(object@input)==0 )
                {
                    message( "Info: two-sample analysis (Input only)." )
                    stop( "control data does not exist. Cannot proceed!\n" )
                }
                
                # error treatment: If M or GC does not exist, TS analysis is not available
                
                if ( analysisType=="OS" & ( length(object@mappability)==0 | length(object@gcContent)==0 ) )
                {
                    message( "Info: one-sample analysis." )
                    stop( "mappability or GC content does not exist. Cannot proceed!\n" )
                }
                
                # error treatment: If input data does not exist, TS analysis is not available
                
                if ( analysisType=="TS" & length(object@input)==0 )
                {
                    message( "Info: two-sample analysis (with mappability & GC content)." )
                    message( "Info: control data does not exist." )
                    message( "Info: one-sample analysis will be implemented instead." )
                    analysisType <- "OS"
                }
                
                # error treatment: If M or GC does not exist, TS analysis is not available
                
                if ( analysisType=="TS" & ( length(object@mappability)==0 | length(object@gcContent)==0 ) )
                {
                    message( "Info: two-sample analysis (with mappability & GC content)." )
                    message( "Info: mappability or GC content data does not exist." )
                    message( "Info: two-sample analysis (Input only) will be implemented instead." )
                    analysisType <- "IO"
                }
        }
        
        # check validity of "bgEst"
        
        if ( bgEst == "automatic" ) {
            message( "Info: background estimation method is determined based on data." )
            
            Y_freq <- table( object@tagCount )
                        
            if ( sum(Y_freq[ as.numeric(names(Y_freq))<=2 ]) / sum(Y_freq) > 0.5 ) {
                message( "Info: background estimation based on bins with low tag counts." )
                bgEst <- "matchLow"
            } else {
                message( "Info: background estimation based on robust method of moment" )
                bgEst <- "rMOM"
            }
        } else if ( bgEst == "matchLow" ) {
            message( "Info: background estimation based on bins with low tag counts." )
        } else if ( bgEst == "rMOM" ) {
            message( "Info: background estimation based on robust method of moment." )
        } else {
            stop( "Incorrect specification for 'bgEst'! 'bgEst' should be one of 'matchLow', 'rMOM', or 'automatic'!" )
        }
        
        # default meanThres for each of "OS" & "TS"
        
        if ( is.na(meanThres) )
        {
            switch( analysisType,
                OS = {
                    meanThres <- 0
                },
                TS = {
                    meanThres <- 1
                },
                IO = {
                    meanThres <- NA     # meanThres is not used for analysisType=="IO"
                }
            )
        }
        
        # MOSAiCS model fit
        
        switch( analysisType,
            OS = {
                # one-sample analysis
                
                message( "Info: one-sample analysis." )
                fit <- .mosaicsFit_OS( object, bgEst=bgEst, k=k, meanThres=meanThres,
                    parallel=parallel, nCore=nCore )
            },
            TS = {
                # two-sample analysis (with M & GC)
                
                message( "Info: two-sample analysis (with mappability & GC content)." )
                fit <- .mosaicsFit_TS( object, bgEst=bgEst, 
                    k=k, meanThres=meanThres, s=s, d=d,
                    parallel=parallel, nCore=nCore )            
            },
            IO = {
                # two-sample analysis (Input only)
                
                message( "Info: two-sample analysis (Input only)." )
                fit <- .mosaicsFit_IO( object, bgEst=bgEst, 
                    k=k, d=d, trans=trans, truncProb=truncProb,
                    parallel=parallel, nCore=nCore )    
            }
        )
        
        message( "Info: done!" )
        
        return(fit)
    }
)

# MOSAiCS one-sample analysis

.mosaicsFit_OS <- function( binData, bgEst, k=3, meanThres=0, 
    parallel=FALSE, nCore=8 )
{        
    message( "Info: use adaptive griding." )
    message( "Info: fitting background model..." )    
    fitParam <- .adapGridMosaicsZ0_OS(
        Y=binData@tagCount, M=binData@mappability, GC=binData@gcContent, 
        bgEst=bgEst, min_n_MGC=50, grids_MGC=c(0.01,0.02,0.04,0.10,0.20,0.50),
        parallel=parallel, nCore=nCore )
    fitZ0 <- .rlmFit_OS( parEst=fitParam, mean_thres=meanThres, bgEst=bgEst,
        Y=binData@tagCount, M=binData@mappability, GC=binData@gcContent )
    pNfit <- .calcPN( Y=binData@tagCount, k=k, a=fitZ0$a, mu_est=fitZ0$muEst ) 
    
    rm( fitParam )
    gc()
    message( "Info: done!" )
    
    Y_bd_all <- .calcYbdAll( fitZ0, k=k )
    message( "Info: fitting one-signal-component model..." )
    fitZ1_1S <- .mosaicsZ1_1S( fitZ0, Y=binData@tagCount, 
        pNfit=pNfit, Y_bd_all=Y_bd_all, k=k )
    message( "Info: fitting two-signal-component model..." )
    fitZ1_2S <- .mosaicsZ1_2S( fitZ0, Y=binData@tagCount, 
        pNfit=pNfit, Y_bd_all=Y_bd_all, k=k )
    
    #message( "Info: calculating BIC of fitted models..." )
    #fitBIC_1S <- .calcModelBIC( fitZ1=fitZ1_1S, Y=binData@tagCount, 
    #    pNfit=pNfit, k=k, model="1S", type="BIC", npar=9 )
    #fitBIC_2S <- .calcModelBIC( fitZ1=fitZ1_2S, Y=binData@tagCount, 
    #    pNfit=pNfit, k=k, model="2S", type="BIC", npar=12 )
    
    mosaicsEst <- new( "MosaicsFitEst",
        pi0=fitZ0$pi0, a=fitZ0$a,
        betaEst=fitZ0$betaEst, muEst=fitZ0$muEst, pNfit=pNfit,
        b=fitZ1_1S$b, c=fitZ1_1S$c,
        p1=fitZ1_2S$p1, b1=fitZ1_2S$b1, c1=fitZ1_2S$c1, b2=fitZ1_2S$b2, c2=fitZ1_2S$c2,
        analysisType="OS" )
    
    rm( fitZ0, fitZ1_1S, fitZ1_2S )
    gc()
    
    message( "Info: calculating BIC of fitted models..." )
    loglik_1S <- .logLik( mosaicsEst=mosaicsEst, tagCount=binData@tagCount, 
        pNfit=pNfit, k=k, signalModel="1S" )
    loglik_2S <- .logLik( mosaicsEst=mosaicsEst, tagCount=binData@tagCount, 
        pNfit=pNfit, k=k, signalModel="2S" )
    fitBIC_1S <- .calcModelBIC( 
    	loglik=loglik_1S, n=length(binData@tagCount), nChr=1,
    	method="mosaics", analysisType="OS", signalModel="1S", type="BIC" )
    fitBIC_2S <- .calcModelBIC( 
    	loglik=loglik_2S, n=length(binData@tagCount), nChr=1,
    	method="mosaics", analysisType="OS", signalModel="2S", type="BIC" )
       
    mosaicsParam <- new( "MosaicsFitParam", k=k, meanThres=meanThres )
    
    new( "MosaicsFit",
        mosaicsEst=mosaicsEst, mosaicsParam=mosaicsParam,
        chrID=binData@chrID, coord=binData@coord, tagCount=binData@tagCount, 
        mappability=binData@mappability, gcContent=binData@gcContent,
        bic1S=fitBIC_1S, bic2S=fitBIC_2S )
}

# MOSAiCS two-sample analysis (with M & GC)

.mosaicsFit_TS <- function( binData, bgEst, k=3, meanThres=1, s=2, d=0.25, 
    parallel=FALSE, nCore=8 )
{    
    message( "Info: use adaptive griding." )
    message( "Info: fitting background model..." )     
    
    # warning if there are insufficient # of bins with 0, 1, 2 counts in control sample
    # -> input-only analysis is preferred.
    
    X_freq <- table(binData@input)
    if( sum(X_freq[ as.numeric(names(X_freq))<=2 ])/sum(X_freq) < 0.5 ) {
        message( paste("Info: insufficient # of bins with counts <=", s, "in control sample.") )
        message( "Info: Model fit can be unstable. Input-only analysis might be preferred." )
    }
    
    fitParam <- .adapGridMosaicsZ0_TS( 
        Y=binData@tagCount, M=binData@mappability, GC=binData@gcContent, X=binData@input,
        bgEst=bgEst, min_n_MGC=50, grids_MGC=c(0.01,0.02,0.04,0.10,0.20,0.50), min_n_X=200,
        parallel=parallel, nCore=nCore )
    fitZ0 <- .rlmFit_TS( parEst=fitParam, mean_thres=meanThres, s=s, d=d, bgEst=bgEst,
        Y=binData@tagCount, M=binData@mappability, GC=binData@gcContent, X=binData@input )
    pNfit <- .calcPN( Y=binData@tagCount, k=k, a=fitZ0$a, mu_est=fitZ0$muEst ) 
    
    rm( fitParam )
    gc()
    message( "Info: done!" )
    
    Y_bd_all <- .calcYbdAll( fitZ0, k=k )
    message( "Info: fitting one-signal-component model..." )
    fitZ1_1S <- .mosaicsZ1_1S( fitZ0, Y=binData@tagCount, 
        pNfit=pNfit, Y_bd_all=Y_bd_all, k=k )
    message( "Info: fitting two-signal-component model..." )
    fitZ1_2S <- .mosaicsZ1_2S( fitZ0, Y=binData@tagCount, 
        pNfit=pNfit, Y_bd_all=Y_bd_all, k=k )
    
    #message( "Info: calculating BIC of fitted models..." )
    #fitBIC_1S <- .calcModelBIC( fitZ1=fitZ1_1S, Y=binData@tagCount, 
    #    pNfit=pNfit, k=k, model="1S", type="BIC", npar=11 )
    #fitBIC_2S <- .calcModelBIC( fitZ1=fitZ1_2S, Y=binData@tagCount, 
    #    pNfit=pNfit, k=k, model="2S", type="BIC", npar=14 )
    
    mosaicsEst <- new( "MosaicsFitEst",
        pi0=fitZ0$pi0, a=fitZ0$a,
        betaEst=fitZ0$betaEst, muEst=fitZ0$muEst, pNfit=pNfit,
        b=fitZ1_1S$b, c=fitZ1_1S$c,
        p1=fitZ1_2S$p1, b1=fitZ1_2S$b1, c1=fitZ1_2S$c1, b2=fitZ1_2S$b2, c2=fitZ1_2S$c2,
        analysisType="TS" )
    
    rm( fitZ0, fitZ1_1S, fitZ1_2S )
    gc()
    
    message( "Info: calculating BIC of fitted models..." )
    loglik_1S <- .logLik( mosaicsEst=mosaicsEst, tagCount=binData@tagCount, 
        pNfit=pNfit, k=k, signalModel="1S" )
    loglik_2S <- .logLik( mosaicsEst=mosaicsEst, tagCount=binData@tagCount, 
        pNfit=pNfit, k=k, signalModel="2S" )
    fitBIC_1S <- .calcModelBIC( 
    	loglik=loglik_1S, n=length(binData@tagCount), nChr=1,
    	method="mosaics", analysisType="TS", signalModel="1S", type="BIC" )
    fitBIC_2S <- .calcModelBIC( 
    	loglik=loglik_2S, n=length(binData@tagCount), nChr=1,
    	method="mosaics", analysisType="TS", signalModel="2S", type="BIC" )
        
    mosaicsParam <- new( "MosaicsFitParam", k=k, meanThres=meanThres, s=s, d=d )
    
    new( "MosaicsFit",
        mosaicsEst=mosaicsEst, mosaicsParam=mosaicsParam,
        chrID=binData@chrID, coord=binData@coord, 
        tagCount=binData@tagCount, input=binData@input, 
        mappability=binData@mappability, gcContent=binData@gcContent,
        bic1S=fitBIC_1S, bic2S=fitBIC_2S )
}

# MOSAiCS two-sample analysis (Input only)

.mosaicsFit_IO <- function( binData, bgEst, k=3, d=0.25, trans="log",
	truncProb=0.999, parallel=FALSE, nCore=8 )
{    
    message( "Info: use adaptive griding." )
    message( "Info: fitting background model..." )
    
    inputTrunc <- quantile( binData@input, truncProb )
    
    #fitParam <- .adapGridMosaicsZ0_IO( Y=binData@tagCount, X=binData@input, 
    #    min_n_X=50 )
    fitParam <- .adapGridMosaicsZ0_IO( Y=binData@tagCount, X=binData@input, 
        bgEst=bgEst, inputTrunc=inputTrunc, min_n_X=50,
        parallel=parallel, nCore=nCore )
    fitZ0 <- .rlmFit_IO( parEst=fitParam, d=d, trans=trans, bgEst=bgEst,
		Y=binData@tagCount, X=binData@input, inputTrunc=inputTrunc )
    pNfit <- .calcPN( Y=binData@tagCount, k=k, a=fitZ0$a, mu_est=fitZ0$muEst ) 
    
    rm( fitParam )
    gc()
    message( "Info: done!" )
    
    Y_bd_all <- .calcYbdAll( fitZ0, k=k )
    message( "Info: fitting one-signal-component model..." )
    fitZ1_1S <- .mosaicsZ1_1S( fitZ0, Y=binData@tagCount, 
        pNfit=pNfit, Y_bd_all=Y_bd_all, k=k )
    message( "Info: fitting two-signal-component model..." )
    fitZ1_2S <- .mosaicsZ1_2S( fitZ0, Y=binData@tagCount, 
        pNfit=pNfit, Y_bd_all=Y_bd_all, k=k )
    
    #message( "Info: calculating BIC of fitted models..." )
    #fitBIC_1S <- .calcModelBIC( fitZ1=fitZ1_1S, Y=binData@tagCount, 
    #    pNfit=pNfit, k=k, model="1S", type="BIC", npar=6 )
    #fitBIC_2S <- .calcModelBIC( fitZ1=fitZ1_2S, Y=binData@tagCount, 
    #    pNfit=pNfit, k=k, model="2S", type="BIC", npar=9 )
    
    mosaicsEst <- new( "MosaicsFitEst",
        pi0=fitZ0$pi0, a=fitZ0$a,
        betaEst=fitZ0$betaEst, muEst=fitZ0$muEst, pNfit=pNfit,
        b=fitZ1_1S$b, c=fitZ1_1S$c,
        p1=fitZ1_2S$p1, b1=fitZ1_2S$b1, c1=fitZ1_2S$c1, b2=fitZ1_2S$b2, c2=fitZ1_2S$c2,
        inputTrunc=inputTrunc, analysisType="IO" )
    
    rm( fitZ0, fitZ1_1S, fitZ1_2S )
    gc()
    
    message( "Info: calculating BIC of fitted models..." )
    loglik_1S <- .logLik( mosaicsEst=mosaicsEst, tagCount=binData@tagCount, 
        pNfit=pNfit, k=k, signalModel="1S" )
    loglik_2S <- .logLik( mosaicsEst=mosaicsEst, tagCount=binData@tagCount, 
        pNfit=pNfit, k=k, signalModel="2S" )
    fitBIC_1S <- .calcModelBIC( 
    	loglik=loglik_1S, n=length(binData@tagCount), nChr=1,
    	method="mosaics", analysisType="IO", signalModel="1S", type="BIC" )
    fitBIC_2S <- .calcModelBIC( 
    	loglik=loglik_2S, n=length(binData@tagCount), nChr=1, 
    	method="mosaics", analysisType="IO", signalModel="2S", type="BIC" )
        
    mosaicsParam <- new( "MosaicsFitParam", k=k, d=d )
    
    new( "MosaicsFit",
        mosaicsEst=mosaicsEst, mosaicsParam=mosaicsParam,
        chrID=binData@chrID, coord=binData@coord, 
        tagCount=binData@tagCount, input=binData@input,
        bic1S=fitBIC_1S, bic2S=fitBIC_2S )
}
