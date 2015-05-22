
# generic methods for "MosaicsFit" class

setMethod(
    f="show",
    signature="MosaicsFit",
    definition=function( object ) {      
        mosaicsParam <- object@mosaicsParam    
        
        bic1S <- object@bic1S
        bic2S <- object@bic2S
        
        analysisType <- object@mosaicsEst@analysisType         
        
        cat( "Summary: MOSAiCS model fitting (class: MosaicsFit)\n" )
        cat( "--------------------------------------------------\n" )
        switch( analysisType,
            OS = {        
                k <- mosaicsParam@k      
                meanThres <- mosaicsParam@meanThres  
                
                cat( "analysis type: one-sample analysis\n" )
                cat( "parameters used: k = ", k, ", meanThres = ", meanThres, "\n", sep="" )
                cat( "BIC of one-signal-component model = ", bic1S, "\n", sep="" )
                cat( "BIC of two-signal-component model = ", bic2S, "\n", sep="" )
            },
            TS = {
                k <- mosaicsParam@k      
                meanThres <- mosaicsParam@meanThres  
                s <- mosaicsParam@s      
                d <- mosaicsParam@d      
            
                cat( "analysis type: two-sample analysis (with mappability & GC content)\n" )
                cat( "parameters used: k = ", k, ", meanThres = ", meanThres,
                    ", s = ", s, ", d = ", d, "\n", sep="" )
                cat( "BIC of one-signal-component model = ", bic1S, "\n", sep="" )
                cat( "BIC of two-signal-component model = ", bic2S, "\n", sep="" )                
            },
            IO = {
                k <- mosaicsParam@k      
                d <- mosaicsParam@d
            
                cat( "analysis type: two-sample analysis (Input only)\n" )
                cat( "parameters used: k = ", k, ", d = ", d, "\n", sep="" )
                cat( "BIC of one-signal-component model = ", bic1S, "\n", sep="" )
                cat( "BIC of two-signal-component model = ", bic2S, "\n", sep="" )                
            }
        )
        
        cat( "--------------------------------------------------\n" )
    }
)

setMethod(
    f="print",
    signature="MosaicsFit",
    definition=function( x ) {
        warning( "'print' method for 'MosaicsFit' class is not supported yet." )
    }
)

setMethod(
    f="plot",
    signature="MosaicsFit",
    definition=function( x ) {
        if ( x@mosaicsEst@analysisType=="OS" ) {
            .plotGOF( mosaicsEst=x@mosaicsEst, tagCount=x@tagCount, k=3 )
        } else {
            .plotGOF( mosaicsEst=x@mosaicsEst, tagCount=x@tagCount, input=x@input, k=3 )
        }
    }
)

setMethod(
    f="estimates",
    signature="MosaicsFit",
    definition=function( object ) {
        mosaicsEst <- object@mosaicsEst
        
        list(            
            pi0=mosaicsEst@pi0,
            a=mosaicsEst@a,
            betaEst=mosaicsEst@betaEst,
            b=mosaicsEst@b,
            c=mosaicsEst@c,
            p1=mosaicsEst@p1,
            b1=mosaicsEst@b1,
            c1=mosaicsEst@c1,
            b2=mosaicsEst@b2,
            c2=mosaicsEst@c2
        )            
    }
)

setMethod(
    f="seqDepth",
    signature="MosaicsFit",
    definition=function( object ) {
      return(object@seqDepth)
    }
)
