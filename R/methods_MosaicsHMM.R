
# generic methods for "MosaicsFit" class

setMethod(
    f="show",
    signature="MosaicsHMM",
    definition=function( object ) {    
              
        peakParam <- object@peakParam
        
        analysisType <- peakParam@analysisType
        signalModel <- peakParam@signalModel
        FDR <- peakParam@FDR
        maxgap <- peakParam@maxgap
        minsize <- peakParam@minsize
        thres <- peakParam@thres
        
        cat( "Summary: MOSAiCS-HMM model fitting (class: MosaicsHMM)\n" )
        cat( "--------------------------------------------------\n" )
        cat( "final model: " )
        switch( analysisType,
            "OS" = {
                cat( "one-sample analysis " )
            },
            "TS" = {
                cat( "two-sample analysis (with M & GC) " )
            },
            "IO" = {
                cat( "two-sample analysis (input only) " )
            }
        )
        switch( signalModel,
            "1S" = {
                cat( "with one signal component\n" )
            },
            "2S" = {
                cat( "with two signal components\n" )
            },
            "2Sr" = {
                cat( "with two robust signal components\n" )
            }
        )        
        cat( "setting for initialization:\n" )
        cat( "FDR = ", FDR, ", maxgap = ", maxgap,
            ", minsize = ", minsize, ", thres = ", thres, "\n", sep="" )
        cat( "--------------------------------------------------\n" )
    }
)
