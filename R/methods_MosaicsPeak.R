
# generic methods for "MosaicsPeak" class

setMethod(
    f="show",
    signature="MosaicsPeak",
    definition=function( object ) {  
              
        peakParam <- object@peakParam
        
        analysisType <- peakParam@analysisType
        signalModel <- peakParam@signalModel
        FDR <- peakParam@FDR
        maxgap <- peakParam@maxgap
        minsize <- peakParam@minsize
        thres <- peakParam@thres
        
        cat( "Summary: MOSAiCS peak calling (class: MosaicsPeak)\n" )
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
            }
        )        
        cat( "setting: FDR = ", FDR, ", maxgap = ", maxgap,
            ", minsize = ", minsize, ", thres = ", thres, "\n", sep="" )
        if ( length( object@peakList$peakStart ) > 0 ) {
            cat( "# of peaks = ", nrow(object@peakList), "\n", sep="" )
            cat( "median peak width = ", median(object@peakList$peakSize), "\n", sep="" )
            cat( "empirical FDR = ", round(empFDR(object)*10000)/10000, "\n", sep="" )
        } else {
            # exception handling (no peak case)
                    
            cat( "(no peak detected)\n" )
        }
        cat( "--------------------------------------------------\n" )
    }
)

setMethod(
    f="print",
    signature="MosaicsPeak",
    definition=function( x ) {    
        peakList <- x@peakList
        
        if ( nrow(peakList) > 0 ) {
            printForm <- x@peakList
            return(printForm)
        } else {
            # exception handling (no peak case)
            
            message( "Info: no peak detected. Nothing can be returned." )           
        }
    }
)

setMethod(
    f="export",
    signature="MosaicsPeak",
    definition=function( object, type=NA, filename=NA ) {
        # error treatment: check invalid type
        
        if ( is.na(type) )
        {
            message( "Info: 'type' is not specified by the user." )
            message( "Info: 'type' is specified as 'txt' instead." )
            type <- "txt"        
        }
        
        allType <- c("txt","gff","bed")
        invalidType <- TRUE
        for ( i in 1:length(type) )
        {
            if ( length(which(!is.na(match(type[i],allType))))==0 )
            {
                invalidType <- FALSE
            }
        }
        if ( !invalidType )
        {
            message( "Info: 'type' incorrect." )
            message( "Info: 'type' is specified as 'txt' instead." )
            type <- "txt"
        }
        
        # error treatment: 'filename' not specified
        
        if ( is.na(filename) )
        {
            message( "Info: 'filename' is not specified by the user." )
            message( "Info: 'filename' is specified as 'peakList' instead." )
            filename <- paste("peakList.",type,sep="")
        }
        
        # export peak lists
        
        peakList <- object@peakList
        
        if ( nrow(object@peakList) > 0 ) {        
                message( "Info: exporting the peak list..." )
                switch( type,
                    "gff" = {
                        .exportGFF( peakList=peakList, filename=filename )
                    },
                    "bed" = {
                        .exportBED( peakList=peakList, filename=filename )
                    },
                    "txt" = {
                        .exportTXT( peakList=peakList, filename=filename )
                    }
                )
            } else {
                # exception handling (no peak case)
                
                message( "Info: no peak identifed. Nothing exported." )
            }
    }
)

setMethod(
    f="empFDR",
    signature="MosaicsPeak",
    definition=function( object ) {
        return(object@empFDR)
    }
)

setMethod(
    f="bdBin",
    signature="MosaicsPeak",
    definition=function( object ) {
        return(object@bdBin)
    }
)
