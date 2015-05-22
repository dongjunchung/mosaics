
# MOSAiCS peak calling

setMethod(
    f="mosaicsPeak",
    signature="MosaicsFit",
    definition=function( object, signalModel="2S", FDR=0.05, 
        binsize=NA, maxgap=200, minsize=50, thres=10 )
    {
        #mosaicsEst <- object@mosaicsEst
        #tagCount <- object@tagCount
        analysisType <- object@mosaicsEst@analysisType
        
        # error treatment: invalid signal model specification
        
        if ( signalModel!="1S" && signalModel!="2S" )
        {
            stop( "Invalid signal model: should be either '1S' or '2S'!" )
        }
        
        # fit marginal distribution & posterior distribution
        
        switch( signalModel,
            "1S" = {
                # one signal model
                
                pi0 <- object@mosaicsEst@pi0
                
                message( "Info: use one-signal-component model." )
                message( "Info: calculating posterior probabilities..." )
                fitMD <- .margDist_1S( mosaicsEst=object@mosaicsEst, 
                    tagCount=object@tagCount, pNfit=object@mosaicsEst@pNfit, k=3 )
                fitPH <- .getPH_1S( margDensity=fitMD, pi0=pi0 )
            },
            "2S" = {
                # two signal model
                
                pi0 <- object@mosaicsEst@pi0
                p1 <- object@mosaicsEst@p1
                
                message( "Info: use two-signal-component model." )
                message( "Info: calculating posterior probabilities..." )
                fitMD <- .margDist_2S( mosaicsEst=object@mosaicsEst, 
                    tagCount=object@tagCount, pNfit=object@mosaicsEst@pNfit, k=3 )
                fitPH <- .getPH_2S( margDensity=fitMD, pi0=pi0, p1=p1 )        
            }
        )
        
        rm(fitMD)
        gc()
        
        # peak calling
        
        message( "Info: calling peaks..." )
        switch( analysisType,
            OS = {
                dataSet <- list( chrID=object@chrID, coord=object@coord, Y=object@tagCount,
                    M=object@mappability, GC=object@gcContent )
                nRatio <- 1
            },
            TS = {
                dataSet <- list( chrID=object@chrID, coord=object@coord, Y=object@tagCount,
                    X=object@input, M=object@mappability, GC=object@gcContent )
                nRatio <- object@seqDepth[1] / object@seqDepth[2]
            },
            IO = {
                dataSet <- list( chrID=object@chrID, coord=object@coord, Y=object@tagCount, 
                    X=object@input )
                nRatio <- object@seqDepth[1] / object@seqDepth[2]
            }        
        )
        fitPeak <- .peakCall( postProb=fitPH, dataSet=dataSet,
            FDR=FDR, binsize=binsize, maxgap=maxgap, minsize=minsize, thres=thres, 
            analysisType=analysisType, nRatio=nRatio )
        peakSet <- fitPeak$peakSet
        
        #rm( fitPH, dataSet )
        rm( dataSet )
        gc()
        
        message( "Info: done!" )
        
        # construct "MosaicsPeak" class object
                
        peakParam <- new( "MosaicsPeakParam",
            analysisType=(object@mosaicsEst)@analysisType, signalModel=signalModel, 
            FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres, decoding="posterior" )
        
        if ( !is.null(peakSet) ) {
            peakList <- fitPeak$peakSet 
        } else {
            # exception handling (no peak case)
            
            peakList <- data.frame()
        }
        
        postProb <- data.frame( object@chrID, object@coord, fitPH$pH0,
        	stringsAsFactors=FALSE )
        colnames(postProb) <- c( "chrID", "coord", "PostProb" )
      
        tagDataEmpty <- new( "TagData", 
            read=list(), coverage=list() )   
        
        new( "MosaicsPeak",         
            peakList=peakList, 
            chrID=object@chrID, coord=object@coord, 
            tagCount=object@tagCount, input=object@input, 
            mappability=object@mappability, gcContent=object@gcContent,
            peakParam=peakParam, bdBin=fitPeak$bdBin, postProb=postProb, empFDR=fitPeak$empFDR,
            tagLoaded=FALSE, tagData=tagDataEmpty, seqDepth=object@seqDepth )
    }
)
