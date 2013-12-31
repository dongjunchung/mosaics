
.reportSummary <- function( summaryFile, resultList,
    chipFile=NULL, chipFileFormat=NULL, 
    controlFile=NULL, controlFileFormat=NULL, 
    binfileDir=NULL, 
    peakFile=NULL, peakFileFormat=NULL, 
    byChr=FALSE, FDR=0.05, fragLen=200, binSize=fragLen, capping=0, 
    analysisType="IO", d=0.25, 
    signalModel="BIC", maxgap=fragLen, minsize=50, thres=10 ) {
    
    cat( "MOSAiCS: Summary of model fitting and peak calling\n", file=summaryFile )
    cat( "\n", file=summaryFile, append=TRUE )
    cat( "------------------------------------------------------------\n", 
        file=summaryFile, append=TRUE )
    cat( "Input/output file settings\n", file=summaryFile, append=TRUE )
    cat( "------------------------------------------------------------\n", 
        file=summaryFile, append=TRUE )
    cat( "\n", file=summaryFile, append=TRUE )
    
    cat( "Name of aligned read file (ChIP): ",chipFile,"\n", sep="",
        file=summaryFile, append=TRUE )
    cat( "Aligned read file format (ChIP):", chipFileFormat, "\n",
        file=summaryFile, append=TRUE )
            
    cat( "\n", file=summaryFile, append=TRUE )
    
    cat( "Name of aligned read file (control): ",controlFile,"\n", sep="",
        file=summaryFile, append=TRUE )
    cat( "Aligned read file format (control):", controlFileFormat, "\n",
        file=summaryFile, append=TRUE )
                 
    cat( "\n", file=summaryFile, append=TRUE )
    
    for ( ff in 1:length(peakFileFormat) ) {   
        cat( "Name of peak result file: ",peakFile[ff],"\n", sep="",
            file=summaryFile, append=TRUE )
        cat( "Peak result file format:", peakFileFormat[ff], "\n",
            file=summaryFile, append=TRUE )
        
        cat( "\n", file=summaryFile, append=TRUE )
    }
               
    #cat( "\n", file=summaryFile, append=TRUE )    
    
    cat( "------------------------------------------------------------\n", 
        file=summaryFile, append=TRUE )
    cat( "Parameter settings\n", file=summaryFile, append=TRUE )
    cat( "------------------------------------------------------------\n", 
        file=summaryFile, append=TRUE )
    cat( "\n", file=summaryFile, append=TRUE )
    
    if ( byChr ) {    
        cat( "Genome-wide or chromosome-wise analysis? Chromosome-wise analysis\n",
            file=summaryFile, append=TRUE )
    } else  {    
        cat( "Genome-wide or chromosome-wise analysis? Genome-wide analysis\n",
            file=summaryFile, append=TRUE )
    }      
    
    if ( analysisType=="OS" ) {
        analysisTypeOut <- "One-sample analysis"
    } else if ( analysisType=="TS" ) {
        analysisTypeOut <- "Two-sample analysis (with mappability & GC content)"
    } else if ( analysisType=="IO" ) {
        analysisTypeOut <- "Two-sample analysis (Input only)"
    }
    
    if ( signalModel=="BIC" ) {
        signalModelOut <- "Automatic signal model selection based on BIC"
    } else if ( signalModel=="1S" ) {
        signalModelOut <- "One-signal-component model"
    } else if ( signalModel=="2S" ) {
        signalModelOut <- "Two-signal-component model"
    }
      
    cat( "False discovery rate (FDR):", FDR, "\n", file=summaryFile, append=TRUE )
    cat( "Fragment length:", fragLen, "\n", file=summaryFile, append=TRUE )
    cat( "Bin size:", binSize, "\n", file=summaryFile, append=TRUE )
    if ( capping > 0 ) {
        cat( "Maximum number of reads allowed in each nucleotide:", capping, "\n",
            file=summaryFile, append=TRUE )
    }
    cat( "Analysis type:", analysisTypeOut, "\n", file=summaryFile, append=TRUE )
    cat( "d:", d, "\n", file=summaryFile, append=TRUE )
    cat( "Signal model:", signalModelOut, "\n", file=summaryFile, append=TRUE )
    cat( "maxgap:", maxgap, "\n", file=summaryFile, append=TRUE )
    cat( "minsize:", minsize, "\n", file=summaryFile, append=TRUE )
    cat( "thres:", thres, "\n", file=summaryFile, append=TRUE )
       
    cat( "\n", file=summaryFile, append=TRUE )
    cat( "------------------------------------------------------------\n", 
        file=summaryFile, append=TRUE )
    cat( "Peak calling summary\n", file=summaryFile, append=TRUE )
    cat( "------------------------------------------------------------\n", 
        file=summaryFile, append=TRUE )
    cat( "\n", file=summaryFile, append=TRUE )
    
    outFormat <- data.frame( 
        resultList$chrID, resultList$n_peaks, 
        resultList$peak_width, resultList$opt_sig_model,
        stringsAsFactors=FALSE )
    colnames(outFormat) <- c( "chrID", "# peaks", 
        "Median peak width", "Optimal/specified signal model" )
      
    cat( as.character(colnames(outFormat)), file=summaryFile, sep="\t", append=TRUE )
    cat( "\n", file=summaryFile, append=TRUE )
    cat( rep("-----",3), file=summaryFile, sep="\t", append=TRUE )
    cat( "\t\t\t-----", file=summaryFile, sep="\t", append=TRUE )
    cat( "\n", file=summaryFile, append=TRUE )
    
    # peak list
     
    for ( i in 1:nrow(outFormat) )
    {
        cat( as.character(outFormat[i,])[1:3], file=summaryFile, sep="\t", append=TRUE )
        cat( "\t\t", as.character(outFormat[i,])[4], 
            file=summaryFile, sep="\t", append=TRUE )
        cat( "\n", file=summaryFile, append=TRUE )
    }
}
