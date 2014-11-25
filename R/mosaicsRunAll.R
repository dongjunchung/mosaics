mosaicsRunAll <- function( 
    chipFile=NULL, chipFileFormat=NULL, 
    controlFile=NULL, controlFileFormat=NULL, 
    binfileDir=NULL, 
    peakFile=NULL, peakFileFormat=NULL,
    reportSummary=FALSE, summaryFile=NULL, 
    reportExploratory=FALSE, exploratoryFile=NULL, 
    reportGOF=FALSE, gofFile=NULL, 
    PET=FALSE, byChr=FALSE, useChrfile=FALSE, chrfile=NULL, excludeChr=NULL, 
    FDR=0.05, fragLen=200, binSize=200, capping=0, bgEst="rMOM", d=0.25, 
    signalModel="BIC", maxgap=200, minsize=50, thres=10, parallel=FALSE, nCore=8 ) {
    
    analysisType <- "IO"
    
    # check options: input & output (required)
    
    if ( is.null(chipFile) ) { stop( "Please specify 'chipFile'!" ) }
    if ( is.null(chipFileFormat) ) { stop( "Please specify 'chipFileFormat'!" ) }
    
    if ( is.null(controlFile) ) { stop( "Please specify 'controlFile'!" ) }
    if ( is.null(controlFileFormat) ) { stop( "Please specify 'controlFileFormat'!" ) }
    
    if ( is.null(peakFile) ) { stop( "Please specify 'peakFile'!" ) }
    if ( is.null(peakFileFormat) ) { stop( "Please specify 'peakFileFormat'!" ) }
    
    if ( is.null(binfileDir) ) { stop( "Please specify 'binfileDir'!" ) }
    
    # check options: peak list
    
    if ( length(peakFile) != length(peakFileFormat) ) {
        stop( "Lengths of 'peakFileName' and 'peakFileFormat' should be same!" )
    }
    
    # check options: reports (optional)
    
    if ( reportSummary ) {
        if ( is.null(summaryFile) ) {
            stop( "Please specify 'summaryFile'!" )
        }
    }
  
    if ( reportGOF ) {
        if ( is.null(gofFile) ) {
            stop( "Please specify 'gofFile'!" )
        }
    }
  
    if ( reportExploratory ) {
        if ( is.null(exploratoryFile) ) {
            stop( "Please specify 'exploratoryFile'!" )
        }
    }
    
    # check options: parallel computing (optional)
    
    if ( parallel == TRUE ) {
        message( "Use 'parallel' package for parallel computing." )        
        if ( length(find.package('parallel',quiet=TRUE)) == 0 ) {
            stop( "Please install 'parallel' package!" )
        }
    }
  
    # construction of bin-level files

    cat( "Info: constructing bin-level files...\n" )
    
    processSet <- list()
    processSet[[1]] <- c( chipFile, chipFileFormat )
    processSet[[2]] <- c( controlFile, controlFileFormat )
    
    
    if ( parallel == TRUE ) {
        # if "parallel" package exists, utilize parallel computing with "mclapply"
        
        mclapply( processSet, function(x) {
            constructBins( 
                infile = x[1], fileFormat = x[2], outfileLoc = binfileDir, 
                byChr = byChr, useChrfile = useChrfile, chrfile = chrfile, excludeChr = excludeChr,
                PET = PET, fragLen = fragLen, binSize = binSize, capping = capping )    
        }, mc.cores=nCore )
    } else {
        # otherwise, use usual "lapply"
        
        lapply( processSet, function(x) {
            constructBins( 
                infile = x[1], fileFormat = x[2], outfileLoc = binfileDir, 
                byChr = byChr, useChrfile = useChrfile, chrfile = chrfile, excludeChr = excludeChr,
                PET = PET, fragLen = fragLen, binSize = binSize, capping = capping )    
        } )
    }
    
    if ( byChr ) {
        ###############################################################
        #                                                             #
        #                chromosome-wise analysis                     #
        #                                                             #
        ###############################################################
    
        # read in bin-level files                           
        
        cat( "Info: analyzing bin-level files...\n" )      
          
        setwd( binfileDir )
        if ( PET == TRUE ) {
            list_chip <- list.files( path=binfileDir,
                paste(basename(chipFile),"_bin",binSize,"_.*.txt",sep="")  ) 
            list_control <- list.files( path=binfileDir,
                paste(basename(controlFile),"_bin",binSize,"_.*.txt",sep="")  )
        } else {
            list_chip <- list.files( path=binfileDir,
                paste(basename(chipFile),"_fragL",fragLen,"_bin",binSize,"_.*.txt",sep="")  ) 
            list_control <- list.files( path=binfileDir,
                paste(basename(controlFile),"_fragL",fragLen,"_bin",binSize,"_.*.txt",sep="")  )
        }
        
        # check list of chromosomes & analyze only chromosomes 
        # that bin-level files for both chip & control exist
        
        #chrID_chip <- unlist( lapply( strsplit( list_chip, paste("_",basename(chipFile),sep="") ), 
        #    function(x) x[1] ) )
        #chrID_control <- unlist( lapply( strsplit( list_control, paste("_",controlFile,sep="") ), 
        #    function(x) x[1] ) )        
        chrID_chip <- unlist( lapply( list_chip, function(x) {
            splitvec <- strsplit( x, "_" )[[1]]
            IDtxt <- splitvec[ length(splitvec) ]
            return( strsplit( IDtxt, ".txt" )[[1]][1] )
        } ) )
        chrID_control <- unlist( lapply( list_control, function(x) {
            splitvec <- strsplit( x, "_" )[[1]]
            IDtxt <- splitvec[ length(splitvec) ]
            return( strsplit( IDtxt, ".txt" )[[1]][1] )
        } ) )      
        index_chip <- which( !is.na( match( chrID_chip, chrID_control ) ) )
        index_control <- match( chrID_chip, chrID_control )
        index_list <- list()
        for ( i in 1:length(index_chip) ) {
            index_list[[i]] <- c( index_chip[i], index_control[i] )
        }
        
        # model fitting & peak calling
        # check whether rparallel is available. if so, use it.
        
        cat( "Info: fitting MOSAiCS model & call peaks...\n" )
      
        if ( length(index_chip) < nCore ) {
            nCore <- length(index_chip)
        }
        
        if ( parallel == TRUE ) {
            # if "parallel" package exists, utilize parallel computing with "mclapply"
            
            out <- mclapply( index_list, function(x) {    
                # read in bin-level file
                
                chip_file <- list_chip[ x[1] ]
                input_file <- list_control[ x[2] ]                
                bin <- readBins( 
                    type=c("chip","input"), fileName=c(chip_file,input_file),
                    parallel=parallel, nCore=nCore )
          
                # fit model
                
                fit <- mosaicsFit( bin, analysisType=analysisType, bgEst=bgEst, d=d,
                    parallel=parallel, nCore=nCore )    
                
                # call peaks
                
                if ( signalModel=="BIC" ) {
                    # if not specified, use BIC
                    
                    if ( fit@bic1S < fit@bic2S ) {
                        peak <- mosaicsPeak( fit, signalModel="1S", 
                            FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )      
                        opt_sig_model <- "One-signal-component model"
                    } else {        
                        peak <- mosaicsPeak( fit, signalModel="2S", 
                            FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )      
                        opt_sig_model <- "Two-signal-component model"
                    }
                } else {                                                  
                    peak <- mosaicsPeak( fit, signalModel=signalModel, 
                        FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )
                    if ( signalModel=="1S" ) {
                        opt_sig_model <- "One-signal-component model"
                    } else {
                        opt_sig_model <- "Two-signal-component model"
                    }
                }  
          
                # keep results
                
                peakPrint <- print(peak)
                
                return( list( chrID=as.character( chrID_chip[ x[1] ] ), 
                    bin=bin, fit=fit, peak=peak, peakPrint=peakPrint,
                    n_peaks=nrow(peak@peakList), peak_width=median(peak@peakList$peakSize),
                    opt_sig_model=opt_sig_model ) )
            }, mc.cores=nCore )
        } else {
            # otherwise, use usual "lapply"
            
            out <- lapply( index_list, function(x) {    
                # read in bin-level file
                
                chip_file <- list_chip[ x[1] ]
                input_file <- list_control[ x[2] ]                            
                bin <- readBins( 
                    type=c("chip","input"), fileName=c(chip_file,input_file),
                    parallel=parallel, nCore=nCore )
          
                # fit model
                
                fit <- mosaicsFit( bin, analysisType=analysisType, bgEst=bgEst, d=d,
                    parallel=parallel, nCore=nCore )    
                
                # call peaks
                
                if ( signalModel=="BIC" ) {
                    # if not specified, use BIC
                    
                    if ( fit@bic1S < fit@bic2S ) {
                        peak <- mosaicsPeak( fit, signalModel="1S", 
                            FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )      
                        opt_sig_model <- "One-signal-component model"
                    } else {        
                        peak <- mosaicsPeak( fit, signalModel="2S", 
                            FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )      
                        opt_sig_model <- "Two-signal-component model"
                    }
                } else {                                                  
                    peak <- mosaicsPeak( fit, signalModel=signalModel, 
                        FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )
                    if ( signalModel=="1S" ) {
                        opt_sig_model <- "One-signal-component model"
                    } else {
                        opt_sig_model <- "Two-signal-component model"
                    }
                }  
          
                # keep results
                
                peakPrint <- print(peak)
                
                return( list( chrID=as.character( chrID_chip[ x[1] ] ), 
                    bin=bin, fit=fit, peak=peak, peakPrint=peakPrint,
                    n_peaks=nrow(peak@peakList), peak_width=median(peak@peakList$peakSize),
                    opt_sig_model=opt_sig_model ) )
            } )
        }
        
        # summarize results
        
        peakSetFinal <- c()
        for ( i in 1:length(out) ) {
            peakSetFinal <- rbind( peakSetFinal, out[[i]]$peakPrint )
        }
        
        resultList <- list()
        resultList$chrID <- resultList$n_peaks <- 
            resultList$peak_width <- resultList$opt_sig_model <- rep( NA, length(out) )
        for ( i in 1:length(out) ) {
            resultList$chrID[i] <- out[[i]]$chrID
            resultList$n_peaks[i] <- out[[i]]$n_peaks
            resultList$peak_width[i] <- out[[i]]$peak_width
            resultList$opt_sig_model[i] <- out[[i]]$opt_sig_model
        }
    } else {
        ###############################################################
        #                                                             #
        #                   genome-wide analysis                      #
        #                                                             #
        ###############################################################
        
        # read in bin-level files                           
        
        cat( "Info: analyzing bin-level files...\n" )
        
        setwd( binfileDir )
        if ( PET == TRUE ) {
            chip_file <- list.files( path=binfileDir,
                paste(basename(chipFile),"_bin",binSize,".txt",sep="")  ) 
            input_file <- list.files( path=binfileDir,
                paste(basename(controlFile),"_bin",binSize,".txt",sep="")  ) 
        } else {
            chip_file <- list.files( path=binfileDir,
                paste(basename(chipFile),"_fragL",fragLen,"_bin",binSize,".txt",sep="")  ) 
            input_file <- list.files( path=binfileDir,
                paste(basename(controlFile),"_fragL",fragLen,"_bin",binSize,".txt",sep="")  ) 
        }
        
        # model fitting & peak calling
        # check whether rparallel is available. if so, use it.
        
        cat( "Info: fitting MOSAiCS model & call peaks...\n" )
        
        out <- list()
        
        # read in bin-level file
        
        out$bin <- readBins( type=c("chip","input"), fileName=c(chip_file,input_file),
            parallel=parallel, nCore=nCore )
  
        # fit model
        
        out$fit <- mosaicsFit( out$bin, analysisType=analysisType, bgEst=bgEst, d=d,
            parallel=parallel, nCore=nCore )    
        
        # call peaks
        
        if ( signalModel=="BIC" ) {
            # if not specified, use BIC
            
            if ( out$fit@bic1S < out$fit@bic2S ) {
                out$peak <- mosaicsPeak( out$fit, signalModel="1S", 
                    FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )      
                out$opt_sig_model <- "One-signal-component model"
            } else {        
                out$peak <- mosaicsPeak( out$fit, signalModel="2S", 
                    FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )      
                out$opt_sig_model <- "Two-signal-component model"
            }
        } else {                                                  
            out$peak <- mosaicsPeak( out$fit, signalModel=signalModel, 
                FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )
            if ( signalModel=="1S" ) {
                out$opt_sig_model <- "One-signal-component model"
            } else {
                out$opt_sig_model <- "Two-signal-component model"
            }
        }  
  
        # keep results
        
        peakSetFinal <- print(out$peak)
        
        peakPrint <- split( peakSetFinal, peakSetFinal$chrID )        
        out$chrID <- names(peakPrint)
        out$n_peaks <- unlist( lapply( peakPrint, nrow ) )         
        out$peak_width <- unlist( lapply( peakPrint, function(x) { median(x$peakSize) } ) )
        
        resultList <- list()
        resultList$chrID <- out$chrID
        resultList$n_peaks <- out$n_peaks
        resultList$peak_width <- out$peak_width
        resultList$opt_sig_model <- rep( out$opt_sig_model, length(resultList$chrID) )
    }
  
    # write peak calling results                      
  
    cat( "Info: writing the peak list...\n" )
  
    for ( ff in 1:length(peakFileFormat) ) {
        if ( peakFileFormat[ff] == "txt" ) {
            .exportTXT( peakList=peakSetFinal, filename=peakFile[ff] )
        } else if ( peakFileFormat[ff] == "bed" ) {
            .exportBED( peakList=peakSetFinal, filename=peakFile[ff] )
        } else if ( peakFileFormat[ff] == "gff" ) {
            .exportGFF( peakList=peakSetFinal, filename=peakFile[ff] )
        } else {
            stop( "Inappropriate peak file format!" )  
        }
    }
  
    # report: summary
  
    cat( "Info: generating reports...\n" )
  
    if ( reportSummary ) {    
        .reportSummary( summaryFile=summaryFile, resultList=resultList,
            chipFile=chipFile, chipFileFormat=chipFileFormat, 
            controlFile=controlFile, controlFileFormat=controlFileFormat, 
            binfileDir=binfileDir, 
            peakFile=peakFile, peakFileFormat=peakFileFormat, 
            byChr=byChr, FDR=FDR, fragLen=fragLen, binSize=binSize, capping=capping, 
            analysisType=analysisType, d=d, 
            signalModel=signalModel, maxgap=maxgap, minsize=minsize, thres=thres )
    }
        
    # GOF
  
    if ( reportGOF ) {    
        pdf(gofFile)
        
        if ( byChr ) {
            for ( i in 1:length(out) ) {
                chrID <- out[[i]]$chrID
                fit <- out[[i]]$fit   
                
                # chrID
                
                plot( 0, 0, type="n", axes=F, ann=F )
                text( 0, 0, chrID, cex=4 )
                
                # GOF
                
                plot(fit)
            }
        } else {
            fit <- out$fit
            plot(fit)
        }        
        
        dev.off()
    }
        
    # exploratory analysis
  
    if ( reportExploratory ) {
        pdf(exploratoryFile)
        
        if ( byChr ) {
            for ( i in 1:length(out) ) {
                chrID <- out[[i]]$chrID
                bin <- out[[i]]$bin
                
                # chrID
                
                plot( 0, 0, type="n", axes=F, ann=F )
                text( 0, 0, chrID, cex=4 )
                
                # exploratory plots 
                
                plot( bin )
                if ( analysisType=="IO" ) {
                    plot( bin, plotType="input" )
                }
                if ( analysisType=="OS" ) {
                    plot( bin, plotType="M" )       
                    plot( bin, plotType="GC" )            
                }
                if ( analysisType=="TS" ) {
                    plot( bin, plotType="M" )       
                    plot( bin, plotType="GC" )          
                    plot( bin, plotType="M|input" )         
                    plot( bin, plotType="GC|input" )        
                }
            }
        } else {
            bin <- out$bin
            
            # exploratory plots 
            
            plot( bin )
            if ( analysisType=="IO" ) {
                plot( bin, plotType="input" )
            }
            if ( analysisType=="OS" ) {
                plot( bin, plotType="M" )       
                plot( bin, plotType="GC" )            
            }
            if ( analysisType=="TS" ) {
                plot( bin, plotType="M" )       
                plot( bin, plotType="GC" )          
                plot( bin, plotType="M|input" )         
                plot( bin, plotType="GC|input" )        
            }
        }
        
        dev.off()
    }
}
