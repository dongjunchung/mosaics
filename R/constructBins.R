
# read alignment files and construct bin-level files

constructBins <- function( infile=NULL, fileFormat=NULL, outfileLoc="./", 
    byChr=FALSE, useChrfile=FALSE, chrfile=NULL, excludeChr=NULL, 
    PET=FALSE, fragLen=200, binSize=200, capping=0, perl = "perl" )
{   
    # preprocessing perl script embedded in "mosaics/inst/Perl/"
    
    if ( PET == TRUE ) {
        script <- "process_readfiles_PET.pl"
        allFormat <- c( "eland_result", "sam", "bam" )
        allFormatName <- c( "Eland result", "SAM", "BAM" )
    } else {
        script <- "process_readfiles_SET.pl"
        allFormat <- c( "eland_result", "eland_extended", "eland_export", 
            "bowtie", "sam", "bam", "bed", "csem" )
        allFormatName <- c( "Eland result", "Eland extended", "Eland export", 
            "Bowtie default", "SAM", "BAM", "BED", "CSEM" )
    }
  
    
    # check whether minimal options are missing
    
    if ( length(infile) != 1 || is.null(infile) )
    {
        stop( "Please specify the name of the aligned read file!" )
    }       
    
    if ( length(fileFormat) != 1 || is.null(fileFormat) )
    {
        stop( "Please specify aligned read file format! Read '?constructBins' for supported file formats" )
    }   
    
    # check file format specification
    
    if ( length(which(!is.na(match( fileFormat, allFormat )))) == 0 )
    {
        stop( "Unsupported aligned read file format! Read '?constructBins' for supported file formats" )
    }
    
    # if useChrfile is TRUE & excludeChr is NOT null, then ignore excludeChr
    
    if ( useChrfile & !is.null(excludeChr) ) {
        message( "User set 'useChrfile' as TRUE and also provided 'excludeChr'." )
        message( "'excludeChr' argument will be ignored." )
        excludeChr <- NULL
    }
    
    # capping is currently not supported for BAM file format
    
    if ( fileFormat == "BAM" & capping > 0 ) {
        message( "Capping is currently not allowed for the BAM file format." )        
    }
  
    
    # print out processing settings:
    # by default, set fragment length = 200, bin size = fragment length, capping = 0.
    
    fileFormatName <- allFormatName[ match( fileFormat, allFormat ) ]
    
    cat( "------------------------------------------------------------\n" )
    cat( "Info: setting summary\n" )
    cat( "------------------------------------------------------------\n" )
    cat( "Name of aligned read file:", infile, "\n" )
    cat( "Aligned read file format:", fileFormatName, "\n" )
    cat( "Directory of processed bin-level files:", outfileLoc, "\n" )
    cat( "Construct bin-level files by chromosome?", ifelse(byChr,"Y","N"), "\n" ) 
    cat( "Is file for chromosome info provided?", ifelse(useChrfile,"Y","N"), "\n" )   
    if ( useChrfile == TRUE ) { 
        cat( "Name of file for chromosome info: ", chrfile, "\n" )
    }
    if ( !is.null(excludeChr) ) {
        cat( "List of chromosomes to be excluded:", paste(excludeChr,collapse=", "), "\n" )
    }
    if ( PET == FALSE ) {
        cat( "Data type: Single-end tag (SET)\n" )
        cat( "Average fragment length:", fragLen, "\n" )
    } else {
        cat( "Data type: Paired-end tag (PET)\n" )
    }
    cat( "Bin size:", binSize, "\n" )
    if ( fileFormat != "BAM" & capping > 0 ) {
        cat( "Maximum number of reads allowed in each nucleotide:", capping, "\n" )
    }
    cat( "------------------------------------------------------------\n" )

  
    if ( fileFormat == "bam" ) {
      
      # check whether BAM index exists. Otherwise, generate BAM index
      
      bamName <- list.files( path=dirname(infile), pattern=basename(infile) )
      if ( length(grep( "bai", bamName )) > 0 ) {
        message( "Use the provided BAM index file." )
      } else {
        message( "BAM index file does not exist. Generating BAM index file..." )
        indexBam( infile )
      }
      
      # chromosome information
      
      if ( useChrfile ) {
        # if chrfile is provided, use it
        
        message( "Use the chromosome information provided in 'chrfile'." )
        
        chrdata <- read.table( chrfile, header=FALSE, stringsAsFactors=FALSE )
        chrnames <- chrdata[,1]
        chrlen <- chrdata[,2]
        names(chrlen) <- chrnames
        
      } else {      
        
        # check BAM file for chromosome information
        
        message( "Chromosome information is extracted from the BAM file." )
        
        bf <- BamFile( infile )
        baminfo <- seqinfo( bf )
        chrnames <- seqnames( baminfo )
        chrlen <- seqlengths( baminfo )
        
      }
      
      # calculate sequencing depth
      
      seqDepth <- countBam(infile)$records
      
      # file name for genome-wide file
        
      infilename <- basename( infile )
      
      if ( byChr == FALSE ) {
        
        if ( PET == TRUE ) {
          outfile <- file.path( outfileLoc, paste( infilename,"_bin",binSize,".txt", sep="" ) )
        } else {
          outfile <- file.path( outfileLoc, paste( infilename,"_fragL",fragLen,"_bin",binSize,".txt", sep="" ) )
        }
        
        # record sequencing depth
        
        #cat( file=outfile )
        cat( "# sequencing depth: ",seqDepth,"\n", sep="", file=outfile )
      }
      
      # process chromosome by chromosome
        
      message( "Info: reading the aligned read file and processing it into bin-level files..." )
      
      for ( chr in chrnames ) {
        
        # skip if in the list of excludeChr or not in the chrfile
        
        if ( !useChrfile && !is.null(excludeChr) && !is.na(match( chr, excludeChr )[1]) ) { next }
        
        # generate bins
        
        bindata <- tileGenome( seqlengths=chrlen[chr], tilewidth=binSize,
          cut.last.tile.in.chrom=TRUE )
        suppressWarnings( start(bindata) <- start(bindata) - 1 )
        suppressWarnings( end(bindata) <- end(bindata) - 1 )
          # to match with perl scripts
        
        # load BAM file
        
        param <- ScanBamParam( which=GRanges( seqnames = chr, IRanges( 1, chrlen[[chr]] ) ) )
    	  
    	  if ( PET == FALSE ) {
    		  suppressWarnings( greads <- readGAlignments( infile, param = param, use.names = FALSE ) )
    		  suppressWarnings( greads <- as( greads, "GRanges" ) )
    		  suppressWarnings( greads <- resize( greads, fragLen ) )
    	  } else {
    		  #suppressWarnings( greads <- readAlignmentsPairsFromBam( infile, param = param ) )
    		  #suppressWarnings( greads <- GRanges( seqnames = seqnames(greads),
    		  #	ranges = IRanges( start=start(left(greads)), end=end(right(greads)) ),
    			#  strand = Rle( "*", length(greads) ) )
    			#)
      
          suppressWarnings( greads <- readGAlignmentPairsFromBam( readfile, param = param ) )
    
          snms = seqnames(greads)
          starts = start(left(greads))
          ends = end(right(greads))
          
          # remove reads with negative widths         
          idx = (starts >= ends)
          if(any(idx)){
            warning("Removing ",sum(idx)," reads, due to negative read lengths")
            snms = snms[!idx]
            starts = starts[!idx]
            ends = ends[!idx]
          }
          suppressWarnings( greads <- GRanges( seqnames = snms,ranges = IRanges( start=starts, end=ends),strand = "*"))
          rm(snms,starts,ends)
    	  }
        
        # summarize counts
        
        counts <- countOverlaps( bindata, greads )
        #bindata$counts <- counts
        
        # file names for chromosome-wise files
        
        if ( byChr ) {
          
          if ( PET == TRUE ) {
            outfile <- file.path( outfileLoc, paste( infilename,"_bin",binSize,"_",chr,".txt", sep="" ) )
          } else {
            outfile <- file.path( outfileLoc, paste( infilename,"_fragL",fragLen,"_bin",binSize,"_",chr,".txt", sep="" ) )
          }
        
          # record sequencing depth
          
          #cat( file=outfile )
          cat( "# sequencing depth: ",seqDepth,"\n", sep="", file=outfile )
        }
      
        # no scientific notation
      
        scipenOrg <- options()$scipen
        options( scipen=999 )
        
        # write bin-level data file
     
        #for ( i in 1:length(bindata) )
        #{
        #    cat( chr, start(bindata)[i], counts[i], file=outfile, sep="\t", append=TRUE )
        #    cat( "\n", file=outfile, append=TRUE )
        #}
        write.table( data.frame( chr, start(bindata), counts ),
          file=outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE,
          append=TRUE )
    
        # recover original setting for scientific notation
      
        options( scipen=scipenOrg )
        
        rm( bindata, greads, counts )
        gc()
      }
       
    } else {
       
      # Check whether perl exists
      
      CMD <- paste( perl, "-v" )
      res <- system( CMD, intern = TRUE, ignore.stderr = TRUE )
    
      if ( length(res) == 0 ) {
          # cannot proceed if perl does not exist
          
          stop( "Perl is not found on your system! Either check $PATH if installed or please install Perl." )
      } else {
          # process read files into bin-level files if perl exists
        
        
          # get path to the perl code (unified script for all file formats)
          
          Fn.Path <- system.file( file.path("Perl",script), package="mosaics")
  
          
          # process read file to bin-level files using perl codes
          
          message( "Info: reading the aligned read file and processing it into bin-level files..." )
          
          if ( capping <= 0 ) {
              capping <- 0
          }
          if ( is.null(excludeChr) ) {
              excludeChrVec <- ""
          } else {
              excludeChrVec <- paste( excludeChr, collapse=" " )
          }
          
          if ( PET == TRUE ) {
              CMD <- paste( perl, 
                  " ", Fn.Path,
                  " ", infile, 
                  " ", outfileLoc, 
                  " ", fileFormat,
                  " ", binSize, 
                  " ", capping,
                  " ", ifelse(byChr,"Y","N"),
                  " ", ifelse(useChrfile,"Y","N"),
                  " ", ifelse(!is.null(chrfile),chrfile,"-"),
                  " ", paste(excludeChrVec,collpase=" "), sep="" )
          } else {    
              CMD <- paste( perl, 
                  " ", Fn.Path,
                  " ", infile, 
                  " ", outfileLoc, 
                  " ", fileFormat, 
                  " ", fragLen, 
                  " ", binSize, 
                  " ", capping, 
                  " ", ifelse(byChr,"Y","N"),
                  " ", ifelse(useChrfile,"Y","N"),
                  " ", ifelse(!is.null(chrfile),chrfile,"-"),
                  " ", paste(excludeChrVec,collpase=" "), sep="" )
          }
          
          res <- system( CMD, intern = TRUE )
      }
    }
  
    message( "Info: done!" ) 
  
  
    # print out processing results
          
    infilename <- basename( infile )
        # extract only filename from infile
    
    if ( PET == TRUE ) {    
        if ( byChr ) {
            outfileName <- list.files( path=outfileLoc,
                pattern=paste(infilename,"_bin",binSize,"_.*.txt",sep="") )
        } else {
            outfileName <- paste(infilename,"_bin",binSize,".txt",sep="")
        }
    } else {
        if ( byChr ) {
            outfileName <- list.files( path=outfileLoc,
                pattern=paste(infilename,"_fragL",fragLen,"_bin",binSize,"_.*.txt",sep="") )
        } else {
            outfileName <- paste(infilename,"_fragL",fragLen,"_bin",binSize,".txt",sep="")
        }
    }
  
    seqDepth <- read.table( file=file.path(outfileLoc,outfileName[1]), nrows=1, comment.char="" )[1,4]
    
    cat( "------------------------------------------------------------\n" )
    cat( "Info: processing summary\n" )
    cat( "------------------------------------------------------------\n" )    
    cat( "Directory of processed bin-level files:", outfileLoc, "\n" )
    if ( byChr ) {
        cat( "List of processed bin-level files:\n" )
        for ( i in 1:length(outfileName) ) {                
            cat( "- ",outfileName[i],"\n", sep="" )
        }
    } else {
        cat( "Processed bin-level file: ",outfileName,"\n", sep="" )   
    }
    cat( "Sequencing depth:",seqDepth, "\n" )
    cat( "------------------------------------------------------------\n" )
}
