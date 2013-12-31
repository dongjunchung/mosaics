
# read alignment files and construct bin-level files

generateWig <- function( infile=NULL, fileFormat=NULL, outfileLoc="./", 
    byChr=FALSE, useChrfile=FALSE, chrfile=NULL, excludeChr=NULL, 
    PET=FALSE, fragLen=200, span=200, capping=0, normConst=1, perl = "perl" )
{   
    # preprocessing perl script embedded in "mosaics/inst/Perl/"
    
    if ( PET == TRUE ) {
        script <- "readfile2wig_PET.pl"
        allFormat <- c( "eland_result", "sam" )
        allFormatName <- c( "Eland result", "SAM" )
    } else {
        script <- "readfile2wig_SET.pl"
        allFormat <- c( "eland_result", "eland_extended", "eland_export", 
            "bowtie", "sam", "bed", "csem" )
        allFormatName <- c( "Eland result", "Eland extended", "Eland export", 
            "Bowtie default", "SAM", "BED", "CSEM" )
    }
    
    # Check whether perl exists
    
    CMD <- paste( perl, "-v" )
    res <- system( CMD, intern = TRUE, ignore.stderr = TRUE )
  
    if ( length(res) == 0 ) {
        # cannot proceed if perl does not exist
        
        stop( "Perl is not found on your system! Either check $PATH if installed or please install Perl." )
    } else {
        # process read files into bin-level files if perl exists
        
        # check whether minimal options are missing
        
        if ( length(infile) != 1 || is.null(infile) )
        {
            stop( "Please specify the name of the aligned read file!" )
        }       
        
        if ( length(fileFormat) != 1 || is.null(fileFormat) )
        {
            stop( "Please specify aligned read file format! Read '?generateWig' for supported file formats" )
        }   
        
        # check file format specification
        
        if ( length(which(!is.na(match( fileFormat, allFormat )))) == 0 )
        {
            stop( "Unsupported aligned read file format! Read '?generateWig' for supported file formats" )
        }
        
        # if useChrfile is TRUE & excludeChr is NOT null, then ignore excludeChr
        
        if ( useChrfile & !is.null(excludeChr) ) {
            message( "User set 'useChrfile' as TRUE and also provided 'excludeChr'." )
            message( "'excludeChr' argument will be ignored." )
            excludeChr <- NULL
        }
        
        # print out processing settings:
        # by default, set fragment length = 200, bin size = fragment length, capping = 0.
        
        fileFormatName <- allFormatName[ match( fileFormat, allFormat ) ]
        
        cat( "------------------------------------------------------------\n" )
        cat( "Info: setting summary\n" )
        cat( "------------------------------------------------------------\n" )
        cat( "Name of aligned read file:", infile, "\n" )
        cat( "Aligned read file format:", fileFormatName, "\n" )
        cat( "Directory of processed wig files:", outfileLoc, "\n" )
        cat( "span of the wig files:", span, "\n" )
        cat( "Normalizing constant:", normConst, "\n" )
        cat( "Construct wig files by chromosome?", ifelse(byChr,"Y","N"), "\n" ) 
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
        if ( capping > 0 ) {
            cat( "Maximum number of reads allowed in each nucleotide:", capping, "\n" )
        }
        cat( "------------------------------------------------------------\n" )
        
        
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
                " ", span, 
                " ", normConst,
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
                " ", span,
                " ", normConst,
                " ", fragLen,
                " ", capping, 
                " ", ifelse(byChr,"Y","N"),
                " ", ifelse(useChrfile,"Y","N"),
                " ", ifelse(!is.null(chrfile),chrfile,"-"),
                " ", paste(excludeChrVec,collpase=" "), sep="" )
        }
        
        res <- system( CMD, intern = TRUE )
        
        message( "Info: done!" )   
            
        
        # print out processing results
        
        infilename <- basename( infile )
            # extract only filename from infile
        
        if ( PET == TRUE ) {    
            if ( byChr ) {
                outfileName <- list.files( path=outfileLoc,
                    pattern=paste(infilename,"_span",span,"_.*.wig",sep="") )
            } else {
                outfileName <- paste(infilename,"_span",span,".wig",sep="")
            }
        } else {
            if ( byChr ) {
                outfileName <- list.files( path=outfileLoc,
                    pattern=paste(infilename,"_fragL",fragLen,"_span",span,"_.*.wig",sep="") )
            } else {
                outfileName <- paste(infilename,"_fragL",fragLen,"_span",span,".wig",sep="")
            }
        }
        
        cat( "------------------------------------------------------------\n" )
        cat( "Info: processing summary\n" )
        cat( "------------------------------------------------------------\n" )    
        cat( "Directory of processed wig files:", outfileLoc, "\n" )
        if ( byChr ) {
            cat( "List of processed wig files:\n" )
            for ( i in 1:length(outfileName) ) {                
                cat( "- ",outfileName[i],"\n", sep="" )
            }
        } else {
            cat( "Processed wig file: ",outfileName,"\n", sep="" )   
        }
        cat( "------------------------------------------------------------\n" )
    }
}
