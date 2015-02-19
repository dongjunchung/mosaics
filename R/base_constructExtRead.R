
# read alignment files and construct bin-level files

.constructExtRead <- function( infile=NULL, outfile=NULL, summaryfile=NULL,
    fileFormat="eland_result", PET=FALSE, fragLen=200, perl = "perl" )
{   
    # preprocessing perl script embedded in "dpeak/inst/Perl/"
    
    if ( PET ) {
        script <- "extract_coord_PET.pl"
        #allFormat <- "eland_result"
        #allFormatName <- "Eland result"
    } else {
        script <- "extract_coord_SET.pl"
        #allFormat <- c( "eland_result", "eland_extended", "eland_export", "bowtie", "sam", "bed" )
        #allFormatName <- c( "Eland result", "Eland extended", "Eland export", "Bowtie default", "SAM", "BED" )
    }
    
    # Check whether perl exists
    
    CMD <- paste( perl, "-v" )
    res <- system( CMD, intern = TRUE, ignore.stderr = TRUE )
  
    if ( length(res) == 0 ) {
        # cannot proceed if perl does not exist
        
        stop( "Perl is not found on your system! Either check $PATH if installed or please install Perl." )
    } else {
        # process read files into bin-level files if perl exists
            
        
        # get path to the perl code (unified script for all file formats)
        
        Fn.Path <- system.file( file.path("Perl",script), package="mosaics" )
                
        
        # process read file to bin-level files using perl codes
        
        if ( PET ) {
            CMD <- paste( perl, 
                " ", Fn.Path,         
                " ", infile, 
                " ", outfile, 
                " ", summaryfile, 
                " ", fileFormat, sep="" )
        } else {
            CMD <- paste( perl, 
                " ", Fn.Path,         
                " ", infile, 
                " ", outfile, 
                " ", summaryfile, 
                " ", fileFormat, 
                " ", fragLen, sep="" )
        }
        
        res <- system( CMD, intern = TRUE )
    }
}
