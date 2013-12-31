\name{constructBins}
\alias{constructBins}
\title{
Construct bin-level ChIP-sep data from an aligned read file
}
\description{
Preprocess and construct bin-level ChIP-sep data from an aligned read file.
}
\usage{
constructBins( infile=NULL, fileFormat=NULL, outfileLoc="./", 
    byChr=FALSE, useChrfile=FALSE, chrfile=NULL, excludeChr=NULL, 
    PET=FALSE, fragLen=200, binSize=200, capping=0, perl = "perl" )
}
\arguments{
  \item{infile}{
    Name of the aligned read file to be processed.
}
  \item{fileFormat}{
    Format of the aligned read file to be processed.
    Currently, \code{constructBins} permits the following aligned read file formats
    for SET data (\code{PET = FALSE}):
        \code{"eland_result"} (Eland result), \code{"eland_extended"} (Eland extended),
        \code{"eland_export"} (Eland export), \code{"bowtie"} (default Bowtie),
        \code{"sam"} (SAM), \code{"bed"} (BED), and \code{"csem"} (CSEM).
    For PET data (\code{PET = TRUE}), the following aligned read file formats are allowed:
        \code{"eland_result"} (Eland result) and \code{"sam"} (SAM).
}
  \item{outfileLoc}{
    Directory of processed bin-level files. 
    By default, processed bin-level files are exported to the current directory. 
}
  \item{byChr}{
    Construct separate bin-level file for each chromosome? 
    Possible values are \code{TRUE} or \code{FALSE}.
    If \code{byChr=FALSE}, bin-level data for all chromosomes are exported to one file. 
    If \code{byChr=TRUE}, bin-level data for each chromosome is exported to a separate file.
    Default is \code{FALSE}.
}
  \item{useChrfile}{
    Is the file for chromosome info provided?
    Possible values are \code{TRUE} or \code{FALSE}.
    If \code{useChrfile=FALSE}, it is assumed that the file for chromosome info is not provided. 
    If \code{useChrfile=TRUE}, it is assumed that the file for chromosome info is provided.
    Default is \code{FALSE}.
}
  \item{chrfile}{
    Name of the file for chromosome info.
    In this file, the first and second columns are ID and size of each chromosome, respectively.
}
  \item{excludeChr}{
    Vector of chromosomes that will be excluded from the analysis.
    This argument is ignored if \code{useChrfile=TRUE}.
}  
  \item{PET}{
    Is the file paired-end tag (PET) data?
    If \code{PET=FALSE}, it is assumed that the file is SET data. 
    If \code{PET=TRUE}, it is assumed that the file is PET data.    
    Default is \code{FALSE} (SET data).
}
  \item{fragLen}{
    Average fragment length. Default is 200.
    This argument is ignored if \code{PET=TRUE}.
}
  \item{binSize}{
    Size of bins. Default is 200.
}
  \item{capping}{
    Maximum number of reads allowed to start at each nucleotide position. 
    To avoid potential PCR amplification artifacts, the maximum number of reads
    that can start at a nucleotide position is capped at \code{capping}. 
    Capping is not applied if non-positive value is used for \code{capping}.
    Default is 0 (no capping).
}
  \item{perl}{
    Name of the perl executable to be called. Default is \code{"perl"}.
}
}
\details{
Bin-level files are constructed from the aligned read file and 
exported to the directory specified in \code{outfileLoc} argument.
If \code{byChr=FALSE}, bin-level files are named 
as \code{[infileName]_fragL[fragLen]_bin[binSize].txt} for SET data (\code{PET = FALSE})
and \code{[infileName]_bin[binSize].txt} for PET data (\code{PET = TRUE}).
If \code{byChr=TRUE}, bin-level files are named 
as \code{[infileName]_fragL[fragLen]_bin[binSize]_[chrID].txt}for SET data (\code{PET = FALSE})
and \code{[infileName]_bin[binSize]_[chrID].txt} for PET data (\code{PET = TRUE}),
where \code{chrID} is chromosome IDs that reads align to.
These chromosome IDs are extracted from the aligned read file.

If the file for chromosome information is provided (\code{useChrfile=TRUE} and \code{chrfile} is not NULL),
only the chromosomes specified in the file will be considered.
Chromosomes that are specified in \code{excludeChr}
will not be included in the processed bin-level files.
\code{excludeChr} argument is ignored if \code{useChrfile=TRUE}.
Constructed bin-level files can be loaded into the R environment using the method \code{readBins}.

\code{constructBins} currently supports the following aligned read file formats
for SET data (\code{PET = FALSE}): 
Eland result (\code{"eland_result"}), Eland extended (\code{"eland_extended"}),
Eland export (\code{"eland_export"}), default Bowtie (\code{"bowtie"}), 
SAM (\code{"sam"}), BED (\code{"bed"}), and CSEM (\code{"csem"}).
For PET data (\code{PET = TRUE}), the following aligned read file formats are allowed:
\code{"eland_result"} (Eland result) and \code{"sam"} (SAM).

If input file format is neither BED nor CSEM BED,
this method retains only reads mapping uniquely to the reference genome. 
}
\value{
Processed bin-level files are exported to the directory specified in \code{outfileLoc}.
}
\references{
Kuan, PF, D Chung, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.
}
\author{ Dongjun Chung, Pei Fen Kuan, Sunduz Keles }
\seealso{
\code{\link{readBins}}, \code{\linkS4class{BinData}}.
}
\examples{
\dontrun{
constructBins( infile="/scratch/eland/STAT1_eland_results.txt", 
    fileFormat="eland_result", outfileLoc="/scratch/eland/", 
    byChr=FALSE, useChrfile=FALSE, chrfile=NULL, excludeChr="chrM", 
    PET = FALSE, fragLen=200, binSize=200, capping=0 )
}
}
\keyword{models}
\keyword{methods}
