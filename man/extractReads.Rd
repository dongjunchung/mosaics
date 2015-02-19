\name{extractReads}
\alias{extractReads}
\alias{extractReads,MosaicsPeak-method}
\title{
Load read-level data and extract reads corresponding to each peak region
}
\description{
Load read-level data and extract reads corresponding to each peak region in the \code{MosaicsPeak} class object, which is a peak calling result.
}
\usage{
extractReads( object, ... )
\S4method{extractReads}{MosaicsPeak}( object, chipFile=NULL, chipFileFormat=NULL,
    chipPET=FALSE, chipFragLen=200,
    controlFile=NULL, controlFileFormat=NULL, 
    controlPET=FALSE, controlFragLen=200,
    parallel=FALSE, nCore=8, tempDir=NULL, perl="perl" )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{Object of class \code{MosaicsPeak},
        a peak list object obtained using either functions \code{mosaicsPeak} or \code{mosaicsPeakHMM}. }
  \item{chipFile}{
    Name of the aligned read file for ChIP sample.
}
  \item{chipFileFormat}{
    Format of the aligned read file for ChIP sample.
    Currently, \code{constructBins} permits the following aligned read file formats
    for SET data (\code{chipPET = FALSE}):
        \code{"eland_result"} (Eland result), \code{"eland_extended"} (Eland extended),
        \code{"eland_export"} (Eland export), \code{"bowtie"} (default Bowtie),
        \code{"sam"} (SAM), \code{"bam"} (BAM), \code{"bed"} (BED), and \code{"csem"} (CSEM).
    For PET data (\code{chipPET = TRUE}), the following aligned read file formats are allowed:
        \code{"eland_result"} (Eland result), \code{"sam"} (SAM), and \code{"bam"} (BAM).
}
  \item{chipPET}{
    Is the file of ChIP sample paired-end tag (PET) data?
    If \code{chipPET=FALSE}, it is assumed that the file is SET data. 
    If \code{chipPET=TRUE}, it is assumed that the file is PET data.    
    Default is \code{FALSE} (SET data).
}
  \item{chipFragLen}{
    Average fragment length of ChIP sample. Default is 200.
    This argument is ignored if \code{chipPET=TRUE}.
}
  \item{controlFile}{
    Name of the aligned read file for matched control sample.
}
  \item{controlFileFormat}{
    Format of the aligned read file for matched control sample.
    Currently, \code{constructBins} permits the following aligned read file formats
    for SET data (\code{controlPET = FALSE}):
        \code{"eland_result"} (Eland result), \code{"eland_extended"} (Eland extended),
        \code{"eland_export"} (Eland export), \code{"bowtie"} (default Bowtie),
        \code{"sam"} (SAM), \code{"bam"} (BAM), \code{"bed"} (BED), and \code{"csem"} (CSEM).
    For PET data (\code{controlPET = TRUE}), the following aligned read file formats are allowed:
        \code{"eland_result"} (Eland result), \code{"sam"} (SAM), and \code{"bam"} (BAM).
}
  \item{controlPET}{
    Is the file of matched control sample paired-end tag (PET) data?
    If \code{controlPET=FALSE}, it is assumed that the file is SET data. 
    If \code{controlPET=TRUE}, it is assumed that the file is PET data.    
    Default is \code{FALSE} (SET data).
}
  \item{controlFragLen}{
    Average fragment length of matched control sample. Default is 200.
    This argument is ignored if \code{controlPET=TRUE}.
}
    \item{parallel}{Utilize multiple CPUs for parallel computing 
        using \code{"parallel"} package?
        Possible values are \code{TRUE} (utilize multiple CPUs) 
        or \code{FALSE} (do not utilize multiple CPUs).
        Default is \code{FALSE} (do not utilize multiple CPUs). }
    \item{nCore}{Number of CPUs when parallel computing is utilized. }
  \item{tempDir}{
    Directory to store temporary files. If \code{tempDir=NULL}, 
    \code{extractReads()} will use the temporary directory used by R.
}
  \item{perl}{
    Name of the perl executable to be called. Default is \code{"perl"}.
}
    \item{...}{ Other parameters to be passed through to generic \code{extractReads}.}
}
\details{
Read-level data is first loaded from aligned read files, and then the reads corresponding to each peak region are extracted and incorporated into the \code{MosaicsPeak} class object.
\code{extractReads} currently supports the following aligned read file formats
for SET data (\code{chipPET = FALSE} and \code{controlPET = FALSE}): 
Eland result (\code{"eland_result"}), Eland extended (\code{"eland_extended"}),
Eland export (\code{"eland_export"}), default Bowtie (\code{"bowtie"}), 
SAM (\code{"sam"}), BAM (\code{"bam"}), BED (\code{"bed"}), and CSEM (\code{"csem"}).
For PET data (\code{chipPET = FALSE} and \code{controlPET = FALSE}), 
the following aligned read file formats are allowed:
\code{"eland_result"} (Eland result), \code{"sam"} (SAM), and \code{"bam"} (BAM).

If input file format is neither BED nor CSEM BED,
this method retains only reads mapping uniquely to the reference genome. 
}
\value{
Construct \code{MosaicsPeak} class object.
}
\references{
Kuan, PF, D Chung, G Pan, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.

Chung, D, Zhang Q, and Keles S (2014), "MOSAiCS-HMM: A model-based approach for detecting regions of histone modifications from ChIP-seq data", Datta S and Nettleton D (eds.), \emph{Statistical Analysis of Next Generation Sequencing Data}, Springer.
}
\author{ Dongjun Chung, Pei Fen Kuan, Sunduz Keles}
\seealso{
\code{\link{mosaicsFit}}, \code{\link{mosaicsFitHMM}}, 
\code{\linkS4class{MosaicsFit}}, \code{\linkS4class{MosaicsHMM}}, \code{\linkS4class{MosaicsPeak}}.
}
\examples{
\dontrun{
library(mosaicsExample)
data(exampleFit)

exampleFitHMM <- mosaicsFitHMM( exampleFit, signalModel = "2S", init.FDR = 0.05 )
examplePeakHMM1 <- mosaicsPeakHMM( exampleFitHMM, decoding="viterbi" )
examplePeakHMM2 <- mosaicsPeakHMM( exampleFitHMM, FDR = 0.05, decoding="posterior" )
}
}
\keyword{models}
\keyword{methods}
