\name{export}
\alias{export}
\alias{export,MosaicsPeak-method}
\title{
Export peak calling results to text files
}
\description{
Export peak calling results to text files in TXT, BED, or GFF file formats.
}
\usage{
export(object, ...)
\S4method{export}{MosaicsPeak}( object, type=NA, filename=NA )
}
\arguments{
    \item{object}{Object of class \code{MosaicsPeak},
        peak calling results obtained using method \code{mosaicsPeak}. }
    \item{type}{Format of the exported file.
        Possible values are \code{"txt"}, \code{"bed"}, \code{"gff"},
        \code{"narrowPeak"}, and \code{"broadPeak"}.
        See Details. }
    \item{filename}{Name of the exported file. }    
    \item{...}{ Other parameters to be passed through to generic \code{export}.}
}
\details{
TXT file format (\code{type="txt"}) exports peak calling results in the most informative way.
Columns include chromosome ID, peak start position, peak end position, peak width, 
-log10 transformed average posterior probability, -log10 transformed minimum posterior probability, 
average of -log10 transformed posterior probability,
average ChIP tag count, maximum ChIP tag count (always), 
average input tag count, average input tag count scaled by sequencing depth, 
average log base 2 ratio of ChIP over input tag counts 
(if matched control sample is also provided), 
average mappability score, and average GC content score
(when mappability and GC content scores are used in the analysis) in each peak. 
\code{type="bed"} and \code{type="gff"} export peak calling results in standard BED and GFF file formats,
respectively, where score is the average ChIP tag counts in each peak.
\code{type="narrowPeak"} and \code{type="broadPeak"} export peak calling results in ENCODE narrowPeak and broadPeak file formats,
respectively, where score, signalValue, pValue, and qValue are average log base 2 ratio of ChIP over input tag counts (or average ChIP tag count, if matched control sample is not provided), ChIP signal at the summit, -log10 transformation of minimum posterior probability (logMinP), and average of -log10 transformed posterior probability (aveLogP), respectively, in each peak.
If no peak is detected, \code{export} method will not generate any file.
}
\references{
Kuan, PF, D Chung, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.

Chung, D, Zhang Q, and Keles S (2014), "MOSAiCS-HMM: A model-based approach for detecting regions of histone modifications from ChIP-seq data", Datta S and Nettleton D (eds.), \emph{Statistical Analysis of Next Generation Sequencing Data}, Springer.
}
\author{ Dongjun Chung, Pei Fen Kuan, Rene Welch, Sunduz Keles }
\seealso{
\code{\link{mosaicsPeak}}, \code{\linkS4class{MosaicsPeak}}.
}
\examples{
\dontrun{
library(mosaicsExample)

constructBins( infile=system.file( file.path("extdata","wgEncodeBroadHistoneGm12878H3k4me3StdAlnRep1_chr22_sorted.bam"), package="mosaicsExample"), 
    fileFormat="bam", outfileLoc="~/", 
    byChr=FALSE, useChrfile=FALSE, chrfile=NULL, excludeChr=NULL, 
    PET=FALSE, fragLen=200, binSize=200, capping=0 )
constructBins( infile=system.file( file.path("extdata","wgEncodeBroadHistoneGm12878ControlStdAlnRep1_chr22_sorted.bam"), package="mosaicsExample"), 
    fileFormat="bam", outfileLoc="~/", 
    byChr=FALSE, useChrfile=FALSE, chrfile=NULL, excludeChr=NULL, 
    PET=FALSE, fragLen=200, binSize=200, capping=0 )

binHM <- readBins( type=c("chip","input"),
    fileName=c( "~/wgEncodeBroadHistoneGm12878H3k4me3StdAlnRep1_chr22_sorted.bam_fragL200_bin200.txt",
    "~/wgEncodeBroadHistoneGm12878ControlStdAlnRep1_chr22_sorted.bam_fragL200_bin200.txt" ) )
fitHM <- mosaicsFit( binHM, analysisType="IO", bgEst="rMOM" )
hmmHM <- mosaicsFitHMM( fitHM, signalModel = "2S", 
  init="mosaics", init.FDR = 0.05, parallel=TRUE, nCore=8 )
peakHM <- mosaicsPeakHMM( hmmHM, FDR = 0.05, decoding="posterior",
  thres=10, parallel=TRUE, nCore=8 )

export( peakHM, type = "txt", filename = "./peakHM.txt" )
export( peakHM, type = "bed", filename = "./peakHM.bed" )
export( peakHM, type = "gff", filename = "./peakHM.gff" )

# read-level data is needed to loaded and peak summits need to be identified
# to export in narrowPeak and broadPeak file formats

peakHM <- extractReads( peakHM,
  chipFile=system.file( file.path("extdata","wgEncodeBroadHistoneGm12878H3k4me3StdAlnRep1_chr22_sorted.bam"), package="mosaicsExample"),
  chipFileFormat="bam", chipPET=FALSE, chipFragLen=200,
  controlFile=system.file( file.path("extdata","wgEncodeBroadHistoneGm12878ControlStdAlnRep1_chr22_sorted.bam"), package="mosaicsExample"), 
  controlFileFormat="bam", controlPET=FALSE, controlFragLen=200, parallel=TRUE, nCore=8 )
peakHM

peakHM <- findSummit( peakHM, parallel=TRUE, nCore=8 )
peakHM <- adjustBoundary( peakHM, parallel=TRUE, nCore=8 )
peakHM <- filterPeak( peakHM, parallel=TRUE, nCore=8 )
export( peakHM, type = "narrowPeak", filename = "./peakHM.narrowPeak" )
export( peakHM, type = "broadPeak", filename = "./peakHM.broadPeak" )
}
}
\keyword{models}
\keyword{methods}
