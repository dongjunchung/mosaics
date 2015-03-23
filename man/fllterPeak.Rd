\name{filterPeak}
\alias{filterPeak}
\alias{filterPeak,MosaicsPeak-method}
\title{
Filter out potentially false positive peaks 
}
\description{
Filter out potentially false positive peaks from \code{MosaicsPeak} class object, which is a peak calling result.
}
\usage{
filterPeak( object, ... )
\S4method{filterPeak}{MosaicsPeak}( object, 
	minSummitQuantile=0.01, minRead=10, FC=50, minSignal = 20, minLen = 200, 
	normC=NA, parallel=FALSE, nCore=8 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{Object of class \code{MosaicsPeak},
        a peak list object obtained using either functions \code{mosaicsPeak} or \code{mosaicsPeakHMM}. }
    \item{minSummitQuantile}{Parameter for the filtering step #1. Peaks are filtered out if signal strengths at summits < \code{minSummitQuantile} quantile of signal strengths at summits across all the peak regions.}  
    \item{minRead}{Parameter for the filtering step #1. Peaks are filtered out if the number of reads in the peak region < \code{minRead}.}  
    \item{FC}{Parameter for the filtering step #1. Peaks are filtered out if the improvement of ChIP over matched control < \code{FC}. }  
    \item{minSignal}{Parameter for the filtering step #2. Peaks are filtered out if signal strengths at summits <= \code{minSignal}.}  
    \item{minLen}{Parameter for the filtering step #2. Peaks are filtered out if their lengths <= \code{minLen}.} 
    \item{normC}{Normalizing constant. If not provided,
        \code{normC} is estimated as ratio of sequencing depth of ChIP over matched control samples.} 
    \item{parallel}{Utilize multiple CPUs for parallel computing 
        using \code{"parallel"} package?
        Possible values are \code{TRUE} (utilize multiple CPUs) 
        or \code{FALSE} (do not utilize multiple CPUs).
        Default is \code{FALSE} (do not utilize multiple CPUs). }
    \item{nCore}{Number of CPUs when parallel computing is utilized. }
    \item{...}{ Other parameters to be passed through to generic \code{filterPeak}.}
}
\details{
\code{filterPeak} filters out potentially false positive peaks, based on signal strengths and peak lengths. While \code{filterPeak} can be applied to a peak list object obtained using either functions \code{mosaicsPeak} or \code{mosaicsPeakHMM}, \code{filterPeak} is developed and tested mainly for peak lists from MOSAiCS-HMM model (i.e., from function \code{mosaicsPeakHMM}). Note that \code{extractReads} should be run first because \code{filterPeak} is used.
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
\author{ Dongjun Chung, Pei Fen Kuan, Rene Welch, Sunduz Keles}
\seealso{
\code{\link{mosaicsPeak}}, \code{\link{mosaicsPeakHMM}}, 
\code{\link{extractReads}}, \code{\link{findSummit}}, \code{\link{adjustBoundary}}, 
\code{\linkS4class{MosaicsPeak}}.
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

peakHM <- extractReads( peakHM,
  chipFile=system.file( file.path("extdata","wgEncodeBroadHistoneGm12878H3k4me3StdAlnRep1_chr22_sorted.bam"), package="mosaicsExample"),
  chipFileFormat="bam", chipPET=FALSE, chipFragLen=200,
  controlFile=system.file( file.path("extdata","wgEncodeBroadHistoneGm12878ControlStdAlnRep1_chr22_sorted.bam"), package="mosaicsExample"), 
  controlFileFormat="bam", controlPET=FALSE, controlFragLen=200, parallel=TRUE, nCore=8 )
peakHM <- findSummit( peakHM, parallel=TRUE, nCore=8 )
peakHM <- adjustBoundary( peakHM, parallel=TRUE, nCore=8 )
peakHM <- filterPeak( peakHM, parallel=TRUE, nCore=8 )
}
}
\keyword{models}
\keyword{methods}
