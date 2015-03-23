\name{mosaicsPeakHMM}
\alias{mosaicsPeakHMM}
\alias{mosaicsPeakHMM,MosaicsHMM-method}
\title{
Call broad peaks using fitted MOSAiCS-HMM model 
}
\description{
Call broad peaks using \code{MosaicsHMM} class object, which is a fitted MOSAiCS-HMM model.
}
\usage{
mosaicsPeakHMM( object, ... )
\S4method{mosaicsPeakHMM}{MosaicsHMM}( object, FDR=0.05, decoding="posterior",
    binsize=NA, maxgap=0, minsize=0, thres=0, 
    parallel=FALSE, nCore=8 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{Object of class \code{MosaicsHMM},
        a fitted MOSAiCS model obtained using function \code{mosaicsFitHMM}. }
    \item{FDR}{False discovery rate. Default is 0.05. Not relevant when \code{decoding="viterbi"}. }  
    \item{decoding}{Approach to determine the undelying state.
        Possible values are "viterbi" (Viterbi algorithm) and
        "posterior" (posterior decoding). Default is "posterior". }  
    \item{binsize}{Size of each bin. Value should be positive integer.
        If \code{binsize=NA}, \code{mosaicsPeakHMM} function calcuates the value from data.
        Default is \code{NA}. }        
    \item{maxgap}{Initial nearby peaks are merged if the distance (in bp) between them is less than \code{maxgap}. Default is 0. }     
    \item{minsize}{An initial peak is removed if its width is narrower than \code{minsize}. Default is 0. }
    \item{thres}{A bin within initial peak is removed if its ChIP tag counts are less than \code{thres}. Default is 0. }  
    \item{parallel}{Utilize multiple CPUs for parallel computing 
        using \code{"parallel"} package?
        Possible values are \code{TRUE} (utilize multiple CPUs) 
        or \code{FALSE} (do not utilize multiple CPUs).
        Default is \code{FALSE} (do not utilize multiple CPUs). }
    \item{nCore}{Number of CPUs when parallel computing is utilized. }
    \item{...}{ Other parameters to be passed through to generic \code{mosaicsHMM}.}
}
\details{
\code{mosaicsFitHMM} and \code{mosaicsPeakHMM} are developed to identify broad peaks such as histone modifications,
using Hidden Markov Model (HMM) approach, as proposed in Chung et al. (2014).
If you are interested in identifying narrow peaks such as transcription factor binding sites,
please use \code{mosaicsPeak} instead of \code{mosaicsFitHMM} and \code{mosaicsPeakHMM}.

\code{maxgap}, \code{minsize}, and \code{thres} are for refining initial peaks called
using specified \code{decoding} (and \code{FDR} if \code{decoding="posterior"}).
If you use a bin size shorter than the average fragment length of the experiment, 
we recommend to set \code{maxgap} to the average fragment length 
and \code{minsize} to the bin size.
If you set the bin size to the average fragment length or if bin size is larger than the average fragment length,
set \code{maxgap} to the average fragment length and
\code{minsize} to a value smaller than the average fragment length. See the vignette for further details.

Parallel computing can be utilized for faster computing
if \code{parallel=TRUE} and \code{parallel} package is loaded.
\code{nCore} determines number of CPUs used for parallel computing.
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
\code{\link{mosaicsFitHMM}}, 
\code{\link{extractReads}}, \code{\link{findSummit}}, \code{\link{adjustBoundary}}, \code{\link{filterPeak}}, 
\code{\linkS4class{MosaicsHMM}}, \code{\linkS4class{MosaicsPeak}}.
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
binHM
plot(binHM)
plot( binHM, plotType="input" )

fitHM <- mosaicsFit( binHM, analysisType="IO", bgEst="rMOM" )
fitHM
plot(fitHM)

hmmHM <- mosaicsFitHMM( fitHM, signalModel = "2S", 
  init="mosaics", init.FDR = 0.05, parallel=TRUE, nCore=8 )
hmmHM
plot(hmmHM)

peakHM <- mosaicsPeakHMM( hmmHM, FDR = 0.05, decoding="posterior",
  thres=10, parallel=TRUE, nCore=8 )

peakHM <- extractReads( peakHM,
  chipFile=system.file( file.path("extdata","wgEncodeBroadHistoneGm12878H3k4me3StdAlnRep1_chr22_sorted.bam"), package="mosaicsExample"),
  chipFileFormat="bam", chipPET=FALSE, chipFragLen=200,
  controlFile=system.file( file.path("extdata","wgEncodeBroadHistoneGm12878ControlStdAlnRep1_chr22_sorted.bam"), package="mosaicsExample"), 
  controlFileFormat="bam", controlPET=FALSE, controlFragLen=200, parallel=TRUE, nCore=8 )
peakHM

peakHM <- findSummit( peakHM, parallel=TRUE, nCore=8 )
head(print(peakHM))
plot( peakHM, filename="~/peakplot_HM.pdf" )

peakHM <- adjustBoundary( peakHM, parallel=TRUE, nCore=8 )
peakHM
head(print(peakHM))

peakHM <- filterPeak( peakHM, parallel=TRUE, nCore=8 )
peakHM
head(print(peakHM))

export( peakHM, type = "txt", filename = "./peakHM.txt" )
export( peakHM, type = "bed", filename = "./peakHM.bed" )
export( peakHM, type = "gff", filename = "./peakHM.gff" )
export( peakHM, type = "narrowPeak", filename = "./peakHM.narrowPeak" )
export( peakHM, type = "broadPeak", filename = "./peakHM.broadPeak" )
}
}
\keyword{models}
\keyword{methods}
