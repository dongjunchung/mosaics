\name{MosaicsHMM-class}
\Rdversion{1.1}
\docType{class}
\alias{MosaicsHMM-class}
\alias{estimates,MosaicsHMM-method}
\alias{plot,MosaicsHMM,missing-method}
\alias{print,MosaicsHMM-method}
\alias{show,MosaicsHMM-method}

\title{Class "MosaicsHMM" }
\description{
This class represents MOSAiCS-HMM model fit.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MosaicsHMM", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{HMMfit}:}{Object of class \code{"list"}, 
    	representing the fitted MOSAiCS-HMM model. }
    \item{\code{inputdata}:}{Object of class \code{"list"}, 
      representing the bin-level data. }
    \item{\code{mosaicsEst}:}{Object of class \code{"MosaicsFitEst"},
        representing estimates of MOSAiCS model fit. }
    \item{\code{init}:}{Object of class \code{"character"}, 
    		representing the approach to initialize MOSAiCS-HMM. }
    \item{\code{initPiMat}:}{Object of class \code{"numeric"}, 
    		representing initial transition matrix. }
    \item{\code{peakParam}:}{Object of class \code{"MosaicsPeakParam"},
        representing parameters for peak calling. }
    \item{\code{binsize}:}{Object of class \code{"numeric"}, 
    		representing size of a bin. }
    \item{\code{nRatio}:}{Object of class \code{"numeric"}, 
    		representing ratio of sequencing depth of ChIP vs. control. }
    \item{\code{bicMosaics}:}{Object of class \code{"numeric"}, 
    		representing the BIC value of MOSAiCS fit. }
    \item{\code{bicMosaicsHMM}:}{Object of class \code{"numeric"}, 
    		representing the BIC value of MOSAiCS-HMM fit. }
  }
}
\section{Methods}{
  \describe{
    \item{estimates}{\code{signature(object = "MosaicsHMM")}: 
        extract estimates from MOSAiCS-HMM model fit. }=
    \item{plot}{\code{signature(x = "MosaicsHMM", y = "missing", seed=12345, parallel=FALSE, nCore=8 )}: 
        draw Goodness of Fit (GOF) plot.
		You can specify random seed in \code{seed}.
		If \code{parallel=TRUE}, parallel computing is utilized to simulate data,
		where \code{nCore} indicates CPUs used for parallel computing. }
    \item{print}{\code{signature(x = "MosaicsHMM")}: 
        (not supported yet) }
    \item{show}{\code{signature(object = "MosaicsHMM")}: provide brief summary of the object. }
     }
}
\references{
Kuan, PF, D Chung, G Pan, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.

Chung, D, Zhang Q, and Keles S (2014), "MOSAiCS-HMM: A model-based approach for detecting regions of histone modifications from ChIP-seq data", Datta S and Nettleton D (eds.), \emph{Statistical Analysis of Next Generation Sequencing Data}, Springer.
}
\author{ Dongjun Chung, Pei Fen Kuan, Rene Welch, Sunduz Keles }
\seealso{
\code{\link{mosaicsFitHMM}}, \code{\link{mosaicsPeakHMM}}.
}
\examples{
showClass("MosaicsHMM")
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
\keyword{classes}
