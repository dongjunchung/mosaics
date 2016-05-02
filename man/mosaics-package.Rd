\name{mosaics-package}
\alias{mosaics-package}
\alias{mosaics}
\docType{package}
\title{
MOSAiCS (MOdel-based one and two Sample Analysis and Inference for ChIP-Seq)
}
\description{
This package provides functions for fitting MOSAiCS,
a statistical framework to analyze one-sample or two-sample ChIP-seq data.
}
\details{
\tabular{ll}{
Package: \tab mosaics\cr
Type: \tab Package\cr
Version: \tab 2.9.9\cr
Date: \tab 2016-03-15\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
This package contains four main classes, \code{BinData}, \code{MosaicsFit}, \code{MosaicsHMM}, and \code{MosaicsPeak},
which represent bin-level ChIP-seq data, MOSAiCS model fit, MOSAiCS-HMM model fit, and MOSAiCS peak calling results, respectively.
This package contains ten main methods, 
\code{constructBins}, \code{readBins}, \code{mosaicsFit}, \code{mosaicsPeak}, \code{mosaicsFitHMM}, \code{mosaicsPeakHMM}, \code{extractReads}, \code{findSummit}, \code{adjustBoundary}, \code{filterPeak}.
\code{constructBins} method constructs bin-level files from the aligned read file.
\code{readBins} method imports bin-level data and construct \code{BinData} class object.
\code{mosaicsFit} method fits a MOSAiCS model using \code{BinData} class object and
 constructs \code{MosaicsFit} class object.
\code{mosaicsPeak} method calls peaks using \code{MosaicsFit} class object and
 construct \code{MosaicsPeak} class object.
\code{mosaicsFitHMM} and \code{mosaicsPeakHMM} are designed to identify broad peaks
and their functions correspond to \code{mosaicsFit} and \code{mosaicsPeak}, respectively.
\code{mosaicsFitHMM} method fits a MOSAiCS-HMM model using \code{MosaicsFit} class object and
 constructs \code{MosaicsHMM} class object.
\code{mosaicsPeakHMM} method calls MOSAiCS-HMM peaks using \code{MosaicsHMM} class object and
 construct \code{MosaicsPeak} class object.
\code{extractReads}, \code{findSummit}, \code{adjustBoundary}, \code{filterPeak} methods 
postprocess MOSAiCS and MOSAiCS-HMM peaks using \code{MosaicsPeak} class object 
(incorporate read-level data, identify peak summits, adjust peak boundaries, and filter potential false positive peaks, respectively)
and construct \code{MosaicsPeak} class object.
\code{MosaicsPeak} class object can be exported as text files or transformed into data frame,
which can be used for the downstream analysis.
This package also provides methods for simple exploratory analysis.

The \code{mosaics} package companion website, \url{http://www.stat.wisc.edu/~keles/Software/mosaics/},
provides preprocessing scripts, preprocessed files for diverse reference genomes,
and easy-to-follow instructions.
We encourage questions or requests regarding \code{mosaics} package to be posted on
our Google group, \url{http://groups.google.com/group/mosaics_user_group}.
Please check the vignette for further details on the \code{mosaics} package and these websites.
}
\author{
Dongjun Chung, Pei Fen Kuan, Rene Welch, Sunduz Keles

Maintainer: Dongjun Chung <dongjun.chung@gmail.com>
}
\references{
Kuan, PF, D Chung, G Pan, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.

Chung, D, Zhang Q, and Keles S (2014), "MOSAiCS-HMM: A model-based approach for detecting regions of histone modifications from ChIP-seq data", Datta S and Nettleton D (eds.), \emph{Statistical Analysis of Next Generation Sequencing Data}, Springer.
}
\keyword{ package }
\seealso{
\code{\link{constructBins}}, \code{\link{readBins}}, \code{\link{mosaicsFit}}, \code{\link{mosaicsPeak}},
\code{\link{mosaicsFitHMM}}, \code{\link{mosaicsPeakHMM}},
\code{\link{extractReads}}, \code{\link{findSummit}}, \code{\link{adjustBoundary}}, \code{\link{filterPeak}},
\code{\linkS4class{BinData}}, \code{\linkS4class{MosaicsFit}}, \code{\linkS4class{MosaicsHMM}}, \code{\linkS4class{MosaicsPeak}}.
}
\examples{
\dontrun{
library(mosaicsExample)

# example analysis workflow for ChIP-seq data of transcription factor binding
# (STAT1 factor binding in GM12878 cell line, from ENCODE database)

generateWig( infile=system.file( file.path("extdata","wgEncodeSydhTfbsGm12878Stat1StdAlnRep1_chr22_sorted.bam"), package="mosaicsExample"), 
    fileFormat="bam", outfileLoc="~/", 
    PET=FALSE, fragLen=200, span=200, capping=0, normConst=1 )

constructBins( infile=system.file( file.path("extdata","wgEncodeSydhTfbsGm12878Stat1StdAlnRep1_chr22_sorted.bam"), package="mosaicsExample"), 
    fileFormat="bam", outfileLoc="~/", 
    PET=FALSE, fragLen=200, binSize=200, capping=0 )
constructBins( infile=system.file( file.path("extdata","wgEncodeSydhTfbsGm12878InputStdAlnRep1_chr22_sorted.bam"), package="mosaicsExample"), 
    fileFormat="bam", outfileLoc="~/", 
    PET=FALSE, fragLen=200, binSize=200, capping=0 )
    
binTFBS <- readBins( type=c("chip","input"),
    fileName=c( "~/wgEncodeSydhTfbsGm12878Stat1StdAlnRep1_chr22_sorted.bam_fragL200_bin200.txt",
    "~/wgEncodeSydhTfbsGm12878InputStdAlnRep1_chr22_sorted.bam_fragL200_bin200.txt" ) )
binTFBS
head(print(binTFBS))
plot(binTFBS)
plot( binTFBS, plotType="input" )

fitTFBS <- mosaicsFit( binTFBS, analysisType="IO", bgEst="rMOM" )
fitTFBS
plot(fitTFBS)

peakTFBS <- mosaicsPeak( fitTFBS, signalModel="2S", FDR=0.05, 
maxgap=200, minsize=50, thres=10 )
peakTFBS
head(print(peakTFBS))

peakTFBS <- extractReads( peakTFBS,
  chipFile=system.file( file.path("extdata","wgEncodeSydhTfbsGm12878Stat1StdAlnRep1_chr22_sorted.bam"), package="mosaicsExample"),
  chipFileFormat="bam", chipPET=FALSE, chipFragLen=200,
  controlFile=system.file( file.path("extdata","wgEncodeSydhTfbsGm12878InputStdAlnRep1_chr22_sorted.bam"), package="mosaicsExample"), 
  controlFileFormat="bam", controlPET=FALSE, controlFragLen=200, parallel=TRUE, nCore=8 )
peakTFBS

peakTFBS <- findSummit( peakTFBS, parallel=FALSE, nCore=8 )

export( peakTFBS, type = "txt", filename = "./peakTFBS.txt" )
export( peakTFBS, type = "bed", filename = "./peakTFBS.bed" )
export( peakTFBS, type = "gff", filename = "./peakTFBS.gff" )
export( peakTFBS, type = "narrowPeak", filename = "./peakTFBS.narrowPeak" )

# example analysis workflow for ChIP-seq data of histone modification
# (H3K4me3 modification in GM12878 cell line, from ENCODE database)

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
