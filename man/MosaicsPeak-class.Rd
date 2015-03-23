\name{MosaicsPeak-class}
\Rdversion{1.1}
\docType{class}
\alias{MosaicsPeak-class}
\alias{bdBin,MosaicsPeak-method}
\alias{empFDR,MosaicsPeak-method}
%\alias{export,MosaicsPeak-method}
\alias{print,MosaicsPeak-method}
\alias{plot,MosaicsPeak,missing-method}
\alias{show,MosaicsPeak-method}
\alias{coverage,MosaicsPeak-method}
\alias{read,MosaicsPeak-method}
\alias{seqDepth,MosaicsPeak-method}
\alias{bdBin}
\alias{empFDR}
\alias{coverage}
\alias{read}
\alias{seqDepth}

\title{Class "MosaicsPeak" }
\description{
This class represents peak calling results.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MosaicsPeak", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{peakList}:}{Object of class \code{"data.frame"}, representing peak list. }
    \item{\code{peakParam}:}{Object of class \code{"MosaicsPeakParam"},
        representing parameters for peak calling. }
    \item{\code{bdBin}:}{Object of class \code{"numeric"}, 
    		representing a vector of bounded bins. }
    \item{\code{empFDR}:}{Object of class \code{"numeric"}, 
    		representing empirical FDR. }
    \item{\code{postProb}:}{Object of class \code{"numeric"}, 
      	representing posterior probability that a bin belongs to background. }
    \item{\code{tagLoaded}:}{Object of class \code{"logical"}, 
    		representing whether read-level data is loaded. }
    \item{\code{tagData}:}{Object of class \code{"TagData"}, 
    		representing read-level data. }
  }
}
\section{Methods}{
  \describe{
    \item{export}{\code{signature(object = "MosaicsPeak")}: export peak list into text files. }
    \item{print}{\code{signature(x = "MosaicsPeak")}: return peak list in data frame format. }
    \item{plot}{\code{signature( x = "MosaicsPeak", y = "missing", filename=NA, peakNum=NA )}: 
      plot ChIP profile in each peak region.
      If file name is specified in \code{filename}, plots are exported to a PDF file named \code{filename}.
      Oterwise, ChIP profiles are plotted to standard output.
      If users are interested in specific peak regions, users can specify them as a vector in \code{peakNum},
      where elements indicate which rows peak regions of interest are located.
      If \code{peakNum} is specified, only ChIP profiles in specified peak regions are plotted.
      Otherwise, ChIP profiles for all the peak regions are plotted.}
    \item{show}{\code{signature(object = "MosaicsPeak")}: provide brief summary of the object. }
    \item{empFDR}{\code{signature(object = "MosaicsPeak")}: return estimated empirical false discovery rate. }
    \item{bdBin}{\code{signature(object = "MosaicsPeak")}: return a vector of bin-level binary indicator of peak calling,
      where each element is 1 if the corresponding bin belongs to a peak and zero otherwise. }
	 \item{coverage}{\code{signature(object = "MosaicsPeak")}: export a list of coverage matrices of ChIP and matched control samples, where each matrix consists of two columns, genomic coordinates and read count.
	 Each element of this list is a list of length 2, which are matrices for each of ChIP and matched control samples.
   Elements of this list are sorted to match the rows of peak list.
	 Users can use \code{coverage} method only after \code{extractReads} method is applied to the \code{MosaicsPeak} object. }
  \item{read}{\code{signature(object = "MosaicsPeak")}: return a list of read-level data of ChIP and matched control samples, where each element is \code{GenomicRanges} class.
	 Each element of this list is a list of length 2, which are \code{GenomicRanges} objects for each of ChIP and matched control samples.
	 Elements of this list are sorted to match the rows of peak list.
	 Users can use \code{coverage} method only after \code{extractReads} method is applied to the \code{MosaicsPeak} object. }
	 \item{seqDepth}{\code{signature(object = "MosaicsPeak")}: return a vector of sequencing depth, where the first and second elements are for ChIP and matched control sample, respectively. If matched control sample is not provided, the second element is NA.
	 Users can use \code{coverage} method only after \code{extractReads} method is applied to the \code{MosaicsPeak} object. }
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
\code{\link{mosaicsPeak}}, \code{\link{mosaicsPeakHMM}}, \code{\link{export}}, 
\code{\link{extractReads}}, \code{\link{findSummit}}, \code{\link{adjustBoundary}}, \code{\link{filterPeak}}.
}
\examples{
showClass("MosaicsPeak")
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
\keyword{classes}
