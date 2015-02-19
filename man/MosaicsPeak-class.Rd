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
	 \item{coverage}{\code{signature(object = "MosaicsPeak")}: export a list of coverage matrices of ChIP and matched control samples, where each matrix consists of two columns, genomic coordinates and read count.
	 Each element of this list is a list of length 2, which are matrices for each of ChIP and matched control samples.
   Elements of this list are sorted to match the rows of peak list.
	 Users can use \code{coverage} method only after \code{extractReads} method is applied to the \code{MosaicsPeak} object. }
  \item{read}{\code{signature(object = "MosaicsPeak")}: export a list of read-level data of ChIP and matched control samples, where each element is \code{GenomicRanges} class.
	 Each element of this list is a list of length 2, which are \code{GenomicRanges} objects for each of ChIP and matched control samples.
	 Elements of this list are sorted to match the rows of peak list.
	 Users can use \code{coverage} method only after \code{extractReads} method is applied to the \code{MosaicsPeak} object. }
	 \item{seqDepth}{\code{signature(object = "MosaicsPeak")}: export a vector of sequencing depth, where the first and second elements are for ChIP and matched control sample, respectively. If matched control sample is not provided, the second element is NA.
	 Users can use \code{coverage} method only after \code{extractReads} method is applied to the \code{MosaicsPeak} object. }
}
}
\references{
Kuan, PF, D Chung, G Pan, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.  

Chung, D, Zhang Q, and Keles S (2014), "MOSAiCS-HMM: A model-based approach for detecting regions of histone modifications from ChIP-seq data", Datta S and Nettleton D (eds.), \emph{Statistical Analysis of Next Generation Sequencing Data}, Springer.
}
\author{ Dongjun Chung, Pei Fen Kuan, Sunduz Keles }
\seealso{
\code{\link{mosaicsPeak}}, \code{\link{export}}.
}
\examples{
showClass("MosaicsPeak")
\dontrun{
library(mosaicsExample)
data(exampleFit)
examplePeak <- mosaicsPeak( exampleFit, signalModel = "2S", FDR = 0.05 )

examplePeak
print(examplePeak)[1:10, ]
export( examplePeak, type = "txt", filename = "./TSpeakList.txt" )
export( examplePeak, type = "bed", filename = "./TSpeakList.bed" )
export( examplePeak, type = "gff", filename = "./TSpeakList.gff" )
}
}
\keyword{classes}
