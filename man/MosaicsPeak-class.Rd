\name{MosaicsPeak-class}
\Rdversion{1.1}
\docType{class}
\alias{MosaicsPeak-class}
\alias{bdBin,MosaicsPeak-method}
\alias{empFDR,MosaicsPeak-method}
%\alias{export,MosaicsPeak-method}
\alias{print,MosaicsPeak-method}
\alias{show,MosaicsPeak-method}

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
  }
}
\section{Methods}{
  \describe{
    \item{export}{\code{signature(object = "MosaicsPeak")}: export peak list into text files. }
    \item{print}{\code{signature(x = "MosaicsPeak")}: return peak list in data frame format. }
    \item{show}{\code{signature(object = "MosaicsPeak")}: provide brief summary of the object. }
     }
}
\references{
Kuan, PF, D Chung, G Pan, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.
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
