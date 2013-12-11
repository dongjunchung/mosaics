\name{MosaicsHMM-class}
\Rdversion{1.1}
\docType{class}
\alias{MosaicsHMM-class}
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
    	representing the data used to fit MOSAiCS-HMM model. }
    \item{\code{peakParam}:}{Object of class \code{"MosaicsPeakParam"},
        representing parameters for peak calling. }
    \item{\code{binsize}:}{Object of class \code{"numeric"}, 
    		representing size of a bin. }
    \item{\code{nRatio}:}{Object of class \code{"numeric"}, 
    		representing ratio of sequencing depth of ChIP vs. control. }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "MosaicsHMM")}: provide brief summary of the object. }
     }
}
\references{
Kuan, PF, D Chung, G Pan, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.
}
\author{ Dongjun Chung, Pei Fen Kuan, Sunduz Keles }
\seealso{
\code{\link{mosaicsFitHMM}}, \code{\link{mosaicsPeakHMM}}.
}
\examples{
showClass("MosaicsHMM")
\dontrun{
library(mosaicsExample)
data(exampleFit)

exampleFitHMM <- mosaicsFitHMM( exampleFit, signalModel = "2S", init.FDR = 0.05 )
exampleFitHMM

examplePeakHMM1 <- mosaicsPeakHMM( exampleFitHMM, decoding="viterbi" )
examplePeakHMM2 <- mosaicsPeakHMM( exampleFitHMM, FDR = 0.05, decoding="posterior" )
}
}
\keyword{classes}
