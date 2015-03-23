\name{MosaicsPeakParam-class}
\Rdversion{1.1}
\docType{class}
\alias{MosaicsPeakParam-class}

\title{Class "MosaicsPeakParam" }
\description{
This class represents parameters for mosaicsPeak function.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MosaicsPeakParam", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{analysisType}:}{Object of class \code{"character"}, 
        representing type of data analysis. }
    \item{\code{signalModel}:}{Object of class \code{"character"}, 
        representing signal model for peak calling. }
    \item{\code{FDR}:}{Object of class \code{"numeric"}, 
        representing specified false discovery rate. }
    \item{\code{maxgap}:}{Object of class \code{"numeric"}, 
        representing maximum gap size to combine two neighboring initial peaks. }
    \item{\code{minsize}:}{Object of class \code{"numeric"}, 
        representing minimum peak size to be remained in the final peak list. }
    \item{\code{thres}:}{Object of class \code{"numeric"}, 
        representing threshold of ChIP tag count to filter initial peaks. }
    \item{\code{decoding}:}{Object of class \code{"character"}, 
        representing deconding approach for MOSAiCS-HMM peak calling. }
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
\code{\link{mosaicsPeak}}, \code{\link{mosaicsPeakHMM}},
\code{\linkS4class{MosaicsPeak}}.
}
\examples{
showClass("MosaicsPeakParam")
}
\keyword{classes}
