\name{MosaicsFitParam-class}
\Rdversion{1.1}
\docType{class}
\alias{MosaicsFitParam-class}

\title{Class "MosaicsFitParam" }
\description{
This class represents parameters for mosaicsFit function.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MosaicsFitParam", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{k}:}{Object of class \code{"numeric"}, 
        representing a parameter for mosaicsFit function. }
    \item{\code{meanThres}:}{Object of class \code{"numeric"}, 
        representing a parameter for mosaicsFit function. }
    \item{\code{s}:}{Object of class \code{"numeric"}, 
        representing a parameter for mosaicsFit function. }
    \item{\code{d}:}{Object of class \code{"numeric"}, 
        representing a parameter for mosaicsFit function. }
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
\code{\link{mosaicsFit}}.
}
\examples{
showClass("MosaicsFitParam")
}
\keyword{classes}
