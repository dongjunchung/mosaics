\name{TagData-class}
\Rdversion{1.1}
\docType{class}
\alias{TagData-class}

\title{Class "TagData" }
\description{
This class represents read-level ChIP-seq data.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TagData", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{coverage}:}{Object of class \code{"list"}, 
        a list of read coverage, where each element is a list including ChIP and matched control samples. }
    \item{\code{numReads}:}{Object of class \code{"numeric"}, 
        a matrix of number of reads in peak regions, where rows correspond to peak regions and columns correspond to ChIP and matched control samples. }
    \item{\code{read}:}{Object of class \code{"list"}, 
        a list of read-level data, where each element is a list including ChIP and matched control samples. }
    \item{\code{keepReads}:}{Object of class \code{"logical"}, 
        representing whether to keep read-level data. }
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
\code{\link{extractReads}}, \code{\link{findSummit}}, \code{\link{filterPeak}}, 
\code{\linkS4class{MosaicsPeak}}.
}
\examples{
showClass("TagData")
}
\keyword{classes}
