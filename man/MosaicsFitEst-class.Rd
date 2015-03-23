\name{MosaicsFitEst-class}
\Rdversion{1.1}
\docType{class}
\alias{MosaicsFitEst-class}

\title{Class "MosaicsFitEst" }
\description{
This class represents estimates of mosaicsFit function.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MosaicsFitEst", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{pi0}:}{Object of class \code{"numeric"}, 
        representing mixing proportion of background component. }
    \item{\code{a}:}{Object of class \code{"numeric"}, 
        representing parameter for background component. }
    \item{\code{betaEst}:}{Object of class \code{"numeric"}, 
        representing parameter for background component (coefficient estimates). }
    \item{\code{muEst}:}{Object of class \code{"numeric"}, 
        representing parameter for background component. }
    \item{\code{pNfit}:}{Object of class \code{"list"}, 
        representing background model fit. }
    \item{\code{b}:}{Object of class \code{"numeric"}, 
        representing parameter for one-signal-component model. }
    \item{\code{c}:}{Object of class \code{"numeric"}, 
        representing parameter for one-signal-component model. }
    \item{\code{p1}:}{Object of class \code{"numeric"}, 
        representing parameter for two-signal-component model (mixing proportion of signal components). }
    \item{\code{b1}:}{Object of class \code{"numeric"}, 
        representing parameter for two-signal-component model (the first signal component). }
    \item{\code{c1}:}{Object of class \code{"numeric"}, 
        representing parameter for two-signal-component model (the first signal component). }
    \item{\code{b2}:}{Object of class \code{"numeric"}, 
        representing parameter for two-signal-component model (the second signal component). }
    \item{\code{c2}:}{Object of class \code{"numeric"}, 
        representing parameter for two-signal-component model (the second signal component). }
    \item{\code{inputTrunc}:}{Object of class \code{"numeric"}, 
        representing parameter for input-only analysis (threshold for maximum tag count for matched control sample. }
    \item{\code{analysisType}:}{Object of class \code{"character"}, 
        representing type of data analysis. }
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
showClass("MosaicsFitEst")
}
\keyword{classes}
