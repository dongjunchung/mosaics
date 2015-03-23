\name{MosaicsFit-class}
\Rdversion{1.1}
\docType{class}
\alias{MosaicsFit-class}
\alias{plot,MosaicsFit,ANY-method}
\alias{print,MosaicsFit-method}
\alias{show,MosaicsFit-method}

\title{Class "MosaicsFit" }
\description{
This class represents MOSAiCS model fit.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MosaicsFit", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{mosaicsEst}:}{Object of class \code{"MosaicsFitEst"},
        representing estimates of MOSAiCS model fit. }
    \item{\code{mosaicsParam}:}{Object of class \code{"MosaicsFitParam"},
        representing tuning parameters for fitting MOSAiCS model. }
    \item{\code{chrID}:}{Object of class \code{"character"},
        a vector of chromosome IDs. }
    \item{\code{coord}:}{Object of class \code{"numeric"},
        a vector of genomic coordinates. }
    \item{\code{tagCount}:}{Object of class \code{"numeric"},
        a vector of tag counts of ChIP sample. }
    \item{\code{mappability}:}{Object of class \code{"numeric"},
        a vector of mappability score. }
    \item{\code{gcContent}:}{Object of class \code{"numeric"},
        a vector of GC content score. }
    \item{\code{input}:}{Object of class \code{"numeric"},
        a vector of tag counts of matched control sample. }
    \item{\code{bic1S}:}{Object of class \code{"numeric"},
        Bayesian Information Criterion (BIC) value of one-signal-component model. }
    \item{\code{bic2S}:}{Object of class \code{"numeric"},
        Bayesian Information Criterion (BIC) value of two-signal-component model. }
  }
}
\section{Methods}{
  \describe{
    \item{estimates}{\code{signature(object = "MosaicsFit")}: 
        extract estimates from MOSAiCS model fit. }
    \item{mosaicsPeak}{\code{signature(object = "MosaicsFit")}: 
        call peaks using MOSAiCS model fit. }
    \item{plot}{\code{signature(x = "MosaicsFit", y = "missing")}: 
        draw Goodness of Fit (GOF) plot. }
    \item{print}{\code{signature(x = "MosaicsFit")}: 
        (not supported yet) }
    \item{show}{\code{signature(object = "MosaicsFit")}:
        provide brief summary of the object. }
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
\code{\link{mosaicsFit}}, \code{\link{mosaicsPeak}}, \code{\link{mosaicsFitHMM}}, \code{\link{estimates}}.
}
\examples{
showClass("MosaicsFit")
\dontrun{
library(mosaicsExample)

data(exampleBinData)
exampleFit <- mosaicsFit( exampleBinData, analysisType="IO" )

exampleFit
plot(exampleFit)
estimates(exampleFit)

examplePeak <- mosaicsPeak( exampleFit, signalModel = "2S", FDR = 0.05 )
}
}
\keyword{classes}
