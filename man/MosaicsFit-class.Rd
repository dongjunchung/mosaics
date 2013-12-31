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
Kuan, PF, D Chung, JA Thomson, R Stewart, and S Keles (2010), "A Statistical Framework
for the Analysis of ChIP-Seq Data", To appear in \emph{Journal of the American Statistical Association} (\url{http://pubs.amstat.org/doi/abs/10.1198/jasa.2011.ap09706}).
}
\author{ Dongjun Chung, Pei Fen Kuan, Sunduz Keles }
\seealso{
\code{\link{mosaicsFit}}, \code{\link{mosaicsPeak}}, \code{\link{estimates}}.
}
\examples{
showClass("MosaicsFit")
\dontrun{
library(mosaicsExample)
data(exampleFit)

exampleFit
plot(exampleFit)
estimates(exampleFit)

examplePeak <- mosaicsPeak( exampleFit, signalModel = "2S", FDR = 0.05 )
}
}
\keyword{classes}
