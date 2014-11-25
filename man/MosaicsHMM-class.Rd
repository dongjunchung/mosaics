\name{MosaicsHMM-class}
\Rdversion{1.1}
\docType{class}
\alias{MosaicsHMM-class}
\alias{estimates,MosaicsHMM-method}
\alias{plot,MosaicsHMM,missing-method}
\alias{print,MosaicsHMM-method}
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
    \item{\code{mosaicsFit}:}{Object of class \code{"MosaicsFit"}, 
    	representing the fitted MOSAiCS model. }
    \item{\code{init}:}{Object of class \code{"character"}, 
    		representing the approach to initialize MOSAiCS-HMM. }
    \item{\code{initPiMat}:}{Object of class \code{"numeric"}, 
    		representing initial transition matrix. }
    \item{\code{peakParam}:}{Object of class \code{"MosaicsPeakParam"},
        representing parameters for peak calling. }
    \item{\code{binsize}:}{Object of class \code{"numeric"}, 
    		representing size of a bin. }
    \item{\code{nRatio}:}{Object of class \code{"numeric"}, 
    		representing ratio of sequencing depth of ChIP vs. control. }
    \item{\code{bicMosaics}:}{Object of class \code{"numeric"}, 
    		representing the BIC value of MOSAiCS fit. }
    \item{\code{bicMosaicsHMM}:}{Object of class \code{"numeric"}, 
    		representing the BIC value of MOSAiCS-HMM fit. }
  }
}
\section{Methods}{
  \describe{
    \item{estimates}{\code{signature(object = "MosaicsHMM")}: 
        extract estimates from MOSAiCS-HMM model fit. }=
    \item{plot}{\code{signature(x = "MosaicsHMM", y = "missing", seed=12345, parallel=FALSE, nCore=8 )}: 
        draw Goodness of Fit (GOF) plot.
		You can specify random seed in \code{seed}.
		If \code{parallel=TRUE}, parallel computing is utilized to simulate data,
		where \code{nCore} indicates CPUs used for parallel computing. }
    \item{print}{\code{signature(x = "MosaicsHMM")}: 
        (not supported yet) }
    \item{show}{\code{signature(object = "MosaicsHMM")}: provide brief summary of the object. }
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
\code{\link{mosaicsFitHMM}}, \code{\link{mosaicsPeakHMM}}.
}
\examples{
showClass("MosaicsHMM")
\dontrun{
library(mosaicsExample)
data(exampleFit)

exampleFitHMM <- mosaicsFitHMM( 
	exampleFit, signalModel = "2S", init = "mosaics",
	init.FDR=0.05, init.maxgap=200, init.minsize=50, init.thres=10,
	parallel=TRUE, nCore=8 )
exampleFitHMM
plot( exampleFitHMM, parallel=TRUE, nCore=8 )
estimates(exampleFitHMM)

examplePeakHMM1 <- mosaicsPeakHMM( exampleFitHMM, decoding="viterbi" )
examplePeakHMM2 <- mosaicsPeakHMM( exampleFitHMM, FDR = 0.05, decoding="posterior" )
}
}
\keyword{classes}
