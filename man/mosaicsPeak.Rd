\name{mosaicsPeak}
\alias{mosaicsPeak}
\alias{mosaicsPeak,MosaicsFit-method}
\title{
Call peaks using fitted MOSAiCS model
}
\description{
Call peaks using \code{MosaicsFit} class object, which is a fitted MOSAiCS model.
}
\usage{
mosaicsPeak( object, ... )
\S4method{mosaicsPeak}{MosaicsFit}( object, signalModel="2S", FDR=0.05, 
    binsize=NA, maxgap=200, minsize=50, thres=10 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{Object of class \code{MosaicsFit},
        a fitted MOSAiCS model obtained using function \code{mosaicsFit}. }
    \item{signalModel}{Signal model.
        Possible values are "1S" (one-signal-component model) and
        "2S" (two-signal-component model). Default is "2S". }
    \item{FDR}{False discovery rate. Default is 0.05. }    
    \item{binsize}{Size of each bin. Value should be positive integer.
        If \code{binsize=NA}, \code{mosaicsPeak} function calcuates the value from data.
        Default is \code{NA}. }     
    \item{maxgap}{Initial nearby peaks are merged if the distance (in bp) between them is less than \code{maxgap}. Default is 200. }     
    \item{minsize}{An initial peak is removed if its width is narrower than \code{minsize}. Default is 50. }
    \item{thres}{A bin within initial peak is removed if its ChIP tag counts are less than \code{thres}. Default is 10. }  
    \item{...}{ Other parameters to be passed through to generic \code{mosaicsPeak}.}
}
\details{
\code{mosaicsPeak} is developed to identify narrow peaks such as transcription factor binding sites. 
If you are interested in identifying broad peaks such as histone modifications,
please use \code{mosaicsFitHMM} and \code{mosaicsPeakHMM} instead of \code{mosaicsPeak}.

When peaks are called, proper signal model needs to be specified.
The optimal choice for the number of signal components depends on the characteristics of ChIP-seq data.
In order to support users in the choice of optimal signal model, 
Bayesian Information Criterion (BIC) values and Goodness of Fit (GOF) plot are provided
for the fitted MOSAiCS model.
BIC values and GOF plot can be obtained by applying \code{show} and \code{plot} methods,
respectively, to the \code{MosaicsFit} class object, which is a fitted MOSAiCS model.

\code{maxgap}, \code{minsize}, and \code{thres} are for refining initial peaks called
using specified \code{signalModel} and \code{FDR}.
If you use a bin size shorter than the average fragment length of the experiment, 
we recommend to set \code{maxgap} to the average fragment length 
and \code{minsize} to the bin size.
If you set the bin size to the average fragment length or if bin size is larger than the average fragment length,
set \code{maxgap} to the average fragment length and
\code{minsize} to a value smaller than the average fragment length. See the vignette for further details.
}
\value{
Construct \code{MosaicsPeak} class object.
}
\references{
Kuan, PF, D Chung, G Pan, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.

Chung, D, Zhang Q, and Keles S (2014), "MOSAiCS-HMM: A model-based approach for detecting regions of histone modifications from ChIP-seq data", Datta S and Nettleton D (eds.), \emph{Statistical Analysis of Next Generation Sequencing Data}, Springer.
}
\author{ Dongjun Chung, Pei Fen Kuan, Rene Welch, Sunduz Keles }
\seealso{
\code{\link{mosaicsFit}}, 
\code{\linkS4class{MosaicsPeak}}, \code{\linkS4class{MosaicsFit}}.
}
\examples{
\dontrun{
library(mosaicsExample)
data(exampleBinData)
exampleFit <- mosaicsFit( exampleBinData, analysisType="IO" )
examplePeak <- mosaicsPeak( exampleFit, signalModel = "2S", FDR = 0.05 )
}
}
\keyword{models}
\keyword{methods}
