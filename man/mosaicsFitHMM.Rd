\name{mosaicsFitHMM}
\alias{mosaicsFitHMM}
\alias{mosaicsFitHMM,MosaicsFit-method}
\title{
Fit MOSAiCS-HMM model
}
\description{
Fit MOSAiCS-HMM model.
}
\usage{
mosaicsFitHMM( object, ... )
\S4method{mosaicsFitHMM}{MosaicsFit}( object, signalModel="2S", binsize=NA,
	init="mosaics", init.FDR=0.05, 
	init.maxgap=200, init.minsize=50, init.thres=10, init.piMat=as.matrix(NA),
	max.iter=100, eps=1e-20, parallel=FALSE, nCore=8 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{Object of class \code{MosaicsFit},
        a fitted MOSAiCS model obtained using function \code{mosaicsFit}. }
    \item{signalModel}{Signal model.
        Possible values are "1S" (one-signal-component model) and
        "2S" (two-signal-component model). Default is "2S". }
    \item{binsize}{Size of each bin. Value should be positive integer.
        If \code{binsize=NA}, \code{mosaicsFitHMM} function calcuates the value from data.
        Default is \code{NA}. }        
    \item{init}{Approach to initialize MOSAiCS-HMM. Possible values are \code{"mosaics"} (use MOSAiCS peak calling results for initialization) or \code{"specify"} (explicitly specify transition matrix). Default is \code{"mosaics"}.}
	\item{init.FDR}{Parameter for the MOSAiCS-HMM initialization. False discovery rate. Default is 0.05. Related only if \code{init="mosaics"}.}  
    \item{init.maxgap}{Parameter for the MOSAiCS-HMM initialization. Initial nearby peaks are merged if the distance (in bp) between them is less than \code{init.maxgap}. Default is 200. Related only if \code{init="mosaics"}.}     
    \item{init.minsize}{Parameter for the MOSAiCS-HMM initialization. An initial peak is removed if its width is narrower than \code{init.minsize}. Default is 50. Related only if \code{init="mosaics"}.}
    \item{init.thres}{Parameter for the MOSAiCS-HMM initialization. A bin within initial peak is removed if its ChIP tag counts are less than \code{init.thres}. Default is 10. Related only if \code{init="mosaics"}.}  
	\item{init.piMat}{Initial value for transition matrix. The first rows/columns correspond to the non-binding state while the second rows/columns correspond to the binding state. Related only if \code{init="specify"}. If \code{init="specify"} but \code{init.piMat} is not specified, \code{mosaicsFitHMM()} uses its default for the MOSAiCS-HMM initialization.}
    \item{max.iter}{Number of iterations for fitting MOSAiCS-HMM. Default is 100.}     
    \item{eps}{Criterion to stop iterations for fitting MOSAiCS-HMM. Default is 1e-20.}
    \item{parallel}{Utilize multiple CPUs for parallel computing 
        using \code{"parallel"} package?
        Possible values are \code{TRUE} (utilize multiple CPUs) 
        or \code{FALSE} (do not utilize multiple CPUs).
        Default is \code{FALSE} (do not utilize multiple CPUs). }
    \item{nCore}{Number of CPUs when parallel computing is utilized. }
    \item{...}{ Other parameters to be passed through to generic \code{mosaicsFitHMM}.}
}
\details{
\code{mosaicsFitHMM} and \code{mosaicsPeakHMM} are developed to identify broad peaks such as histone modifications,
using Hidden Markov Model (HMM) approach, as proposed in Chung et al. (2014). 
If you are interested in identifying narrow peaks such as transcription factor binding sites,
please use \code{mosaicsPeak} instead of \code{mosaicsFitHMM} and \code{mosaicsPeakHMM}.

When peaks are called, proper signal model needs to be specified.
The optimal choice for the number of signal components depends on the characteristics of ChIP-seq data.
In order to support users in the choice of optimal signal model, 
Bayesian Information Criterion (BIC) values and Goodness of Fit (GOF) plot are provided
for the fitted MOSAiCS model.
BIC values and GOF plot can be obtained by applying \code{show} and \code{plot} methods,
respectively, to the \code{MosaicsFit} class object, which is a fitted MOSAiCS model.

\code{init.FDR}, \code{init.maxgap}, \code{init.minsize}, and \code{init.thres} are the parameters for MOSAiCS-HMM initialization when MOSAiCS peak calling results are used for initialization (\code{init="mosaics"}). If user specifies transition matrix (\code{init="specify"}), only \code{init.piMat} is used for initialization.
If you use a bin size shorter than the average fragment length of the experiment, 
we recommend to set \code{init.maxgap} to the average fragment length 
and \code{init.minsize} to the bin size.
If you set the bin size to the average fragment length or if bin size is larger than the average fragment length,
set \code{init.maxgap} to the average fragment length and
\code{init.minsize} to a value smaller than the average fragment length. See the vignette for further details.

Parallel computing can be utilized for faster computing
if \code{parallel=TRUE} and \code{parallel} package is loaded.
\code{nCore} determines number of CPUs used for parallel computing.
}
\value{
Construct \code{MosaicsHMM} class object.
}
\references{
Kuan, PF, D Chung, G Pan, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.

Chung, D, Zhang Q, and Keles S (2014), "MOSAiCS-HMM: A model-based approach for detecting regions of histone modifications from ChIP-seq data", Datta S and Nettleton D (eds.), \emph{Statistical Analysis of Next Generation Sequencing Data}, Springer.
}
\author{ Dongjun Chung, Pei Fen Kuan, Sunduz Keles }
\seealso{
\code{\link{mosaicsFit}}, \code{\link{mosaicsPeakHMM}}, 
\code{\linkS4class{MosaicsFit}}, \code{\linkS4class{MosaicsHMM}}, \code{\linkS4class{MosaicsPeak}}.
}
\examples{
\dontrun{
library(mosaicsExample)
data(exampleFit)

exampleFitHMM <- mosaicsFitHMM( 
	exampleFit, signalModel = "2S", init = "mosaics",
	init.FDR=0.05, init.maxgap=200, init.minsize=50, init.thres=10,
	parallel=TRUE, nCore=8 )
examplePeakHMM1 <- mosaicsPeakHMM( exampleFitHMM, decoding = "viterbi" )
examplePeakHMM2 <- mosaicsPeakHMM( exampleFitHMM, FDR = 0.05, decoding = "posterior" )
}
}
\keyword{models}
\keyword{methods}
