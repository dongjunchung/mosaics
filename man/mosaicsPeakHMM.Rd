\name{mosaicsPeakHMM}
\alias{mosaicsPeakHMM}
\alias{mosaicsPeakHMM,MosaicsHMM-method}
\title{
Call broad peaks using fitted MOSAiCS-HMM model 
}
\description{
Call broad peaks using \code{MosaicsHMM} class object, which is a fitted MOSAiCS-HMM model.
}
\usage{
mosaicsPeakHMM( object, ... )
\S4method{mosaicsPeakHMM}{MosaicsHMM}( object, FDR=0.05, decoding="viterbi",
    binsize=NA, maxgap=0, minsize=0, thres=0, 
    parallel=FALSE, nCore=8 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{Object of class \code{MosaicsHMM},
        a fitted MOSAiCS model obtained using function \code{mosaicsFitHMM}. }
    \item{FDR}{False discovery rate. Default is 0.05. Not relevant when \code{decoding="viterbi"}. }  
    \item{decoding}{Approach to determine the undelying state.
        Possible values are "viterbi" (Viterbi algorithm) and
        "posterior" (posterior decoding). Default is "viterbi". }  
    \item{binsize}{Size of each bin. Value should be positive integer.
        If \code{binsize=NA}, \code{mosaicsPeakHMM} function calcuates the value from data.
        Default is \code{NA}. }        
    \item{maxgap}{Initial nearby peaks are merged if the distance (in bp) between them is less than \code{maxgap}. Default is 0. }     
    \item{minsize}{An initial peak is removed if its width is narrower than \code{minsize}. Default is 0. }
    \item{thres}{A bin within initial peak is removed if its ChIP tag counts are less than \code{thres}. Default is 0. }  
    \item{parallel}{Utilize multiple CPUs for parallel computing 
        using \code{"parallel"} package?
        Possible values are \code{TRUE} (utilize multiple CPUs) 
        or \code{FALSE} (do not utilize multiple CPUs).
        Default is \code{FALSE} (do not utilize multiple CPUs). }
    \item{nCore}{Number of CPUs when parallel computing is utilized. }
    \item{...}{ Other parameters to be passed through to generic \code{mosaicsHMM}.}
}
\details{
\code{mosaicsFitHMM} and \code{mosaicsPeakHMM} are developed to identify broad peaks such as histone modifications,
using Hidden Markov Model (HMM) approach, as proposed in Chung et al. (2014).
If you are interested in identifying narrow peaks such as transcription factor binding sites,
please use \code{mosaicsPeak} instead of \code{mosaicsFitHMM} and \code{mosaicsPeakHMM}.

\code{maxgap}, \code{minsize}, and \code{thres} are for refining initial peaks called
using specified \code{decoding} (and \code{FDR} if \code{decoding="posterior"}).
If you use a bin size shorter than the average fragment length of the experiment, 
we recommend to set \code{maxgap} to the average fragment length 
and \code{minsize} to the bin size.
If you set the bin size to the average fragment length or if bin size is larger than the average fragment length,
set \code{maxgap} to the average fragment length and
\code{minsize} to a value smaller than the average fragment length. See the vignette for further details.

Parallel computing can be utilized for faster computing
if \code{parallel=TRUE} and \code{parallel} package is loaded.
\code{nCore} determines number of CPUs used for parallel computing.
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
\author{ Dongjun Chung, Pei Fen Kuan, Sunduz Keles}
\seealso{
\code{\link{mosaicsFit}}, \code{\link{mosaicsFitHMM}}, 
\code{\linkS4class{MosaicsFit}}, \code{\linkS4class{MosaicsHMM}}, \code{\linkS4class{MosaicsPeak}}.
}
\examples{
\dontrun{
library(mosaicsExample)
data(exampleFit)

exampleFitHMM <- mosaicsFitHMM( exampleFit, signalModel = "2S", init.FDR = 0.05 )
examplePeakHMM1 <- mosaicsPeakHMM( exampleFitHMM, decoding="viterbi" )
examplePeakHMM2 <- mosaicsPeakHMM( exampleFitHMM, FDR = 0.05, decoding="posterior" )
}
}
\keyword{models}
\keyword{methods}
