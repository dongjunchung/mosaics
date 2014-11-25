\name{mosaicsFit}
\alias{mosaicsFit}
\alias{mosaicsFit,BinData-method}
\title{
Fit MOSAiCS model
}
\description{
Fit one-sample or two-sample MOSAiCS models with one signal component and two signal components.
}
\usage{
mosaicsFit( object, ... )
\S4method{mosaicsFit}{BinData}( object, analysisType="automatic", bgEst="rMOM",
    k=3, meanThres=NA, s=2, d=0.25, trans="power", truncProb=0.999, parallel=FALSE, nCore=8 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{Object of class \code{BinData},
        bin-level ChIP-seq data imported using method \code{readBins}. }
    \item{analysisType}{Analysis type.
        Possible values are "OS" (one-sample analysis),
        "TS" (two-sample analysis using mappability and GC content), and 
        "IO" (two-sample analysis without using mappability and GC content).
        If \code{analysisType="automatic"},
                this method tries to make the best guess for \code{analysisType}, 
                based on the data provided. }
    \item{bgEst}{Parameter to determine background estimation approach.
        Possible values are "matchLow" (estimation using bins with low tag counts) and
        "rMOM" (estimation using robust method of moment (MOM)).
        If \code{bgEst="automatic"},
                this method tries to make the best guess for \code{bgEst}, 
                based on the data provided.
		Default is \code{bgEst="rMOM"}. }    
    \item{k}{Parameter for estimating background distribution.
        It is not recommended for users to change this value. }    
    \item{meanThres}{Parameter for estimating background distribution.
        Default is 1 for \code{analysisType="TS"} and 0 for \code{analysisType="OS"}.
        Not relevant when \code{analysisType="IO"}. }     
    \item{s}{Parameter for estimating background distribution.
        Relevant only when \code{analysisType="TS"}. Default is 2. }  
    \item{d}{Parameter for estimating background distribution.
        Relevant only when \code{analysisType="TS"} or \code{analysisType="IO"}.
        Default is 0.25. }
    \item{trans}{Transformation of matching control tag count.
        Possible values are "log" (logarithm transformation) and
        "power" (power transformation).
		Relevant only when \code{analysisType="IO"}. 
		Default is \code{trans="power"}. }  
    \item{truncProb}{Parameter for estimating background distribution.
    	Relevant only when \code{analysisType="IO"}. }  
    \item{parallel}{Utilize multiple CPUs for parallel computing 
        using \code{"parallel"} package?
        Possible values are \code{TRUE} (utilize multiple CPUs) 
        or \code{FALSE} (do not utilize multiple CPUs).
        Default is \code{FALSE} (do not utilize multiple CPUs). }
    \item{nCore}{Number of CPUs when parallel computing is utilized. }
    \item{...}{ Other parameters to be passed through to generic \code{mosaicsFit}.}  
}
\details{
The imported data type constraints the analysis that can be implemented.
If only data for ChIP sample and matched control sample
(i.e., either \code{type=c("chip", "input")} or \code{type=c("chip", "input", "N")}
was used in method \code{readBins}),
only two-sample analysis without using mappability and GC content
(\code{analysisType="IO"}) is allowed.
If matched control data is available 
with mappability score, GC content score, and sequence ambiguity score,
(i.e., \code{type=c("chip", "input", "M", "GC", "N")} was used in method \code{readBins}),
user can do all of three analysis types
(\code{analysisType="OS"}, \code{analysisType="TS"}, or \code{analysisType="IO"}).
If there is no data for matched control sample
(i.e., \code{type=c("chip", "M", "GC", "N")} was used in method \code{readBins}),
only one-sample analysis (\code{analysisType="OS"}) is permitted.

Parallel computing can be utilized for faster computing
if \code{parallel=TRUE} and \code{parallel} package is loaded.
\code{nCore} determines number of CPUs used for parallel computing.
\code{meanThres}, \code{s}, \code{d}, \code{trans}, and \code{truncProb} are 
the tuning parameters for estimating background distribution.
The vignette and Kuan et al. (2011) provide further details about these tuning parameters.
Please do not try different value for \code{k} argument.
}
\value{
Construct \code{MosaicsFit} class object.
}
\references{
Kuan, PF, D Chung, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.
}
\author{ Dongjun Chung, Pei Fen Kuan, Sunduz Keles }
\seealso{
\code{\link{readBins}}, \code{\linkS4class{MosaicsFit}}.
}
\examples{
\dontrun{
library(mosaicsExample)
data(exampleBinData)

exampleFit <- mosaicsFit( exampleBinData, analysisType="IO", bgEst="automatic" )
}
}
\keyword{models}
\keyword{methods}
