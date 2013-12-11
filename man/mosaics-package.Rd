\name{mosaics-package}
\alias{mosaics-package}
\alias{mosaics}
\docType{package}
\title{
MOSAiCS (MOdel-based one and two Sample Analysis and Inference for ChIP-Seq)
}
\description{
This package provides functions for fitting MOSAiCS,
a statistical framework to analyze one-sample or two-sample ChIP-seq data.
}
\details{
\tabular{ll}{
Package: \tab mosaics\cr
Type: \tab Package\cr
Version: \tab 1.6.0\cr
Date: \tab 2013-07-30\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
This package contains three main classes, \code{BinData}, \code{MosaicsFit}, and \code{MosaicsPeak},
which represent bin-level ChIP-seq data, MOSAiCS model fit, and MOSAiCS peak calling results, respectively.
This package contains four main methods, 
\code{constructBins}, \code{readBins}, \code{mosaicsFit}, and \code{mosaicsPeak}.
\code{constructBins} method constructs bin-level files from the aligned read file.
\code{readBins} method imports bin-level data and construct \code{BinData} class object.
\code{mosaicsFit} method fits a MOSAiCS model using \code{BinData} class object and
 constructs \code{MosaicsFit} class object.
\code{mosaicsPeak} method calls peaks using \code{MosaicsFit} class object and
 construct \code{MosaicsPeak} class object.
\code{MosaicsPeak} class object can be exported as text files or transformed into data frame,
which can be used for the downstream analysis.
This package also provides methods for simple exploratory analysis.

The \code{mosaics} package companion website, \url{http://www.stat.wisc.edu/~keles/Software/mosaics/},
provides preprocessing scripts, preprocessed files for diverse reference genomes,
and easy-to-follow instructions.
We encourage questions or requests regarding \code{mosaics} package to be posted on
our Google group, \url{http://groups.google.com/group/mosaics_user_group}.
Please check the vignette for further details on the \code{mosaics} package and these websites.
}
\author{
Dongjun Chung, Pei Fen Kuan, Sunduz Keles

Maintainer: Dongjun Chung <chungdon@stat.wisc.edu>
}
\references{
Kuan, PF, D Chung, G Pan, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.
}
\keyword{ package }
\seealso{
\code{\link{constructBins}}, \code{\link{readBins}}, \code{\link{mosaicsFit}}, \code{\link{mosaicsPeak}} ,
\code{\linkS4class{BinData}}, \code{\linkS4class{MosaicsFit}}, \code{\linkS4class{MosaicsPeak}}.
}
\examples{
\dontrun{
library(mosaicsExample)

exampleBinData <- readBins( type=c("chip","input"),
    fileName=c( system.file("extdata/chip_chr21.txt", package="mosaicsExample"),
    system.file("extdata/input_chr21.txt", package="mosaicsExample") ) )
exampleBinData
print(exampleBinData)[1:10, ]
plot(exampleBinData)
plot( exampleBinData, plotType="input" )

exampleFit <- mosaicsFit( exampleBinData, analysisType="IO" )
exampleFit
plot(exampleFit)
estimates(exampleFit)

examplePeak <- mosaicsPeak( exampleFit, signalModel = "2S", FDR = 0.05 )
examplePeak
print(examplePeak)[1:10, ]
export( examplePeak, type = "txt", filename = "./TSpeakList.txt" )
export( examplePeak, type = "bed", filename = "./TSpeakList.bed" )
export( examplePeak, type = "gff", filename = "./TSpeakList.gff" )
}
}
