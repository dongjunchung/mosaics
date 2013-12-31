\name{export}
\alias{export}
\alias{export,MosaicsPeak-method}
\title{
Export peak calling results to text files
}
\description{
Export peak calling results to text files in TXT, BED, or GFF file formats.
}
\usage{
export(object, ...)
\S4method{export}{MosaicsPeak}( object, type=NA, filename=NA )
}
\arguments{
    \item{object}{Object of class \code{MosaicsPeak},
        peak calling results obtained using method \code{mosaicsPeak}. }
    \item{type}{Format of the exported file.
        Possible values are \code{"txt"}, \code{"bed"}, and \code{"gff"}. See Details. }
    \item{filename}{Name of the exported file. }    
    \item{...}{ Other parameters to be passed through to generic \code{export}.}
}
\details{
TXT file format (\code{type="txt"}) exports peak calling results in the most informative way.
Columns include chromosome ID, peak start position, peak end position, peak width, 
average posterior probability, minimum posterior probability, 
average ChIP tag count, maximum ChIP tag count (always), 
average input tag count, average input tag count scaled by sequencing depth, 
average log base 2 ratio of ChIP over input tag counts 
(if matched control sample is also provided), 
average mappability score, and average GC content score
(when mappability and GC content scores are used in the analysis) in each peak. 
\code{type="bed"} and \code{type="gff"} export peak calling results in standard BED and GFF file formats,
respectively, where score is the average ChIP tag counts in each peak.
If no peak is detected, \code{export} method will not generate any file.
}
\references{
Kuan, PF, D Chung, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.
}
\author{ Dongjun Chung, Pei Fen Kuan, Sunduz Keles }
\seealso{
\code{\link{mosaicsPeak}}, \code{\linkS4class{MosaicsPeak}}.
}
\examples{
\dontrun{
library(mosaicsExample)
data(exampleFit)

examplePeak <- mosaicsPeak( exampleFit, signalModel = "2S", FDR = 0.05 )
export( examplePeak, type = "txt", filename = "./TSpeakList.txt" )
export( examplePeak, type = "bed", filename = "./TSpeakList.bed" )
export( examplePeak, type = "gff", filename = "./TSpeakList.gff" )
}
}
\keyword{models}
\keyword{methods}
