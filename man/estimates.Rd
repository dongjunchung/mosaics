\name{estimates}
\alias{estimates}
\alias{estimates,MosaicsFit-method}
\title{
Extract estimates of the fitted MOSAiCS model
}
\description{
Extract estimates from \code{MosaicsFit} class object, which is a fitted MOSAiCS model.
}
\usage{
estimates( object, ... )
\S4method{estimates}{MosaicsFit}( object )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{Object of class \code{MosaicsFit}, which represents
        a fitted MOSAiCS model obtained using method \code{mosaicsFit}. }
    \item{...}{ Other parameters to be passed through to generic \code{estimates}.}
}
\value{
Returns a list with components:
    \item{pi0}{ Mixing proportion of background component. }
    \item{a}{ Parameter for background component. }
    \item{betaEst}{ Parameter for background component (coefficient estimates). }
    \item{muEst}{ Parameter for background component. }
    \item{b}{ Parameter for one-signal-component model. }
    \item{c}{ Parameter for one-signal-component model. }
    \item{p1}{ Parameter for two-signal-component model (mixing proportion of signal components). }
    \item{b1}{ Parameter for two-signal-component model (the first signal component). }
    \item{c1}{ Parameter for two-signal-component model (the first signal component). }
    \item{b2}{ Parameter for two-signal-component model (the second signal component). }
    \item{c2}{ Parameter for two-signal-component model (the second signal component). }
    \item{analysisType}{ Analysis type.
        Possible values are "OS" (one-sample analysis),
        "TS" (two-sample analysis using mappability and GC content), and 
        "IO" (two-sample analysis without using mappability and GC content). }
}
\references{
Kuan, PF, D Chung, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.
}
\author{ Dongjun Chung, Pei Fen Kuan, Sunduz Keles }
\seealso{
\code{\link{mosaicsFit}}, \code{\linkS4class{MosaicsFit}}.
}
\examples{
\dontrun{
library(mosaicsExample)
data(exampleFit)

estimates(exampleFit)
}
}
\keyword{models}
\keyword{methods}
