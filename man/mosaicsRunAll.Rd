\name{mosaicsRunAll}
\alias{mosaicsRunAll}
\title{
Analyze ChIP-seq data using the MOSAiCS framework
}
\description{
Construct bin-level ChIP-sep data from aligned read files of ChIP and matched control samples,
fit a MOSAiCS model, call peaks, export peak calling results, 
and generate reports for diagnostics.
}
\usage{
mosaicsRunAll( 
    chipFile=NULL, chipFileFormat=NULL, 
    controlFile=NULL, controlFileFormat=NULL, 
    binfileDir=NULL, 
    peakFile=NULL, peakFileFormat=NULL,
    reportSummary=FALSE, summaryFile=NULL, 
    reportExploratory=FALSE, exploratoryFile=NULL, 
    reportGOF=FALSE, gofFile=NULL, 
    PET=FALSE, byChr=FALSE, useChrfile=FALSE, chrfile=NULL, excludeChr=NULL, 
    FDR=0.05, fragLen=200, binSize=200, capping=0, bgEst="rMOM", d=0.25, 
    signalModel="BIC", maxgap=200, minsize=50, 
    thres=10, parallel=FALSE, nCore=8 )
}
\arguments{
  \item{chipFile}{
    Name of the aligned read file of ChIP sample to be processed.
}
  \item{chipFileFormat}{
    Format of the aligned read file of ChIP sample to be processed.
    Currently, \code{mosaicsRunAll} permits the following aligned read file formats:
        \code{"eland_result"} (Eland result), \code{"eland_extended"} (Eland extended),
        \code{"eland_export"} (Eland export), \code{"bowtie"} (default Bowtie),
        \code{"sam"} (SAM), \code{"bed"} (BED), and \code{"csem"} (CSEM BED).
        Note that \code{"csem"} does not mean CSEM output file format, but CSEM BED file format.
}
  \item{controlFile}{
    Name of the aligned read file of matched control sample to be processed.
}
  \item{controlFileFormat}{
    Format of the aligned read file of matched control sample to be processed.
    Currently, \code{mosaicsRunAll} permits the following aligned read file formats:
        \code{"eland_result"} (Eland result), \code{"eland_extended"} (Eland extended),
        \code{"eland_export"} (Eland export), \code{"bowtie"} (default Bowtie),
        \code{"sam"} (SAM), \code{"bed"} (BED), and \code{"csem"} (CSEM BED).
        Note that \code{"csem"} does not mean CSEM output file format, but CSEM BED file format.
}
  \item{binfileDir}{
    Directory to store processed bin-level files.
}
  \item{peakFile}{
    Name of the peak list generated from the analysis.
}
  \item{peakFileFormat}{
    Format of the peak list generated from the analysis.
    Possible values are \code{"txt"}, \code{"bed"}, and \code{"gff"}.
}
  \item{reportSummary}{
    Report the summary of model fitting and peak calling? 
    Possible values are \code{TRUE} (YES) and \code{FALSE} (NO).
    Default is \code{FALSE} (NO). 
}
  \item{summaryFile}{
    File name of the summary report of model fitting and peak calling.
    The summary report is a text file.
}
  \item{reportExploratory}{
    Report the exploratory analysis plots? 
    Possible values are \code{TRUE} (YES) and \code{FALSE} (NO).
    Default is \code{FALSE} (NO). 
}
  \item{exploratoryFile}{
    Name of the file for exploratory analysis plots.
    The exploratory analysis results are exported as a PDF file.
}
  \item{reportGOF}{
    Report the goodness of fit (GOF) plots? 
    Possible values are \code{TRUE} (YES) and \code{FALSE} (NO).
    Default is \code{FALSE} (NO). 
}
  \item{gofFile}{
    Name of the file for goodness of fit (GOF) plots.
    The GOF plots are exported as a PDF file.
}
  \item{PET}{
    Is the file paired-end tag (PET) data?
    If \code{PET=FALSE}, it is assumed that the file is SET data. 
    If \code{PET=TRUE}, it is assumed that the file is PET data.    
    Default is \code{FALSE} (SET data).
}
  \item{byChr}{
    Analyze ChIP-seq data for each chromosome separately or analyze it genome-wide?
    Possible values are \code{TRUE} (chromosome-wise) and \code{FALSE} (genome-wide).
    Default is \code{FALSE} (genome-wide analysis).
}
  \item{useChrfile}{
    Is the file for chromosome info provided?
    Possible values are \code{TRUE} or \code{FALSE}.
    If \code{useChrfile=FALSE}, it is assumed that the file for chromosome info is not provided. 
    If \code{useChrfile=TRUE}, it is assumed that the file for chromosome info is provided.
    Default is \code{FALSE}.
}
  \item{chrfile}{
    Name of the file for chromosome info.
    In this file, the first and second columns are ID and size of each chromosome, respectively.
}
  \item{excludeChr}{
    Vector of chromosomes that will be excluded from the analysis.
}  
  \item{FDR}{
    False discovery rate. Default is 0.05.
}  
  \item{fragLen}{
    Average fragment length. Default is 200.
}
  \item{binSize}{
    Size of bins. Default is 200.
}
  \item{capping}{
    Maximum number of reads allowed to start at each nucleotide position. 
    To avoid potential PCR amplification artifacts, the maximum number of reads
    that can start at a nucleotide position is capped at \code{capping}. 
    Capping is not applied if non-positive \code{capping} is used.
    Default is 0 (no capping).
}  
  \item{bgEst}{Parameter to determine background estimation approach.
    Possible values are "matchLow" (estimation using bins with low tag counts) and
    "rMOM" (estimation using robust method of moment (MOM)).
    If \code{bgEst="automatic"},
            this method tries to make the best guess for \code{bgEst}, 
            based on the data provided.
        Default is \code{bgEst="rMOM"}. }    
  \item{d}{Parameter for estimating background distribution.
    Default is 0.25. } 
  \item{signalModel}{Signal model.    
    Possible values are "BIC" (automatic model selection using BIC), 
    "1S" (one-signal-component model), and
    "2S" (two-signal-component model). Default is "BIC". } 
  \item{maxgap}{Initial nearby peaks are merged if the distance (in bp) 
    between them is less than \code{maxgap}. Default is 200. }     
  \item{minsize}{An initial peak is removed if its width is narrower than \code{minsize}. 
    Default is 50. }
  \item{thres}{A bin within initial peak is removed if its ChIP tag counts are less than \code{thres}.
    Default is 10. }
  \item{parallel}{Utilize multiple CPUs for parallel computing 
    using \code{"parallel"} package?
    Possible values are \code{TRUE} (use multiple CPUs) 
    or \code{FALSE} (do not use multiple CPUs).
    Default is \code{FALSE} (do not use multiple CPUs).}
  \item{nCore}{Number of maximum number of CPUs used for the analysis. 
    Default is 8. }
}
\details{
This method implements the work flow for the two-sample analysis of ChIP-seq data 
using the MOSAiCS framework (without using mappability and GC content scores).
It imports aligned read files of ChIP and matched control samples,
processes them into bin-level files, fits MOSAiCS model, calls peaks, 
exports the peak lists to text files, and generates reports for diagnostics. 
This method is a wrapper function of \code{constructBins}, \code{readBins},
\code{mosaicsFit}, \code{mosaicsPeak}, \code{export} functions, and methods of
\code{BinData}, \code{MosaicsFit}, and \code{MosaicsPeak} classes.

See the vignette of the package for the illustration of the work flow
and the description of employed methods and their options.
Exploratory analysis plots and goodness of fit (GOF) plots are generated 
using the methods \code{plot} of the classes \code{BinData} and \code{MosaicsFit}, respectively.
See the help of \code{constructBins} for details of the options \code{PET}, \code{chipFileFormat},
\code{controlFileFormat}, \code{byChr}, \code{useChrfile}, \code{chrfile}, \code{excludeChr}, 
\code{fragLen}, \code{binSize}, and \code{capping}.
See the help of \code{mosaicsFit} 
for details of the options \code{bgEst} and \code{d}.
See the help of \code{mosaicsPeak} for details of the options \code{FDR}, 
\code{signalModel}, \code{maxgap}, \code{minsize}, and \code{thres}.
See the help of \code{export} for details of the option \code{peakFileFormat}.

When the data contains multiple chromosomes, 
parallel computing can be utilized for faster preprocessing and model fitting
if \code{parallel=TRUE} and \code{parallel} package is loaded.
\code{nCore} determines number of CPUs used for parallel computing.
}
\value{
Processed bin-level files are exported to the directory specified in \code{binfileDir} argument.
If \code{byChr=FALSE} (genome-wide analysis), 
one bin-level file is generated for each of ChIP and matched control samples,
where file names are \code{[chipFile]_fragL[fragLen]_bin[binSize].txt}
and \code{[controlFile]_fragL[fragLen]_bin[binSize].txt}, respectively, for SET data (\code{PET = FALSE}).
For PET data (\code{PET = TRUE}), file names for each of ChIP and matched control samples
are \code{[chipFile]_bin[binSize].txt} and \code{[controlFile]_bin[binSize].txt}, respectively.
If \code{byChr=TRUE} (chromosome-wise analysis),
bin-level files are generated for each chromosome of each of ChIP and matched control samples,
where file names are \code{[chipFile]_fragL[fragLen]_bin[binSize]_[chrID].txt}
and \code{[controlFile]_fragL[fragLen]_bin[binSize]_[chrID].txt}, respectively, for SET data (\code{PET = FALSE})
(\code{[chrID]} is chromosome IDs that reads align to).
For PET data (\code{PET = TRUE}), file names for each of ChIP and matched control samples
are \code{[chipFile]_bin[binSize]_[chrID].txt} 
and \code{[controlFile]_bin[binSize]_[chrID].txt}, respectively.

The peak list generated from the analysis are exported 
to the file with the name specified in \code{peakFile}.
If \code{reportSummary=TRUE}, the summary of model fitting and peak calling is exported 
to the file with the name specified in \code{summaryFile} (text file).
If \code{reportExploratory=TRUE}, the exploratory analysis plots are exported 
to the file with the name specified in \code{exploratoryFile} (PDF file).
If \code{reportGOF=TRUE}, the goodness of fit (GOF) plots are exported 
to the file with the name specified in \code{gofFile} (PDF file).
}
\references{
Kuan, PF, D Chung, G Pan, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.
}
\author{ Dongjun Chung, Pei Fen Kuan, Sunduz Keles }
\seealso{
\code{\link{constructBins}}, \code{\link{readBins}}, 
\code{\link{mosaicsFit}}, \code{\link{mosaicsPeak}}, \code{\link{export}},
\code{\linkS4class{BinData}}, \code{\linkS4class{MosaicsFit}}, \code{\linkS4class{MosaicsPeak}}.
}
\examples{
\dontrun{
# minimal input (without any reports for diagnostics)

mosaicsRunAll( 
    chipFile = "/scratch/eland/STAT1_eland_results.txt", 
    chipFileFormat = "eland_result", 
    controlFile = "/scratch/eland/input_eland_results.txt", 
    controlFileFormat = "eland_result", 
    binfileDir = "/scratch/bin/", 
    peakFile = "/scratch/peak/STAT1_peak_list.bed", 
    peakFileFormat = "bed" )
    
# generate all reports for diagnostics  

library(parallel)    
mosaicsRunAll( 
    chipFile = "/scratch/eland/STAT1_eland_results.txt", 
    chipFileFormat = "eland_result", 
    controlFile = "/scratch/eland/input_eland_results.txt", 
    controlFileFormat = "eland_result", 
    binfileDir = "/scratch/bin/", 
    peakFile = "/scratch/peak/STAT1_peak_list.bed", 
    peakFileFormat = "bed",
    reportSummary = TRUE, 
    summaryFile = "/scratch/reports/mosaics_summary.txt", 
    reportExploratory = TRUE, 
    exploratoryFile = "/scratch/reports/mosaics_exploratory.pdf", 
    reportGOF = TRUE, 
    gofFile = "/scratch/reports/mosaics_GOF.pdf",
    PET = FALSE, byChr = FALSE, useChrfile=FALSE, chrfile=NULL, excludeChr = "chrM", 
    FDR = 0.05,  fragLen = 200, capping = 0, bgEst="automatic", thres=10,
    parallel = TRUE, nCore = 8 )
}
}
\keyword{models}
\keyword{methods}
