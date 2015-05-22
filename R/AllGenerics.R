
# generic methods for "BinData" class

setGeneric( "chrID",
    function( object, ... )
    standardGeneric("chrID")
)

setGeneric( "coord",
    function( object, ... )
    standardGeneric("coord")
)

setGeneric( "tagCount",
    function( object, ... )
    standardGeneric("tagCount")
)

setGeneric( "mappability",
    function( object, ... )
    standardGeneric("mappability")
)

setGeneric( "gcContent",
    function( object, ... )
    standardGeneric("gcContent")
)

setGeneric( "input",
    function( object, ... )
    standardGeneric("input")
)

setGeneric( "mosaicsFit",
    function( object, ... )
    standardGeneric("mosaicsFit")
)

# generic methods for "MosaicsFit" class

setGeneric( "estimates",
    function( object, ... )
    standardGeneric("estimates")
)

setGeneric( "mosaicsPeak",
    function( object, ... )
    standardGeneric("mosaicsPeak")
)

setGeneric( "mosaicsFitHMM",
    function( object, ... )
    standardGeneric("mosaicsFitHMM")
)

# generic methods for "MosaicsHMM" class

setGeneric( "mosaicsPeakHMM",
    function( object, ... )
    standardGeneric("mosaicsPeakHMM")
)

# generic methods for "MosaicsPeak" class

setGeneric( "extractReads",
    function( object, ... )
    standardGeneric("extractReads")
)

setGeneric( "findSummit",
    function( object, ... )
    standardGeneric("findSummit")
)

setGeneric( "adjustBoundary",
    function( object, ... )
    standardGeneric("adjustBoundary")
)

setGeneric( "filterPeak",
    function( object, ... )
    standardGeneric("filterPeak")
)

#setGeneric( "peakList",
#    function( object, ... )
#    standardGeneric("peakList")
#)

setGeneric( "export",
    function( object, ... )
    standardGeneric("export")
)

setGeneric( "bdBin",
    function( object, ... )
    standardGeneric("bdBin")
)

setGeneric( "empFDR",
    function( object, ... )
    standardGeneric("empFDR")
)

setGeneric( "readCoverage",
    function( object, ... )
    standardGeneric("readCoverage")
)

setGeneric( "read",
    function( object, ... )
    standardGeneric("read")
)

setGeneric( "seqDepth",
    function( object, ... )
    standardGeneric("seqDepth")
)

setGeneric( "postProb",
    function( object, ... )
    standardGeneric("postProb")
)
