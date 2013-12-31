
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

setGeneric( "peakList",
    function( object, ... )
    standardGeneric("peakList")
)

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
