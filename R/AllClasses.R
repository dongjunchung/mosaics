
# bin-level data

setClass( Class="BinData",
    representation=representation(
        chrID="character",
        coord="numeric",
        tagCount="numeric",
        mappability="numeric",
        gcContent="numeric",
        input="numeric",
        dataType="character"
    )
)

# MOSAiCS model fit

setClass( Class="MosaicsFitEst",
    representation=representation(
        pi0="numeric",
        a="numeric",
        betaEst="numeric",
        muEst="numeric",    
        pNfit="list",
        b="numeric",
        c="numeric",
        p1="numeric",
        b1="numeric",
        c1="numeric",
        b2="numeric",
        c2="numeric",
        inputTrunc="numeric",
        analysisType="character"
    )
)

setClass( Class="MosaicsFitParam",
    representation=representation(
        k="numeric",
        meanThres="numeric",
        s="numeric",
        d="numeric"
    )
)

setClass( Class="MosaicsFit",
    representation=representation(
        mosaicsEst="MosaicsFitEst",
        mosaicsParam="MosaicsFitParam",
        chrID="character",
        coord="numeric",
        tagCount="numeric",
        mappability="numeric",
        gcContent="numeric",
        input="numeric",
        bic1S="numeric",
        bic2S="numeric"
    )
)

# MOSAiCS peak calling & HMM

setClass( Class="MosaicsPeakParam",
    representation=representation(
        analysisType="character",
        signalModel="character",
        FDR="numeric",
        maxgap="numeric",
        minsize="numeric",
        thres="numeric",
        decoding="character"
    )
)

setClass( Class="MosaicsPeak",
    representation=representation(
        peakList="data.frame",
        peakParam="MosaicsPeakParam",
        bdBin="data.frame",
        postProb="data.frame",
        empFDR="numeric"
    )
)

setClass( Class="MosaicsHMM",
    representation=representation(
        HMMfit="list",
        mosaicsEst="MosaicsFitEst",
		inputdata="list",
        init="character",
        initPiMat="matrix",
        peakParam="MosaicsPeakParam",
        binsize="numeric",
        nRatio="numeric",
        bicMosaics="numeric",
        bicMosaicsHMM="numeric"
    )
)
