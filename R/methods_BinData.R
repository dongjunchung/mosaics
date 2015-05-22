
# generic methods for "BinData" class

setMethod(
    f="show",
    signature="BinData",
    definition=function( object ) {            
        # genome-wide summary
        
        chrList <- sort(unique( chrID(object) ))
        
        cat( "Summary: bin-level data (class: BinData)\n" )
        cat( "----------------------------------------\n" )
        cat( "- # of chromosomes in the data: ", length(chrList), "\n", sep="" )
        cat( "- total effective tag counts: ", sum(tagCount(object)), "\n", sep="" )
        cat( "  (sum of ChIP tag counts of all bins)\n" )
        if ( length(object@input)>0 )
        {
            cat( "- control sample is incorporated\n" )
        } else
        {
            cat( "- control sample is NOT incorporated\n" )        
        }
        if ( length(object@mappability)>0 )
        {
            cat( "- mappability score is incorporated\n" )
        } else
        {
            cat( "- mappability score is NOT incorporated\n" )        
        }
        if ( length(object@gcContent)>0 )
        {
            cat( "- GC content score is incorporated\n" )
        } else
        {
            cat( "- GC content score is NOT incorporated\n" )        
        }
        switch( object@dataType,
            unique = {
                cat( "- uni-reads are assumed\n" )
            },
            multi = {            
                cat( "- both uni-reads and multi-reads are assumed\n" )
            }
        )
        cat( "----------------------------------------\n" )  
        
        # chromosome-wise summary
        
        if ( length(chrList) > 1 ) {
            coordByChr <- split( coord(object), chrID(object) )
            countByChr <- split( tagCount(object), chrID(object) )
            
            cat( "chrID:\t# of bins, total effective tag counts\n" )
            for ( chr in 1:length(chrList) ) {
                nCord <- length( coordByChr[[ chrList[chr] ]] )
                nCount <- sum( countByChr[[ chrList[chr] ]] )
                cat( chrList[chr],":\t",nCord,", ",nCount,"\n", sep="" )
            }
            cat( "----------------------------------------\n" )  
            
            rm( coordByChr, countByChr )
            gc()
        }
    }
)

setMethod(
    f="print",
    signature="BinData",
    definition=function( x ) {
        if ( length(x@input)>0 )
        {
            if ( length(x@mappability)>0 & length(x@gcContent)>0 )
            {
                # two-sample analysis (with M & GC)
                
                printForm <- data.frame(
                    chrID=chrID(x), coord=coord(x), tagCount=tagCount(x),
                    mappability=mappability(x), gcContent=gcContent(x),
                    input=input(x) )   
            } else
            {
                # two-sample analysis (Input only)
                
                printForm <- data.frame(
                    chrID=chrID(x), coord=coord(x), tagCount=tagCount(x), 
                    input=input(x) )             
            }
        } else
        {   
            # one-sample analysis
            
            printForm <- data.frame(
                chrID=chrID(x), coord=coord(x), tagCount=tagCount(x),
                mappability=mappability(x), gcContent=gcContent(x) )  
        }
        return(printForm)
    }
)

setMethod(
    f="plot",
    signature=c("BinData","missing"),
    definition=function( x, y, plotType=NULL, inputGrid=c(-1:10,15,20,30,50,100), ... ) {
        if ( !is.null(plotType) ) {
            if ( plotType=="M" ) {
                # plot mean tag count vs. mappability, with 95% CI
                         
                if ( length(x@mappability)>0 ) {            
                    statM <- .computeStat( Y=tagCount(x), S=mappability(x) ) 
                    plot( statM$uS, statM$meanYall,
                        xlab='Mappability score', ylab='Mean ChIP tag count', 
                        main='Mappability score vs. Mean ChIP tag count',
                        ylim=quantile( statM$meanYall, prob=c(0.05,0.95) ) )
                    segments( statM$uS, statM$meanYall, 
                        statM$uS, statM$meanYall+1.96*sqrt(statM$varYall/statM$nitem) )
                    segments( statM$uS, statM$meanYall, 
                        statM$uS, statM$meanYall-1.96*sqrt(statM$varYall/statM$nitem) )
                } else {
                    stop( "bin-level data does not include mappability score information!" )
                }
            } else if ( plotType=="GC" ) {
                # plot mean tag count vs. GC content, with 95% CI
                                        
                if ( length(x@gcContent)>0 ) {
                    statGC <- .computeStat( Y=tagCount(x), S=gcContent(x) )
                    plot( statGC$uS, statGC$meanYall,
                        xlab='GC content score', ylab='Mean ChIP tag count', 
                        main='GC content score vs. Mean ChIP tag count',
                        ylim=quantile( statGC$meanYall, prob=c(0.05,0.95) ) )
                    segments( statGC$uS, statGC$meanYall, 
                        statGC$uS, statGC$meanYall+1.96*sqrt(statGC$varYall/statGC$nitem) )
                    segments( statGC$uS, statGC$meanYall, 
                        statGC$uS, statGC$meanYall-1.96*sqrt(statGC$varYall/statGC$nitem) ) 
                } else {
                    stop( "bin-level data does not include GC content score information!" )
                }        
            } else if ( plotType=="input" ) {
                # plot mean tag count vs. input mean tag count, with 95% CI
                                        
                if ( length(x@input)>0 ) {
                    statInput <- .computeStat( Y=tagCount(x), S=input(x) )
                    plot( statInput$uS, statInput$meanYall,
                        xlab='Control tag count', ylab='Mean ChIP tag count', 
                        main='Control tag count vs. Mean ChIP tag count',
                        xlim=c(0,quantile(input(x),0.9999)) )
                    segments( statInput$uS, statInput$meanYall, 
                        statInput$uS, 
                        statInput$meanYall+1.96*sqrt(statInput$varYall/statInput$nitem) )
                    segments( statInput$uS, statInput$meanYall, 
                        statInput$uS, 
                        statInput$meanYall-1.96*sqrt(statInput$varYall/statInput$nitem) ) 
                } else {
                    stop( "bin-level data does not include control tag count information!" )
                }        
            } else if ( plotType=="M|input" ) {
                # plot mean tag count vs. mappability, conditional on input
                                        
                if ( length(x@input)>0 ) {
                    #library(lattice)
                      
                    M_all <- X_all <- Y_all <- input_mat <- c()
                        
                    for( i in 1:(length(inputGrid)-1) ) {
                        sub_id <- which( input(x)>inputGrid[i] & input(x)<=inputGrid[(i+1)] )
                        if ( length(sub_id) > 0 ) {
                            statM <- .computeStat( Y=tagCount(x)[sub_id], S=mappability(x)[sub_id] )
                            meanYall <- statM$meanYall
                                
                            MaxY <- 1 
                            if( length(meanYall) > 0 ) {
                                MaxY <- max(meanYall)
                            }
                            if( MaxY==0 ) {
                                MaxY <- 1
                            }
                        
                            input_mat <- rbind( input_mat, inputGrid[i:(i+1)] )
                            M_all <- c( M_all, statM$uS )
                            Y_all <- c( Y_all, meanYall/MaxY )
                            X_all <- c( X_all, rep(mean(inputGrid[i:(i+1)]),length(statM$uS)) )
                        }
                    }
                    Input <- shingle( X_all, intervals=input_mat )
                    print( xyplot( Y_all ~ M_all | Input,
                        xlab='Mappability score', ylab='Mean ChIP tag count',
                        main='Mappability score vs. Mean ChIP tag count,\nconditional on Control tag count', cex=0.5 ) ) 
                } else {
                    stop( "bin-level data does not include control tag count information!" )
                }                   
            } else if ( plotType=="GC|input" ) {
                # plot mean tag count vs. GC content, conditional on input
                                        
                if ( length(x@input)>0 ) {
                    #library(lattice)
                          
                    GC_all <- X_all <- Y_all <- input_mat <- c()
                    
                    for( i in 1:(length(inputGrid)-1) ) {
                        sub_id <- which( input(x)>inputGrid[i] & input(x)<=inputGrid[(i+1)] )
                        if ( length(sub_id) > 0 ) {
                            statGC <- .computeStat( Y=tagCount(x)[sub_id], S=gcContent(x)[sub_id] )
                            meanYall <- statGC$meanYall
                        
                            MaxY <- 1 
                            if( length(meanYall) > 0 ) {
                                MaxY <- max(meanYall)
                            }
                            if( MaxY==0 ) {
                                MaxY <- 1
                            }
                        
                            input_mat <- rbind( input_mat, inputGrid[i:(i+1)] )
                            GC_all <- c( GC_all, statGC$uS )
                            Y_all <- c( Y_all, meanYall/MaxY )
                            X_all <- c( X_all, rep(mean(inputGrid[i:(i+1)]),length(statGC$uS)) )
                        }
                    }
                    Input <- shingle( X_all, intervals=input_mat )
                    print( xyplot( Y_all ~ GC_all | Input,
                        xlab='GC content score', ylab='Mean ChIP tag count',
                        main='GC content score vs. Mean ChIP tag count,\nconditional on Control tag count', cex=0.5 ) ) 
                } else {
                    stop( "bin-level data does not include control tag count information!" )
                }                   
            } 
        } else {
            # if "plotType" is missing, just show histogram of tag counts (& input, if possible)
                 
            YFreq <- table( tagCount(x) )
            YVal <- as.numeric( names(YFreq) )   
            
            if ( length(x@input)>0 ) {                       
                XFreq <- table( input(x) )
                XVal <- as.numeric( names(XFreq) )   
            }
            
            plot( log10(YVal+1), log10(YFreq), type='l', axes=FALSE, 
                ylab='Frequency', xlab='Tag count', main='Histogram of tag count' )
            
        if ( length(x@input)>0 ) {                       
           points( log10(XVal+1), log10(XFreq), type='l', col='darkgray' )
        }
            axis(1,0:6,10^c(0:6)-1)
            axis(2,0:10,10^c(0:10))
            
            if ( length(x@input)>0 ) {
                legend( 1, log10(max(YFreq)+1), c('ChIP','control'),
                    col=c('black','darkgray'),lty=c(1,1),bty='n')
            } else {
                legend( 1, log10(max(YFreq)+1), c('ChIP'),
                    col=c('black'),lty=1,bty='n')
            }
        }
    }
)

setMethod(
    f="chrID",
    signature="BinData",
    definition=function( object ) {
        return(object@chrID)
    }
)

setMethod(
    f="coord",
    signature="BinData",
    definition=function( object ) {
        return(object@coord)
    }
)

setMethod(
    f="tagCount",
    signature="BinData",
    definition=function( object ) {
        return(object@tagCount)
    }
)

setMethod(
    f="mappability",
    signature="BinData",
    definition=function( object ) {
        if ( length(object@mappability)>0 )
        {
            return(object@mappability)
        } else
        {
            warning( "no mappability score information provided!" )
            return(NULL)
        }   
    }
)

setMethod(
    f="gcContent",
    signature="BinData",
    definition=function( object ) {
        if ( length(object@gcContent)>0 )
        {
            return(object@gcContent)
        } else
        {
            warning( "no GC content score information provided!" )
            return(NULL)
        }   
    }
)

setMethod(
    f="input",
    signature="BinData",
    definition=function( object ) {
        if ( length(object@input)>0 )
        {
            return(object@input)
        } else
        {
            warning( "no control sample information provided!" )
            return(NULL)
        }   
    }
)

setMethod(
    f="seqDepth",
    signature="BinData",
    definition=function( object ) {
      return(object@seqDepth)
    }
)
