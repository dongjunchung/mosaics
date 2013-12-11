
#include <Rcpp.h>

RcppExport SEXP cpp_viterbi( SEXP pi, SEXP g, SEXP pi0 ) { 

	using namespace Rcpp;

    Rcpp::NumericMatrix piMat(pi), gMat(g);
    Rcpp::NumericVector pi0Vec(pi0);
    
    int nState = gMat.nrow();
    int nObs = gMat.ncol();
    
    Rcpp::NumericMatrix vMat( nState, nObs );
    Rcpp::IntegerMatrix ptr( nState, nObs );
    Rcpp::IntegerVector path( nObs );
    
    double maxVal, kVal;
    
    // calculate log(pi) for faster calculation
    
    Rcpp::NumericMatrix logPiMat( nState, nState );
    for ( int g1 = 0; g1 < nState; g1++ ) {
        for ( int g2 = 0; g2 < nState; g2++ ) {
            logPiMat(g1,g2) = log( piMat(g1,g2) );
        }
    }
    
    // initialization
    
    for (int g = 0; g < nState; g++) {
        vMat(g,0) = log( gMat(g,0) ) + log( pi0Vec[g] );
        ptr(g,0) = 0;
    }
    
    // iteration
    
    for (int i = 1; i < nObs; i++ ) {
        for (int g = 0; g < nState; g++) {
            maxVal = vMat(0,(i-1)) + logPiMat(0,g);
            for (int g1 = 1; g1 < nState; g1++) {
                kVal = vMat(g1,(i-1)) + logPiMat(g1,g);
                if ( kVal > maxVal ) {
                    maxVal = kVal;
                    ptr(g,i) = g1;
                }
                vMat(g,i) = log( gMat(g,i) ) + maxVal;
            }            
        }
    }
    
    // termination
    
    for (int g = 0; g < nState; g++) {
        maxVal = vMat(0,(nObs-1));
        for (int g1 = 1; g1 < nState; g1++) {
            kVal = vMat(g1,(nObs-1));
            if ( kVal > maxVal ) {
                maxVal = kVal;
                path[(nObs-1)] = g1;
            }
        }            
    }
    
    // traceback
    
    for (int i = nObs-2; i >= 0; i-- ) {
        path[i] = ptr( path[(i+1)], (i+1) );
    }
    
    return path;
}
