
#include <Rcpp.h>

RcppExport SEXP cpp_updatePiMat( SEXP tauL, SEXP tauM ) { 

	using namespace Rcpp;

    Rcpp::List tauList(tauL);
    Rcpp::NumericMatrix tauMat(tauM);
    
    int nState = tauMat.nrow();
    int nObs = tauMat.ncol();
    
    Rcpp::NumericMatrix piMat( nState, nState );
    
    // calculate denominator
    
    Rcpp::NumericVector denom( nState );
    for (int g = 0; g < nState; g++) {   
        denom[g] = 0;
        for (int i = 0; i<(nObs-1); i++ ) {
            denom[g] += tauMat(g,i);
        }
    }
    
    // calculate pi matrix
    
    for (int g1 = 0; g1 < nState; g1++) {   
        for (int g2 = 0; g2 < nState; g2++) {   
            piMat(g1,g2) = 0;
            for (int i = 0; i<(nObs-1); i++ ) {
                Rcpp::NumericMatrix tauCurrent = tauList[i];
                piMat(g1,g2) += tauCurrent(g1,g2);
            }
            piMat(g1,g2) /= denom[g1];
        }
    }
    
    return piMat;
}
