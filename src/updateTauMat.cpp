
#include <Rcpp.h>

RcppExport SEXP cpp_updateTauMat( SEXP tau, SEXP g, SEXP pi0 ) { 

	using namespace Rcpp;

    Rcpp::List tauList(tau);
    Rcpp::NumericMatrix gMat(g);
    Rcpp::NumericVector pi0Vec(pi0);
    
    int nState = gMat.nrow();
    int nObs = gMat.ncol();
    
    Rcpp::NumericMatrix tauMat( nState, nObs );
    double denom;
    
    // first element
    
    denom = 0;
    for (int g = 0; g < nState; g++) {         
        tauMat(g,0) = pi0Vec[g] * gMat(g,0);
        denom += tauMat(g,0);
    }
        
    for (int g = 0; g < nState; g++) {         
        tauMat(g,0) /= denom;
    }
    
    // main iteration
    
    for (int i = 1; i < nObs; i++ ) {
        Rcpp::NumericMatrix tauCurrent = tauList[(i-1)];
        
        for (int g = 0; g < nState; g++) {            
            tauMat(g,i) = 0;
            for (int g1 = 0; g1 < nState; g1++) {
                tauMat(g,i) += tauCurrent(g1,g);
            }
        }
    }
    
    return tauMat;
}
