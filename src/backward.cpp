
#include <Rcpp.h>

RcppExport SEXP cpp_backward( SEXP pi, SEXP g, SEXP s ) { 

	using namespace Rcpp;

    Rcpp::NumericMatrix piMat(pi), gMat(g);
    Rcpp::NumericVector sVec(s);
    
    int nState = gMat.nrow();
    int nObs = gMat.ncol();
    Rcpp::NumericMatrix bMat( nState, nObs );
    
    // initialization
    
    for (int g = 0; g < nState; g++) {
        bMat(g,(nObs-1)) = 1;
    }
    
    // iteration
    
    for (int i = nObs-2; i>=0; i-- ) {
        for (int g = 0; g < nState; g++) {
            bMat(g,i) = 0;
            for (int g1 = 0; g1 < nState; g1++) {
                bMat(g,i) += piMat(g,g1) * gMat(g1,(i+1)) * bMat(g1,(i+1));
            }
            bMat(g,i) /= sVec[i];
        }
    }
    
    return bMat;
}
