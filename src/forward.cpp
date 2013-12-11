
#include <Rcpp.h>

RcppExport SEXP cpp_forward( SEXP pi, SEXP g, SEXP pi0 ) { 

	using namespace Rcpp;

    Rcpp::NumericMatrix piMat(pi), gMat(g);
    Rcpp::NumericVector pi0Vec(pi0);
    
    int nState = gMat.nrow();
    int nObs = gMat.ncol();
    Rcpp::NumericMatrix aMat( (nState+1), nObs );
    // last row is to keep & return sVec
    Rcpp::NumericVector sVec( nObs );
    
    // initialization
    
    sVec[0] = 0;
    for (int g = 0; g < nState; g++) {
        aMat(g,0) = pi0Vec[g] * gMat(g,0);
        sVec[0] += aMat(g,0);
    }
    for (int g = 0; g < nState; g++) {
        aMat(g,0) = aMat(g,0) / sVec[0];
    }
    aMat(nState,0) = sVec[0];
    
    // iteration
    
    for (int i = 1; i < nObs; i++ ) {
        // update scaling factor
        
        sVec[i] = 0;        
        for (int g1 = 0; g1 < nState; g1++) {
            for (int g2 = 0; g2 < nState; g2++) {
                sVec[i] += aMat(g1,(i-1)) * piMat(g1,g2) * gMat(g2,i);
            }
        }
        
        // add s to the last row of aMat
        
        aMat(nState,i) = sVec[i];
        
        // estimate for forward
        
        for (int g = 0; g < nState; g++) {
            aMat(g,i) = 0;
            for (int g1 = 0; g1 < nState; g1++) {
                aMat(g,i) += aMat(g1,(i-1)) * piMat(g1,g);
            }
            aMat(g,i) *= gMat(g,i) / sVec[i];
        }
    }
    
    return aMat;
}
