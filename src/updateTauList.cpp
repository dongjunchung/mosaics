
#include <Rcpp.h>

RcppExport SEXP cpp_updateTauList( SEXP a, SEXP b, SEXP pi, SEXP g ) { 
		
	using namespace Rcpp;

	Rcpp::NumericMatrix aMat(a), bMat(b), piMat(pi), gMat(g);
	int nState = gMat.nrow();
	int nObs = gMat.ncol();
	
	Rcpp::List tauList( nObs );    
	double denom;
	
	for (int i = 0; i < (nObs-1); i++ ) {
	    denom = 0;
	    
	    // update each element
	    
	    Rcpp::NumericMatrix tauEach( nState, nState );
	    
	    for (int g1 = 0; g1 < nState; g1++) {            
	        for (int g2 = 0; g2 < nState; g2++) {
	            tauEach(g1,g2) = aMat(g1,i) * piMat(g1,g2) * gMat(g2,(i+1)) * bMat(g2,(i+1));
	            denom += tauEach(g1,g2);
	        }
	    }
	    
	    // normalize
	    
	    for (int g1 = 0; g1 < nState; g1++) {            
	        for (int g2 = 0; g2 < nState; g2++) {
	            tauEach(g1,g2) /= denom;
	        }
	    }
	    
	    // add updated element to the list
	    
	    tauList[i] = tauEach;
	}
	
	// null element in the end
	
	Rcpp::NumericMatrix tauEach( nState, nState );
	for (int g1 = 0; g1 < nState; g1++) {            
	    for (int g2 = 0; g2 < nState; g2++) {
	        tauEach(g1,g2) = 0;
	    }
	}
	tauList[nObs-1] = tauEach;
	
	return tauList;
}
