
#include <Rcpp.h>

RcppExport SEXP cpp_updatePiMat( SEXP zL ) { 

	using namespace Rcpp;

    Rcpp::List zList(zL);
	Rcpp::NumericMatrix zCurrent = zList[0];
    int nState = zCurrent.nrow();
    int nObs = zList.size();
    double denom = 0.0;
    Rcpp::NumericMatrix piMat( nState, nState );
    
    // calculate pi matrix
    
    for (int g1 = 0; g1 < nState; g1++) {   
        for (int g2 = 0; g2 < nState; g2++) {   
            piMat(g1,g2) = 0;
            for (int i = 0; i<(nObs-1); i++ ) {
                Rcpp::NumericMatrix zCurrent = zList[i];
                piMat(g1,g2) += zCurrent(g1,g2);
            }
        }
    }
	
	// normalize: each row should sum to one
    
    for (int g1 = 0; g1 < nState; g1++) {   
        denom = 0.0;
		for (int g2 = 0; g2 < nState; g2++) { 
            denom += piMat(g1,g2);
        }
		for (int g2 = 0; g2 < nState; g2++) { 
            piMat(g1,g2) /= denom;
        }
    }
    
    return piMat;
}
