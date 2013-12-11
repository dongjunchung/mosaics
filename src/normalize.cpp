
#include <Rcpp.h>

RcppExport SEXP cpp_normalize( SEXP pp0, SEXP pp1 ) { 

	using namespace Rcpp;

    Rcpp::NumericVector pp0Vec(pp0), pp1Vec(pp1); 
    int nObs = pp0Vec.size();       
    Rcpp::NumericVector pp0Final(nObs);
    
    // iteration
    
    for ( int i = 0; i < nObs; i++ ) {
        pp0Final[i] = pp0Vec[i] / ( pp0Vec[i] + pp1Vec[i] );
    }
    
    return pp0Final;
}
