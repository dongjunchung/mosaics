
#include <Rcpp.h>

RcppExport SEXP convolution_1S_cpp( 
	SEXP y, SEXP mu, SEXP mu_U, SEXP pN, SEXP pS ) { 

	using namespace Rcpp;
	
	// initialization
	
	Rcpp::IntegerVector yVec( y );
	Rcpp::NumericVector muVec( mu );
	Rcpp::NumericVector muUVec( mu_U );
	Rcpp::NumericMatrix pNMat( pN );
	Rcpp::NumericVector pSVec( pS );
	
	// result vector
	
	NumericVector conv( yVec.size() );
	
	// convolution
	// - length of y = length of mu
	// - rows of pN match mu_U
	
	for ( int i=0; i<yVec.size(); i++ ) {
		conv[i] = 0.0;
		for ( int j=0; j<muUVec.size(); j++ ) {
			if ( muUVec[j]==muVec[i] ) {
				for ( int k=0; k<=yVec[i]; k++ ) {
					conv[i] += pSVec[(yVec[i]-k)] * pNMat(j,k);
				}
			}
		}
	}
	
	return conv;

}
