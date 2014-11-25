
#include <Rcpp.h>

RcppExport SEXP convolution_2S_cpp( 
	SEXP y, SEXP mu, SEXP mu_U, SEXP pN, SEXP pS1, SEXP pS2 ) { 

	using namespace Rcpp;
	
	// initialization
	
	Rcpp::IntegerVector yVec( y );
	Rcpp::NumericVector muVec( mu );
	Rcpp::NumericVector muUVec( mu_U );
	Rcpp::NumericMatrix pNMat( pN );
	Rcpp::NumericVector pS1Vec( pS1 );
	Rcpp::NumericVector pS2Vec( pS2 );
	
	// result vector
	
	Rcpp::NumericVector conv1(yVec.size());
	Rcpp::NumericVector conv2(yVec.size());
	Rcpp::NumericVector conv(2*yVec.size());
	int conv_index = 0;
	
	// convolution
	// - length of y = length of mu
	// - rows of pN match mu_U
	
	for ( int i=0; i<yVec.size(); i++ ) {
		conv1[i] = 0.0;
		conv2[i] = 0.0;
		for ( int j=0; j<muUVec.size(); j++ ) {
			if ( muUVec[j]==muVec[i] ) {
				for ( int k=0; k<=yVec[i]; k++ ) {
					conv1[i] += pS1Vec[(yVec[i]-k)] * pNMat(j,k);
					conv2[i] += pS2Vec[(yVec[i]-k)] * pNMat(j,k);
				}
			}
		}
	}
	
	// summarize results
	
	for ( int i=0; i<yVec.size(); i++ ) {
		conv[conv_index] = conv1[i];
		conv_index++;
		conv[conv_index] = conv2[i];
		conv_index++;
	}
	
	return conv;

}
