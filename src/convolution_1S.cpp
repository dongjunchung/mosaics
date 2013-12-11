
#include <Rcpp.h>

using namespace Rcpp;

RCPP_FUNCTION_5(NumericVector,convolution_1S_cpp,IntegerVector y,NumericVector mu,NumericVector mu_U,NumericMatrix pN,NumericVector pS) { 
	
	// result vector
	
	NumericVector conv( y.size() );
	
	// convolution
	// - length of y = length of mu
	// - rows of pN match mu_U
	
	for ( int i=0; i<y.size(); i++ ) {
		conv[i] = 0.0;
		for ( int j=0; j<mu_U.size(); j++ ) {
			if ( mu_U[j]==mu[i] ) {
				for ( int k=0; k<=y[i]; k++ ) {
					conv[i] += pS[(y[i]-k)] * pN(j,k);
				}
			}
		}
	}
	
	return conv;

}
