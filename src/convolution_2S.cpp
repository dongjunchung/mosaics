
#include <Rcpp.h>

using namespace Rcpp;

RCPP_FUNCTION_6(NumericVector,convolution_2S_cpp,IntegerVector y,NumericVector mu,NumericVector mu_U,NumericMatrix pN,NumericVector pS1,NumericVector pS2) {
	
	// result vector
	
	NumericVector conv1(y.size());
	NumericVector conv2(y.size());
	NumericVector conv(2*y.size());
	int conv_index = 0;
	
	// convolution
	// - length of y = length of mu
	// - rows of pN match mu_U
	
	for ( int i=0; i<y.size(); i++ ) {
		conv1[i] = 0.0;
		conv2[i] = 0.0;
		for ( int j=0; j<mu_U.size(); j++ ) {
			if ( mu_U[j]==mu[i] ) {
				for ( int k=0; k<=y[i]; k++ ) {
					conv1[i] += pS1[(y[i]-k)] * pN(j,k);
					conv2[i] += pS2[(y[i]-k)] * pN(j,k);
				}
			}
		}
	}
	
	// summarize results
	
	for ( int i=0; i<y.size(); i++ ) {
		conv[conv_index] = conv1[i];
		conv_index++;
		conv[conv_index] = conv2[i];
		conv_index++;
	}
	
	return conv;

}
