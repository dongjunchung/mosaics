// Calculate score function to obtain updated mu (Rcpp)

#include <Rcpp.h> 

RcppExport SEXP cpp_stack( SEXP S, SEXP E, SEXP minx, SEXP maxx ) { 
   using namespace Rcpp;
	
   Rcpp::IntegerVector SVec(S), EVec(E);
   int minxValue = Rcpp::as<int>(minx);
   int maxxValue = Rcpp::as<int>(maxx);
   
   Rcpp::NumericVector yvar( (maxxValue - minxValue + 1) );
   
   for (int i = 0; i < yvar.size(); i++ ) {
	   yvar[i] = 0;
   }
   
   for (int i = 0; i < SVec.size(); i++ ) {
	   for (int j = SVec[i]; j <= EVec[i]; j++ ) {
		   yvar[j-minxValue] += 1;
	   }
   }
   
   return yvar;
}
