
#include <Rcpp.h>

RcppExport SEXP cpp_emHMM( SEXP pi, SEXP g, SEXP pi0, SEXP maxi, SEXP eps ) { 

	using namespace Rcpp;
	
	int maxIter = as<int>( maxi );
	double epsPP0 = as<double>( eps );

    Rcpp::NumericMatrix piMat(pi);
	Rcpp::NumericMatrix gMat(g);
    Rcpp::NumericVector pi0Vec(pi0);

    int nState = gMat.nrow();
    int nObs = gMat.ncol();
	
    Rcpp::NumericMatrix aMat( nState, nObs );
		// forward
    Rcpp::NumericMatrix bMat( nState, nObs );
		// backward
    Rcpp::NumericVector sVec( nObs );
		// normalizer
	Rcpp::NumericVector z0Vec( nState );
	Rcpp::NumericMatrix zMat( nState * nState, (nObs-1) );
		// E step
	Rcpp::NumericVector pp0Vec( nObs );
		// P( Z = 0 | Y )
	Rcpp::NumericVector pp0VecPre( nObs );
		// stopping rule
	
	double denom = 0.0;
	double diffCur = 0.0;
	Rcpp::NumericMatrix zEach( nState, nState );
	Rcpp::NumericVector loglik( maxIter );
		// track changes in log likelihood
	Rcpp::NumericVector maxDiff( maxIter );
		// track changes in pp0
	
	for ( int iter = 0; iter < maxIter; iter++ ) {
		
		//Rcout << iter << std::endl;
		
		/////////////////////////////////////////////
		//                                         //
		//		forward algorithm                  //
		//                                         //
		/////////////////////////////////////////////
		
	
		// forward: initialization
		
		sVec[0] = 0;
		for (int g = 0; g < nState; g++) {
			aMat(g,0) = pi0Vec[g] * gMat(g,0);
			sVec[0] += aMat(g,0);
		}
		for (int g = 0; g < nState; g++) {
			aMat(g,0) = aMat(g,0) / sVec[0];
		}
		
		// forward: iteration
		
		for (int i = 1; i < nObs; i++ ) {
			// update scaling factor
			
			sVec[i] = 0;        
			for (int g1 = 0; g1 < nState; g1++) {
				for (int g2 = 0; g2 < nState; g2++) {
					sVec[i] += aMat(g1,(i-1)) * piMat(g1,g2) * gMat(g2,i);
				}
			}
			
			// estimate for forward
			
			for (int g = 0; g < nState; g++) {
				aMat(g,i) = 0;
				for (int g1 = 0; g1 < nState; g1++) {
					aMat(g,i) += aMat(g1,(i-1)) * piMat(g1,g) * gMat(g,i);
				}
				aMat(g,i) /= sVec[i];
			}
		}
		
		
		/////////////////////////////////////////////
		//                                         //
		//		backward algorithm                 //
		//                                         //
		/////////////////////////////////////////////
		
		
		// backward: initialization
		
		for (int g = 0; g < nState; g++) {
			bMat(g,(nObs-1)) = 1;
		}
		
		// backward: iteration
		
		for (int i = nObs-2; i>=0; i-- ) {
			for (int g = 0; g < nState; g++) {
				bMat(g,i) = 0;
				for (int g1 = 0; g1 < nState; g1++) {
					bMat(g,i) += piMat(g,g1) * gMat(g1,(i+1)) * bMat(g1,(i+1));
				}
				bMat(g,i) /= sVec[i];
			}
		}
		
		
		/////////////////////////////////////////////
		//                                         //
		//		E step                             //
		//                                         //
		/////////////////////////////////////////////
        
		// update prob for the first bin
		
		denom = 0.0;
		for (int g = 0; g < nState; g++) {
			z0Vec[g] = pi0Vec[g] * gMat(g,0);
			denom += z0Vec[g];
		}
		for (int g = 0; g < nState; g++) {
			z0Vec[g] /= denom;
		}
		
		// update transition
	
		for (int i = 0; i < (nObs-1); i++ ) {
			denom = 0.0;
			
			// update each element
			
			for (int g1 = 0; g1 < nState; g1++) {            
				for (int g2 = 0; g2 < nState; g2++) {
					zEach(g1,g2) = aMat(g1,i) * piMat(g1,g2) * gMat(g2,(i+1)) * bMat(g2,(i+1));
					denom += zEach(g1,g2);
				}
			}
			
			// normalize
			
			for (int g1 = 0; g1 < nState; g1++) {            
				for (int g2 = 0; g2 < nState; g2++) {
					zEach(g1,g2) /= denom;
				}
			}
			
			// update transition probabilities
			
			for (int g1 = 0; g1 < nState; g1++) {            
				for (int g2 = 0; g2 < nState; g2++) {
					zMat( (nState*g1+g2), i ) = zEach(g1,g2);
				}
			}
		}
		
		
		/////////////////////////////////////////////
		//                                         //
		//		M step                             //
		//                                         //
		/////////////////////////////////////////////
          
		// update pi0Vec
		
		denom = 0.0;
		for (int g = 0; g < nState; g++) {
			pi0Vec[g] = z0Vec[g];
			denom += pi0Vec[g];
		}
		for (int g = 0; g < nState; g++) {
			pi0Vec[g] /= denom;
		}
    
		// calculate pi matrix
		
		for (int g1 = 0; g1 < nState; g1++) {   
			for (int g2 = 0; g2 < nState; g2++) {   
				piMat(g1,g2) = 0;
				for (int i = 0; i<(nObs-1); i++ ) {
					piMat(g1,g2) += zMat( (nState*g1+g2), i );
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
		   
		//Rcout << pi0Vec[0] << " " << pi0Vec[1] << std::endl;
		//Rcout << piMat(0,0) << " " << piMat(0,1) << std::endl;
		//Rcout << piMat(1,0) << " " << piMat(1,1) << std::endl;
		
		
		/////////////////////////////////////////////
		//                                         //
		// posterior probability: P( Z_j = 0 | Y ) //
		//                                         //
		/////////////////////////////////////////////
                
        for (int i = 0; i < nObs; i++ ) {
			denom = 0.0;
			for (int g = 0; g < nState; g++) {
				denom += aMat(g,i) * bMat(g,i);
			}
			pp0Vec[i] = aMat(0,i) * bMat(0,i) / denom;
		}
		
		
		/////////////////////////////////////////////
		//                                         //
		//		stopping rule                      //
		//                                         //
		/////////////////////////////////////////////     

		loglik[iter] = 0.0;
		for (int i = 0; i < nObs; i++ ) {
			loglik[iter] += log(sVec[i]);
		}
		
		//Rcout << loglik[iter] << std::endl;
				
		// do we stop?
		
		maxDiff[iter] = 0.0;
		for (int i = 0; i < nObs; i++ ) {
			diffCur = pp0Vec[i] - pp0VecPre[i];
			if ( diffCur < 0 ) {
				diffCur = - diffCur;
			}
			if ( maxDiff[iter] < diffCur ) {
				maxDiff[iter] = diffCur;
			}
		}
		
		//Rcout << maxDiff[iter] << std::endl;
        
        if ( maxDiff[iter] < epsPP0 ) {
            break;
        }
		
		// update pp0
		  
        for (int i = 0; i < nObs; i++ ) {
			pp0VecPre[i] = pp0Vec[i];
		}
	}
	
	// output
	
	return Rcpp::List::create(
		Rcpp::Named( "aMat" ) = aMat,
        Rcpp::Named( "bMat" ) = bMat,
		Rcpp::Named( "sVec" ) = sVec,
		Rcpp::Named( "piMat" ) = piMat,
		Rcpp::Named( "pi0Vec" ) = pi0Vec,
		Rcpp::Named( "loglik" ) = loglik,
		Rcpp::Named( "maxDiff" ) = maxDiff
	);
}
