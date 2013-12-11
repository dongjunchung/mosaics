
conv_1S <- function( y, mu_round, mu_round_U, pN, pS ) {
	conv <- .Call( "convolution_1S_cpp", y, mu_round, mu_round_U, pN, pS, PACKAGE="mosaics" )
	return(conv)
}
