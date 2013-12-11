
conv_2S <- function( y, mu_round, mu_round_U, pN, pS1, pS2 ) {
	conv <- .Call( "convolution_2S_cpp", y, mu_round, mu_round_U, pN, pS1, pS2, PACKAGE="mosaics" )
	return(conv)
}
