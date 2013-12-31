
.ff_emHMM <- function( pi, g, pi0, maxi, eps ) {
	out <- .Call( "cpp_emHMM", pi, g, pi0, maxi, eps, PACKAGE="mosaics" )
}
    
.ff_updatePiMat <- function( zL ) {
    out <- .Call( "cpp_updatePiMat", zL, PACKAGE="mosaics" )
    return(out)
} 
    
.ff_viterbi <- function( pi, g, pi0 ) {
    out <- .Call( "cpp_viterbi", pi, g, pi0, PACKAGE="mosaics" )
    return(out)
} 
