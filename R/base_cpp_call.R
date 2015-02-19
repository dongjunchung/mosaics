
# MOSAiCS-HMM

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

# Stack fragments

.ff_stack <- function( S, E, minx, maxx ) {
    stopifnot( length(S) >= 1L, length(E) >= 1L, length(S) == length(E) )
    out <- .Call( "cpp_stack", S, E, minx, maxx, PACKAGE="mosaics" )
    return(out)
}
