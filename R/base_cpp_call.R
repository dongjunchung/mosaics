
.ff_forward <- function( pi, g, pi0 ) {
    out <- .Call( "cpp_forward", pi, g, pi0, PACKAGE="mosaics" )
    return(out)
}
    
.ff_backward <- function( pi, g, s ) {
    out <- .Call( "cpp_backward", pi, g, s, PACKAGE="mosaics" )
    return(out)
} 

.ff_updateTauList <- function( a, b, pi, g ) {
    out <- .Call( "cpp_updateTauList", a, b, pi, g, PACKAGE="mosaics" )
    return(out)
} 
    
.ff_updateTauMat <- function( tau, g, pi0 ) {
    out <- .Call( "cpp_updateTauMat", tau, g, pi0, PACKAGE="mosaics" )
    return(out)
} 
    
.ff_updatePiMat <- function( tauL, tauM ) {
    out <- .Call( "cpp_updatePiMat", tauL, tauM, PACKAGE="mosaics" )
    return(out)
} 
    
.ff_normalize <- function( pp0, pp1 ) {
    out <- .Call( "cpp_normalize", pp0, pp1, PACKAGE="mosaics" )
    return(out)
} 
    
.ff_viterbi <- function( pi, g, pi0 ) {
    out <- .Call( "cpp_viterbi", pi, g, pi0, PACKAGE="mosaics" )
    return(out)
} 
