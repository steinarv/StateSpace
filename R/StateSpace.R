loadModule("kalMod", TRUE)

rcpparma_hello_world <- function(){
	.Call( "rcpparma_hello_world", PACKAGE = "StateSpace" )
}

abso <- function(x){
	.Call( "abso", x, PACKAGE = "StateSpace" )
}