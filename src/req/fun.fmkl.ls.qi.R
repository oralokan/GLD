"fun.fmkl.ls.qi" <-
function(n, param){
	
	#Qi does not work for moderately large n because of the gamma fns going to Inf
	
	L3 <- param[1]
	L4 <- param[2]
	
	g <- gamma
	i <- 1:n
	
	EPL3 <- (g(n+1) * g(i+L3)) / (g(i) * g(n+L3+1))
	EPL4 <- (g(n+1) * g(n-i+L4+1)) / (g(n-i+1) * g(n+L4+1))
	
	l = list(i, EPL3, EPL4)
	
	print(l)
	
	return(EPL3 - EPL4)
}