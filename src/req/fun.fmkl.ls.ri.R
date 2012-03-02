"fun.fmkl.ls.ri" <-
function(n, param){
	
	L3 <- param[1]
	L4 <- param[2]
	
	i <- 1:n
	
	X1 <- (i/(n+1))^L3
	X2 <- ((n-i+1)/(n+1))^L4
	
	return (X1 - X2)
	
}
