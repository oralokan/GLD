"ks.test.gld" <-
function(data, params){
	
	data <- sort(data)
	n <- length(data)

	i <- c(1:n)	
	F0 <- pgl(data, params[1], params[2], params[3], params[4])
	Fn <- i/n
	
	D <- max(abs(Fn-F0))
	
	#d <- data.frame(i, data, F0, Fn, Fn_1, D_plus, D_minus)
	#print(d)
	
	return(D)
}