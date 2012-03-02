"gld.test.method" <-
function(method, n, params, runs=500){
	
	
	res <- c(1:runs)
	
	for(i in 1:runs){
		
		data <- rgl(n, params[1], params[2], params[3], params[4])
		
		sol <- do.call(method, list(data))
	
		res[i] <- ks.test.gld(data, params)
		
	}
	
	return(cbind(mean(res), var(res)))
}