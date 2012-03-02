"ks.test.gld.p" <-
function(D, n){
	
	p <- 0.001
	
	while(sqrt(-0.5 * log(p/2))/sqrt(n) > D[1]){
		p <- p + 0.001
	}
	
	if(p>1){
		p <- 1
	}
	
	return(p)
}