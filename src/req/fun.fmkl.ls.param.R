"fun.fmkl.ls.param" <-
function(pair, data){
	L3 <- pair[1]
	L4 <- pair[2]
	
	n <- length(data)
	
	ri <- (fun.fmkl.ls.ri(n, pair))
	
	reg <- lm(data ~ ri)
	coef <- reg$coefficients[2]
	
	#For some reason L2 is half of the correct value when 1/coef
	L2 <- 1 / coef
	L1 <- mean(data) - coef * mean(ri)
	
	return(cbind(L1, L2, L3, L4))
	
}