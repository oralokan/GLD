"fun.RMFMKL.ls" <-
function(data, fmkl.init = c(-0.25, 1.5), leap=3, FUN="runif.sobol", no=100){
	
	data <- sort(data)
	n <- length(data)
	
	e <- fmkl.init[1]
	d <- fmkl.init[2] - e
	
	ncol.init <- 2
	
	#Generate the initial L3, L4 values
	g.init <- fun.gen.qrn(n=no, dimension=ncol.init, scrambling=leap, FUN=FUN) * d + e
	
	ri <- apply(g.init, 1, function(x) fun.fmkl.ls.ri(n, x))
		
	reg <- apply(ri, 2, function(x) fun.fmkl.ls.obj(data, x))
	
	t <- cbind("L3"=g.init[,1], "L4"=g.init[,2], "R@2"=reg)
	
	#Order descending according to R^2
	o <- order(-t[,3])
	t <- apply(t, 2, function(x) x[o])
	
	init.sol <- t[1,1:2]
	
	optim.result <- optim(init.sol, fun.fmkl.ls.obj.param, data=data, control=list(maxit=200000))

	return(fun.fmkl.ls.param(optim.result$par, data))
}