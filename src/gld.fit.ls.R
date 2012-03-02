#Least Squares Estimation Method Implementation
#For Generalized Lambda Distribution Parameter Estimation
#Galip Oral Okan
#January 10, 2012

#Based on the algorithm developed by Ozturk and Dale (1985)

"gld.fit.ls" <-
function(data, init = c(-0.25, 1.5)){

	#This is the main method in which the Genetic Algorithm for GLD Parameter Estimation is defined.

	#	data		-	a set of data on which the GLD fitting will be done
	#	init		-	the interval on which initial shape parameters will be generated

	data	<-	sort(data)
	n		<-	length(data)

	#GENERATION OF INITIAL SHAPE PARAMETER VALUES
	#This algorithm uses the Sobol sequence to generate values. Other methods can also be used. The Sobol sequence generation algorithm used has been defined in the GLDEX package written by Steven Su.

	#The variables 'e' and 'd' are just used for convenience in describing the interval defined by init. 'e' is the start of the interval and 'd' is the length of the interval.

	e	<-	init[1]
	d	<-	init[2] - e

	#This is where the actual generation takes place. It requires the GLDEX package to function

	g.init <- fun.gen.qrn(n=30, dimension=2, scrambling=3, FUN="runif.sobol") * d + e

	#Determine Q_i for all of the generated initial values
	qi <- apply(g.init, 1, function(x) gld.fit.ls.qi(n, x))

	#Determine the sqare error by regressing the data against Q_i lists corresponding to each initial solution
	sq_err <- apply(qi, 2, function(x) gld.fit.ls.obj(data, x))

	#Order the list of generated solutions in descending order of squares error.
	o <- order(-sq_err)
	init.sol <- apply(g.init, 2, function(x) x[o])

	#Use the Nelder-Mead Simplex Algorithm (builtin R function 'optim') to minimize squared error using the initial values determined above. 

	optim.result <- optim(init.sol, gld.fit.ls.obj.param, data=data, control=list(maxit=200000))

	#Return the optimal solution
	return(gld.fit.ls.param(optim.result$par, data))
}


"gld.fit.ls.qi" <-
function(n, param){
	#This function generates the list of Q_i as defined by Ozturk and Dale (1985).

	#Since the actual definition of Q_i requires solving a number of gamma functions (which get large quickly), there is a chance that several values may be calculated as Inf, preventing the function from returning a useable answer. Ozturk and Dale also present an approximation method (R_i) that can be used to avoid this problem. Since both 'n' and the parameters involved have an effect on whether the problem blows up or not, the approach taken is to first try and use Q_i - and then solve using R_i if NaN values are detected. However, it has been observed that NaN values come up no matter what the parameters are for 'n' values of roughly 100 or more. To prevent trying the Q_i approach when it is sure to fail (hence wasting computer resources) this step is omitted if 'n' is larger than 100.

	#Calculation of Q_i is omitted for n > 100
	omit_qi	<-	FALSE
	if(n > 100){ omit_qi <- TRUE }

	L3 	<-	param[1]
	L4	<-	param[2]

	i	<-	1:n

	if(!omit_qi){
		#This step calculates Q_i as defined by Ozturk and Dale (1985)

		g	<-	gamma

		EPL3 <- (g(n+1) * g(i+L3)) / (g(i) * g(n+L3+1))
		EPL4 <- (g(n+1) * g(n-i+L4+1)) / (g(n-i+1) * g(n+L4+1))

		qi	<-	EPL3 - EPL4
	}

	if(omit_qi || sum(is.nan(qi)) > 0){
		#If Q_i has been omitted, or has been calculated but includes NaN values, it is approximated using R_i as defined by Ozturk and Dale (1985)

		X1 <- (i/(n+1))^L3
		X2 <- ((n-i+1)/(n+1))^L4

		qi <- X1 - X2
	}

	return(qi)
}

"gld.fit.ls.obj" <- function(xi, qi){
	#This function returns the objective function value for a given data set and list of Q_i. This is the R^2 value of the regression between xi and qi.

	#Linear regression of xi as a function of qi
	reg <- lm(xi ~ qi)

	#Return the R^2 value of the regression above
	return(summary(reg)$r.squared)
}

"gld.fit.ls.obj.param" <- function(x, data){
	qi <- gld.fit.ls.qi(length(data), x)
	reg <- gld.fit.ls.obj(data, qi)
	return(-reg)
}

"gld.fit.ls.param" <- function(shape, data){
	#This function determines the location and scale parameters from the pair of shape parameters and the regression coefficients as described by Ozturk and Dale (1985)

	L3 <- shape[1]
	L4 <- shape[2]

	n <- length(data)

	qi <- gld.fit.ls.qi(n, shape)

	reg		<- lm(data ~ qi)
	coef	<- reg$coefficients[2]

	L2	<-	1 / coef
	L1	<-	mean(data) - coef * mean(qi)

	return(unname(c(L1, L2, L3, L4)))
}