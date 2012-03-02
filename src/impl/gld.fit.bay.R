#Bayesian Algorithm Implementation
#For Generalized Lambda Distribution Parameter Estimation
#Galip Oral Okan
#06 January, 2012

#Based on the algorithm by Dr. Canan Gunes Corlu and Dr. Bahar Biller

"gld.fit.bay" <- function(data, tolerance=1, loc_int=c(-1,1), scal_int=c(0.001, 1), shap_int1=c(-2,2), shap_int2=c(-2,2), max_iter=10000){
	
	#This is the main method in which the Bayesian Method Algorithm developed by Dr. Corlu and Dr. Biller is implemented.

	#	data		-	a set of data on which the GLD fitting will be done
	#	tolerance	-	fitness function 
	#	loc_int		-	interval for location parameter
	#	scal_int	-	interval for scale parameter
	#	shap_int	-	interval for first shape parameter
	#	shap_int2	-	interval for second shape parameter
	#	max_iter	-	max number of iterations allowed

	data	<- 	sort(data)
	n		<-	length(data)

	iter		<-	0
	sol.found	<-	FALSE

	#The initial solution is arbitrarily defined with a large fitness value. The format is [1, 2, 3, 4, 9999] where the first four values correspond to parameter values and the last is the fitness value that results.

	sol		<-	c(1:4)
	sol		<- c(sol, 9999)


	#This is the main loop that iterates until a solution is found (fitness is smaller than tolerance), or the maximum number of iterations has been reached
	while(iter < max_iter && !sol.found){
		
		#First, a set of parameters is generated uniformly from the given parameter intervals.
		p <- c( runif(1, loc_int[1], loc_int[2]),
				runif(1, scal_int[1], scal_int[2]),
				runif(1, shap_int1[1], shap_int1[2]),
				runif(1, shap_int2[1], shap_int2[2]) )

		
		#This set of parameters is used to simulate a set of random values from a GLD distribution that has these parameters.
		d <- gld.fit.bay.rgld(n, p[1], p[2], p[3], p[4])
		d <- sort(d)

		#The original data is compared to the simulated data to determine a fitness value.
		t <- gld.fit.bay.fitness(data, d)

		#If the current fitness value is better than that of the incumbent solution, the current solution become the new incumbent solution.
		if(t < sol[5]){
			sol <- c(p, t)

			#If the current solution has just becom the incumbent solution, it is tested to see whether its fitness is better than the tolerance parameter.
			if(t < tolerance){
				sol.found <- TRUE
			}
		}
		iter <- iter + 1
	}

	#The parameter values of the solution are returned (the fitness value is truncated)
	return(sol[1:4])
}

"gld.fit.bay.rgld" <- function(n, lam_1, lam_2, lam_3, lam_4){
	#Returns 'n' random variats distributed according to GLD distribution with given parameters

	u = runif(n)
	return(lam_1 + ( (u^lam_3 - 1)/lam_3 - ((1-u)^lam_4 - 1)/lam_4 )/lam_2)
}

"gld.fit.bay.fitness" <- function(d1, d2){
	res <- d1 - d2
	res <- res^2

	return(sqrt(sum(res)))
}