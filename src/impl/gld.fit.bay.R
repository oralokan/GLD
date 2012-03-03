#Bayesian Algorithm Implementation
#For Generalized Lambda Distribution Parameter Estimation
#Galip Oral Okan
#06 January, 2012

#Based on the algorithm by Dr. Canan Gunes Corlu and Dr. Bahar Biller

"gld.fit.bay" <- function(data, tolerance=50, a=c(-1,0.001,-2,-2), b=c(1,1,2,2), s=c(0.1,0.1,0.1,0.1), K=10000){
	
	#This is the main method in which the Bayesian Method Algorithm developed by Dr. Corlu and Dr. Biller is implemented.

	#	data		-	a set of data on which the GLD fitting will be done
	#	tolerance	-	fitness function 
	#	a			-	lower bounds for parameters
	#	b			-	upper bounds for parameters
	#	s			-	std dev for norm
	#	K			-	max number of iterations

	data	<- 	sort(data)
	n		<-	length(data)

	init_sol.found	<- FALSE

	while(!init_sol.found){

		#First, a set of parameters is generated uniformly from the given parameter intervals.
		p <- c( runif(1, a[1], b[1]),
				runif(1, a[2], b[2]),
				runif(1, a[3], b[3]),
				runif(1, a[4], b[4]) )


		#This set of parameters is used to simulate a set of random values from a GLD distribution that has these parameters.
		d <- gld.fit.bay.rgld(n, p[1], p[2], p[3], p[4])
		d <- sort(d)

		#The original data is compared to the simulated data to determine a fitness value.
		fit <- gld.fit.bay.fitness(data, d)

		if(fit < tolerance){
			init_sol.found <- TRUE

			#Note that since the while loop will not iterate again, the final value of 'p' is our initial solution.
		}
	}

	#Now we initialize the matrix that will hold the iterations and insert the initial solution
	iter <- matrix(0, K+1, 4)
	iter[1,] <- p

	# MAIN ALGORITHM

	count <- 2

	while(count < K+2){
		iter[count,] <- c( rnorm(1,iter[count-1,1], s[1]*s[1]) , rnorm(1,iter[count-1,2], s[2]*s[2]) , rnorm(1,iter[count-1,3], s[3]*s[3]) , rnorm(1,iter[count-1,4], s[4]*s[4]) )

		#Truncating
		for(j in 1:4){
			if(iter[count,j] < a[j]) iter[count,j] <- a[j]
			else if(iter[count,j] > b[j]) iter[count,j] <- b[j]
		}

		d <- gld.fit.bay.rgld(n, iter[count, 1], iter[count, 2], iter[count, 3], iter[count, 4])
		d <- sort(d)

		u <- runif(1)

		proceed <- TRUE

		for(j in 1:4){
			t <- (runif(1, a[j], b[j]) * iter[count,j]   * rnorm(1, iter[count-1,j], s[j]*s[j])) /
				( runif(1, a[j], b[j]) * iter[count-1,j] * rnorm(1, iter[count,j], s[j]*s[j]))

			if (u > t) proceed <- FALSE
		}

		fit <- gld.fit.bay.fitness(data, d)

		if(fit > tolerance) proceed <- FALSE

		if(proceed) count <- count + 1
	}

	#replace this with return(iter[K+1, ]) to get actual estimate
	return(iter)
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