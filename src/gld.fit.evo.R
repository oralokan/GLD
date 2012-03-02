#Genetic Algorithm Implementation
#For Generalized Lambda Distribution Parameter Estimation
#Galip Oral Okan
#18 December, 2011

#Basic genetic algorithm parameters based on 
##	1.	Poojari CA, Lucas C, Mitra G (2007) Robust solutions and risk measures for a supply chain planning problem under uncertainty. J Oper Res Soc 59: 2–12. doi:10.1057/palgrave.jors.2602381


"gld.fit.evo" <- function(data, shape_init=c(-0.25, 1.5), p_cross=0.7, p_mut=0.1, p_repl=0.5, pop_size=20, max_gen=400, decay=0.9, tolerance=0.01){
	#This is the main method in which the Genetic Algorithm for GLD Parameter Estimation is defined.

	#	data		-	a set of data on which the GLD fitting will be done
	#	shape_init	-	the interval on which initial shape parameters will be generated
	#	p_cross		-	crossover probability
	#	p_mut		-	mutation probability
	#	p_repl		-	proportion replaced between generations
	#	pop_size	-	population size
	#	max_gen		-	maximum number of generations before algorithm terminates
	#	decay		-	explained under 'gld.fit.evo.selection'
	#	tolerance	-	goodness of fit value (KS) threshold of an acceptable solution


	data 		<- 	sort(data)
	data.mean	<-	mean(data)
	data.sd 	<-	sd(data)

	
	#POPULATION INITIALIZATION

	#The initial population is created by generating pairs of shape parameters from a Sobol
	#Sequence, and then by determining the corresponding location and scale parameters from
	#the equations derived in the definition of the moment matching method.

	#The first part is implemented in the method 'gld.fit.evo.initiate'
	#The second part is implemented in the method 'gld.fit.evo.complete'

	#The format of the resulting population is a (pop_size x 4) matrix.
	
	population <- gld.fit.evo.initiate(pop_size, shape_init, data.mean, data.sd)

	#Now, the population is ordered according to goodness-of-fit. The 'gld.fit.evo.fitness' method
	#returns the KS goodness-of-fit value for a given solution

	gof 		<- 	apply(population, 1, function(x) gld.fit.evo.fitness(x, data))
	population 	<-	population[order(gof), ]

	#INITIALIZATION OF ITERATION VARIABLES

	#best_sol represents the incumbent solution to the fitting problem. It is a vector with
	#length 5. The first four places represent the four parameters of the fit, and the last
	#place holds the goodness-of-fit (KS) value. Its value is initialized using the first
	#solution in the population, which is the best initial solution due to the ordering done
	#above.
	
	current_gen 	<-	0
	best_sol		<-	c(population[1,], gld.fit.evo.fitness(population[1,], data))


	#START OF MAIN ITERATION LOOP
	repeat{

		#This loop represents a generation in the evolutionary process.
		#There are two parts. The first part generates a number of new offspring (according to
		#p_repl - 'proportion replaced'). The rest of the new generation is filled by selecting
		#members from the old population. The selection algorithm is described under
		#'gld.fit.evo.selection'

		#There are two TERMINATING CONDITIONS.
		#	1 - An acceptable solution is found
		#	2 - The maximum number of generations has been reached

		current_gen <- current_gen + 1

		#Here, we initialize the matrix that will hold the next generation
		next_population <- matrix(0, pop_size, 4)

		#Initializing the number of members added to the new population so far to 0
		new_size <- 0

		#PART 1 - NEW OFFSPRING CREATION
		repeat{
			#A new offspring is created by first selecting an original parent from the population.
			#This prototype can then be subjected to crossover and mutation operations (with
			#probabilities as defined by p_cross and p_mut).
			
			#If crossover is to take place, a second member of the population is also selected to
			#act as the second parent.

			new_size <- new_size + 1

			offspring <- gld.fit.evo.selection(population, decay)

			#CROSSOVER
			if(runif(1) <= p_cross){
				p2 			<-	gld.fit.evo.selection(population, decay)	#Second parent selection
				offspring 	<-	gld.fit.evo.crossover(offspring, p2, data.mean, data.sd)
			}

			#MUTATION
			if(runif(1) <= p_mut){
				offspring 	<-	gld.fit.evo.mutation(offspring, data.mean, data.sd)
			}

			#Now, the newly created offspring is added to the new population
			next_population[new_size, ] <- offspring

			#Check if 'p_repl' has been reached, and if so, terminate new offspring creation
			if(new_size >= floor(pop_size * p_repl)) { break }
		}

		#PART 2 - OLD MEMBER RETENTION
		repeat{
			new_size <- new_size + 1

			#Choose members from the old population to be saved
			next_population[new_size, ] <- gld.fit.evo.selection(population, decay)

			#Check if the population size has reached 'pop_size', and if so, terminate.
			if(new_size >= pop_size) { break }
		}


		#Replace old population by new one and reorder it according to goodness-of-fit
		population 	<-	next_population
		gof 		<- 	apply(population, 1, function(x) gld.fit.evo.fitness(x, data))
		population 	<-	population[order(gof), ]


		#Update best_sol if necessary
		pop_best_gof <- gld.fit.evo.fitness(population[1,] ,data)
		if(pop_best_gof <= best_sol[5]){
			best_sol <- c(population[1,], pop_best_gof)
		}

		#Check for terminating conditions
		if(current_gen >= max_gen) { break }
		if(best_sol[5] <= tolerance) { break }

	} #END OF MAIN ITERATION LOOP

	#Return fit parameters
	return(best_sol[1:4])
}


"gld.fit.evo.initiate" <- function(pop_size, shape_init, data.mean, data.sd){
	#This function initializes the population as described above. The Sobol Sequence is generated
	#using the GLDEX method 'fun.gen.qrn'

	a 	<-	shape_init[1]
	b 	<-	shape_init[2] - shape_init[1]

	#Generating the shape parameters and holding them in a (pop_size x 2) matrix.
	pop.shape <- fun.gen.qrn(n=pop_size, dimension=2, scrambling=3, FUN="runif.sobol") *  b + a

	#Initializing the actual (pop_size x 4) matrix that will hold the initial population
	pop.init <- matrix(0, pop_size, 4)

	#Deriving the location and scale parameters for each generated pair of shape parameters
	#and storing them in a row of 'pop.init'
	for(i in 1:pop_size){
		pop.init[i,] <- gld.fit.evo.complete(pop.shape[i,], data.mean, data.sd)
	}

	return(pop.init)
}

"gld.fit.evo.complete" <- function(shape_params, data.mean, data.sd){
	#This function determines the location and scale parameters that correspond to a pair of shape
	#parameters under given (estimated) values for data mean and data variance.

	L3 <- shape_params[1]
	L4 <- shape_params[2]

	V1 = 1 / (L3 * (L3 + 1)) - 
		 1 / (L4 * (L4 + 1))

	V2 = 1 / (L3^2 * (2*L3 + 1)) +
		 1 / (L4^2 * (2*L4 + 1)) -
		 (2 / (L3 * L4)) * beta(L3 + 1, L4 + 1)

	L2 <- sqrt(V2 - V1^2) / data.sd
	L1 <- data.mean + (1/L2) * (1/(L3 + 1) - 1/(L4 + 1))

	return(cbind(L1,L2,L3,L4))	
}

"gld.fit.evo.fitness" <- function(x, data){
	#This function returns the KS goodness-of-fit value of the candidate solution represented by the
	#variable 'x'. The 'pgl' method requires the GLDEX package to function.
 
	n 	<-	length(data)
	i 	<-	c(1:n)

	F0 	<-	pgl(data, x[1], x[2], x[3], x[4])
	Fn 	<-	i/n

	D  	<-	max( abs( Fn-F0 ) )

	return(D)
}

"gld.fit.evo.selection" <- function(population, decay){
	#This method selects a member from the population and returns it. The method has a
	#higher probability of returning solutions that are have smaller index values.
	#Since the population is ordered according to goodness-of-fit, this correlates
	#to a higher probability of returning solutions with better goodness-of-fit.

	#The way this method works is to generate an exponential random variate with a fixed
	#probability of being greater than the population size (0.01 if decay = 1). If this
	#happens, a solution is selected randomly. Otherwise, the exponential variate determines
	#the solution to be selected.

	#i.e. The solution to be returned is selected using an exponential variate, if the
	#exponential variate yields a number smaller than the population size; and is selected
	#randomly otherwise.

	#Decreasing the 'decay' value reduces the steepness of the exponential decay and increases
	#the probability of selecting a solution randomly rather than depending on fitness.

	#This fixes the overshooting probability to 0.01 if decay = 1
	rate <-	decay * 4.60517 / nrow(population)	

	#Generate exponential variate
	r <- ceiling(rexp(1, rate))

	#If overshot population size, choose randomly
	if(r > nrow(population)) {r <- ceiling(runif(1, 0, nrow(population)))}

	return(population[r,])
}

"gld.fit.evo.crossover" <- function(pX, pY, data.mean, data.sd){
	#This function performs crossover between the two parents pX and pY by interpolating points
	#between the corresponding shape parameters, and then determining the other parameters.

	seedX <- pX[3:4]
	seedY <- pY[3:4]

	seedX[1] <- seedX[1] + runif(1) * (seedY[1] - seedX[1])
	seedX[2] <- seedX[2] + runif(1) * (seedY[2] - seedX[2])

	return(gld.fit.evo.complete(seedX, data.mean, data.sd))
}

"gld.fit.evo.mutation" <- function(cand, data.mean, data.sd, scale){
	#This function performs mutation by randomly increasing or decreasing the shape parameters
	#and then deriving the other parameters

	seed <- cand[3:4]

	seed[1] <- seed[1] + rnorm(1, 0, 0.01)
	seed[2] <- seed[2] + rnorm(1, 0, 0.01)

	return(gld.fit.evo.complete(seed, data.mean, data.sd))
}


