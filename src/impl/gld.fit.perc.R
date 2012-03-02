#Percentile Method Algorithm Implementation
#For Generalized Lambda Distribution Parameter Estimation
#Galip Oral Okan
#06 January, 2012

#Based on the algorithm by Dr. Canan Gunes Corlu and Dr. Bahar Biller

#The rootSolve package is required to numerically solve a system of two nonlinear equations using the function 'multiroot'

require(rootSolve)

"gld.fit.perc" <- function(data, u=0.1){
	
	data <- sort(data)

	roe_1 <- gld.fit.perc.p(data, 0.5)
	roe_2 <- gld.fit.perc.p(data, 1 - u) - gld.fit.perc.p(data, u)
	roe_3 <- (gld.fit.perc.p(data,0.5) - gld.fit.perc.p(data,0.1))/(gld.fit.perc.p(data,1-u) - gld.fit.perc.p(data,0.5))
	roe_4 <- (gld.fit.perc.p(data,0.75) - gld.fit.perc.p(data,0.25)) / roe_2

	dif <- gld.fit.perc.dif

	f <- function(x)	c((dif(0.5, u, x[1]) - dif(0.5, 1-u, x[2])) / (dif(1-u, 0.5, x[1]) - dif(u, 0.5, x[2]))  - roe_3,
						(dif(0.75, 0.25, x[1]) - dif(0.25, 0.75, x[2])) / (dif(1-u, u, x[1]) - dif(u, 1-u, x[2]))  - roe_4)

	
	x <- multiroot(f, c(0.5,0.5))$root

	x <- c( (dif(1-u,u,x[1]) - dif(u,1-u,x[2])) / roe_2 ,   x)

	x <- c( roe_1 - (  dif(0.5,1,x[2]) - dif(0.5,1,x[3])  )  /x[1] ,  x)

	return(x)
}

"gld.fit.perc.p" <- function(data, p){
	
	x <- (length(data) + 1) * p

	ab <- x%%1
	r <- x - ab

	y <- data[r] + ab * (data[r+1] - data[r])

	return(y)
}

"gld.fit.perc.dif" <- function(x, y, z){
	
	return(( x^z - y^z ) / z )

}