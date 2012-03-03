#TESTING ALGORITHM

start_time <- proc.time()

ord = TRUE

#NUMBER OF RUNS
runs = 1


#METHODS TO BE TESTED

methods = matrix(ncol=1)
methods = rbind(methods, "fun.RMFMKL.lm")
methods = rbind(methods, "fun.RMFMKL.ml")
methods = rbind(methods, "fun.RMFMKL.qs")
methods = rbind(methods, "fun.RMFMKL.ls")
methods = rbind(methods, "genetic")
methods = rbind(methods, "bayesian")
methods = rbind(methods, "percentile")
methods = methods[2:nrow(methods),]


#DATA SIZES
data.sizes = c(25, 50, 100, 200, 400)
#data.sizes = 25
#L3, L4 PAIRS
params.shape = matrix(ncol=2)
#										L3		L4
params.shape = rbind(params.shape, c(	0.4,	0.4	))
params.shape = rbind(params.shape, c(	0.6, 	0.1	))
params.shape = rbind(params.shape, c(	0.8, 	0.8	))
params.shape = rbind(params.shape, c(	2,		0.1	))
params.shape = rbind(params.shape, c(	1.5,	1.2	))
params.shape = rbind(params.shape, c(	2.5,	1.5	))
params.shape = rbind(params.shape, c(	2.5,	2.5	))
params.shape = params.shape[2:nrow(params.shape),]

#L1, L2 PAIRS
rows = nrow(params.shape)
params.location = matrix(data=0, nrow=rows)
params.scale = matrix(data=1, nrow=rows)
rm(rows)

#PARAMETER MATRIX
params = cbind(params.location, params.scale, params.shape)
rm(params.shape, params.location, params.scale)

#INITIALIZE RESULTS MATRIX
res.rows = length(methods) * length(data.sizes) * nrow(params)
res.init.char = matrix("na", nrow=res.rows)
res.init.num = matrix(0, nrow=res.rows)
results <- data.frame(method=res.init.char, n=res.init.num, l1=res.init.num, l2=res.init.num, l3=res.init.num, l4=res.init.num, gof=res.init.num, gof_var=res.init.num, p_val=res.init.num,  stringsAsFactors = FALSE)

#Create Progress Bar
pb <- txtProgressBar(min=0, max=res.rows, style=3)


#TEST

i = 0
for(m in 1:length(methods)){
	for(n in 1:length(data.sizes)){
		for(p in 1:nrow(params)){			
			i = i + 1
			
			result.new = gld.test.method(methods[m], data.sizes[n], params[p,], runs)
						
			results[i, 1] = methods[m]
			results[i, 2] = data.sizes[n]
			results[i, 3] = params[p,1]
			results[i, 4] = params[p,2]
			results[i, 5] = params[p,3]
			results[i, 6] = params[p,4]
			results[i, 7] = result.new[1]
			results[i, 8] = result.new[2]
			results[i, 9] = ks.test.gld.p(result.new, data.sizes[n])
			
			setTxtProgressBar(pb, i)
		}
	}
}

close(pb)

if(ord){
	o <- order(results[,"gof"])
	results <- results[o,]
}

print("TEST COMPLETE")
cat("Time: ", (proc.time() - start_time)[3])

print(results)
