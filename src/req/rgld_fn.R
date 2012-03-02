rgld <- function(n, lam_1, lam_2, lam_3, lam_4){
	u = runif(n)
	return(lam_1 + ( (u^lam_3 - 1)/lam_3 - ((1-u)^lam_4 - 1)/lam_4 )/lam_2)
}