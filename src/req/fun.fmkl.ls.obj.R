"fun.fmkl.ls.obj" <-
function(xi, qi){
	reg <- lm(xi ~ qi)
	return(summary(reg)$r.squared)
}