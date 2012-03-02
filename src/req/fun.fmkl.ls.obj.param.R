"fun.fmkl.ls.obj.param" <-
function(x, data){
	
	ri <- fun.fmkl.ls.ri(length(data), x)
	reg <- fun.fmkl.ls.obj(data, ri)
	
	return(-reg)
}