Calc_Roots = function(ra,rb,level) {
# computes Root distribution by vegetation type and returns desired cumilative total for a layer depth (cm)
# INPUT:
#	ra: parameter from Zeng (2001) 
#	rb: parameter from Zeng (2001) 
#	level: depth of cumilative layer, whole number with cm units
	
	z = seq(0,1, by=0.01)
	nz = length(z)

	roots = array(NA, c(nz))
	roots.cum = array(NA, c(nz))

	for(i in 2:nz) roots[i] = 0.5*(exp(-1*ra*z[(i-1)]) + exp(-1*rb*z[(i-1)]) - exp(-1*ra*z[i]) - exp(-1*rb*z[i]))

	roots[1] = 0
	
	for(i in 1:nz) roots.cum[i] = sum(roots[1:i])

	roots.surf = roots.cum[(level+1)]
	roots.surf
	
	}