# iptrack - interpolate track
#
# SYNTAX
# p=iptrack(filename)
#
# INPUT
# filename: data file containing exported tracking data on the standard
# Tracker export format
#
# mass_A
# t	x	y
# 0.0	-1.0686477620876644	42.80071293284619
# 0.04	-0.714777136706708	42.62727536827738
# ...
#
# OUTPUT
# p=iptrack(filename) returns the coefficients of a polynomial of degree 15
# that is the least square fit to the data y(x). Coefficients are given in
# descending powers.

import numpy as np

def iptrack(filename, degree):
	data=np.loadtxt(filename,skiprows=2)
	return np.polyfit(data[:,1],data[:,2],degree)

# trvalues - track values
#
# SYNTAX
# [y,dydx,d2ydx2,alpha,R]=trvalues(p,x)
#
# INPUT
# p: the n+1 coefficients of a polynomial of degree n, given in descending
# order. (For instance the output from p=iptrack(filename).)
# x: ordinate value at which the polynomial is evaluated.
#
# OUTPUT
# [y,dydx,d2ydx2,alpha,R]=trvalues(p,x) returns the value y of the
# polynomial at x, the derivative dydx and the second derivative d2ydx2 in
# that point, as well as the slope alpha(x) and the radius of the
# osculating circle. 
# The slope angle alpha is positive for a curve with a negative derivative. 
# The sign of the radius of the osculating circle is the same as that of 
# the second derivative.

import numpy as np

def trvalues(p,x):
	y=np.polyval(p,x)
	dp=np.polyder(p)
	dydx=np.polyval(dp,x)
	ddp=np.polyder(dp)
	d2ydx2=np.polyval(ddp,x)
	alpha=np.arctan(-dydx)
	R=(1.0+dydx**2)**1.5/d2ydx2
	return [y,dydx,d2ydx2,alpha,R]
