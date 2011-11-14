# -*- coding: utf-8 -*-
from numpy import *

def evaluate(*args):
	global globals().keys()
	return asarray(eval(*args))


	#print "NO NUMEXP AVAILABLE"
	#from numpy import *
	##def pseudoNumexpr(function):
		##
		##return(np.asarray(eval(function)))
	#feval = eval
