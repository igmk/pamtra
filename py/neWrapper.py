# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np


try:
	import numexpr as ne
	feval = ne.evaluate

except:
	print "NO NUMEXP AVAILABLE"
	from numpy import *
	#def pseudoNumexpr(function):
		#global globals().keys()
		#return(np.asarray(eval(function)))
	feval = eval
