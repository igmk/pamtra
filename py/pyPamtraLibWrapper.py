import pyPamtraLib
import os
import sys



def PamtraFortranWrapper(*args):
	#this wrapper is needed because pp cannot work with fortran modules directly. returns results from pamtra AND name of the host (for pp statistics)

	result = pyPamtraLib.pypamtralib(*args)

	return result, os.uname()[1]
