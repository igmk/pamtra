import pyPamtraLib
import os


def PamtraFortranWrapper(*args):
	#is needed because pp cannot work with fortran modules directly
	code = ""
	for ii in range(len(args)):
		code = code + "args["+str(ii)+"],"
	result = eval("pyPamtraLib.pypamtralib("+code[0:-1]+")")

	return result, os.uname()[1]
