import sys
import pyPamtra
import numpy as np

inputFile="../test/referenceProfile/testProfiles.dat"

t = pyPamtra.pyPamtra()
t.readProfile(inputFile)

t.set["data_path"]='/home/mech/models/pamtra/data/'
t.set["verbose"]=0


testNo = sys.argv[1]

if testNo == "1":


	t.runPamtra([24,90,150])
	
elif testNo == "2":

	t.set["passive"]=False
	t.set["active"]=True
	t.set["EM_snow"]='liudb'
	t.set["isnow_n0"]=0
	t.set["liu_type"]=4


	t.runPamtra(35)

elif testNo == "3":

	t.set["passive"]=True
	t.set["active"]=False
	t.set["EM_snow"]='liudb'
	t.set["isnow_n0"]=0
	t.set["liu_type"]=8
	
	t.runPamtra(35)
elif testNo == "4":

	t.set["EM_snow"]='icesf'
	t.set["isnow_n0"]=2
	t.set["EM_grau"]='icesf'
	t.runPamtra(35)
	t.writeResultsToNetCDF("../test/tmp/pythontest4.nc")
	
else:
	sys.exit("unknown test number "+testNo)
	

#uncomment if test should be defined again
#t.writeResultsToNumpy("../test/referenceOutput/python"+testNo+".pickle")

if testNo != "4":

	reference = pyPamtra.pyPamtra()
	reference.loadResultsFromNumpy("../test/referenceOutput/"+testNo+"/python"+testNo+".pickle")

	for key in ["Ze","attenuationAtmo","attenuationHydro","angles","tb","hgt"]:
		if np.any(reference.r[key] != t.r[key]):
			raise IOError(key+" does not match")
		#else:
			#print key, "OK"

			#for key in self.__dict__.keys():
				#print key,type(self.__dict__[key]),self.__dict__[key]
