import sys
import pyPamtra
import numpy as np

inputFile="../test/referenceProfile/testProfiles.dat"

t = pyPamtra.pyPamtra()
t.readProfile(inputFile)

t.settings["data_path"]='/home/mech/models/pamtra/data/'
t.settings["verbose"]=0


testNo = sys.argv[1]

if testNo == "1":


	t.runPamtra([24,90,150])
	
elif testNo == "2":

	t.settings["passive"]=False
	t.settings["active"]=True
	#t.settings["EM_snow"]='liudb'
	#t.settings["isnow_n0"]=0
	#t.settings["liu_type"]=0


	t.runPamtra(35)

elif testNo == "3":

	t.settings["passive"]=True
	t.settings["active"]=False
	#t.settings["EM_snow"]='liudb'
	#t.settings["isnow_n0"]=0
	#t.settings["liu_type"]=3
	
	t.runPamtra(35)
else:
	sys.exit("unknown test number "+testNo)
	

#uncomment if test should be defined again
#t.writeResultsToNumpy("../test/referenceOutput/python"+testNo+".pickle")

reference = pyPamtra.pyPamtra()
reference.loadResultsFromNumpy("../test/referenceOutput/python"+testNo+".pickle")

for key in ["Ze","attenuationAtmo","attenuationHydro","angles","tb","hgt"]:
	if np.any(reference.r[key] != t.r[key]):
		raise IOError(key+" does not match")
	#else:
		#print key, "OK"

		#for key in self.__dict__.keys():
			#print key,type(self.__dict__[key]),self.__dict__[key]
