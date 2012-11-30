import sys
import numpy as np

#ensure that you take the lib from the build directory!!
#sys.path.insert(0,"../py/lib/") #doesn't have impact on pp, so discarded
import pyPamtra



inputFile="../test/referenceProfile/testProfiles.dat"

t = pyPamtra.pyPamtra()
t.readPamtraProfile(inputFile)



#make artificially less levels for one profile!
t.p["nlyrs"][2,0] = 25

#change artificially one timestamp
t.p["unixtime"][0,0] = 600 #1970-01-01 00:10


testNo = sys.argv[1]

t.readNmlFile("../test/nmls/test"+testNo+".nml")
t.set["verbose"]=0
t.set["pyVerbose"]=1



#if testNo == "1":

	#pp_servers=()
	##t.runPamtra([24,90,150])
	#t.runParallelPamtra([24,90,150],pp_local_workers=2,pp_servers=pp_servers,pp_deltaF=1,pp_deltaX=1)
	
#elif testNo == "2":

	#t.runPamtra(35)
#elif testNo == "3":
	
	#t.runPamtra(35)
#elif testNo == "4":
	
	#t.runPamtra(35)
	#t.writeResultsToNetCDF("../test/tmp/pythontest4.nc")
	
#else:
	#sys.exit("unknown test number "+testNo)
	
#if testNo != "4":

	##uncomment if test should be defined again
	###
	##t.writeResultsToNumpy("../test/referenceOutput/"+testNo+"/python"+testNo+".pickle");print "warning, rewrting tests!!"


	#reference = pyPamtra.pyPamtra()
	#reference.loadResultsFromNumpy("../test/referenceOutput/"+testNo+"/python"+testNo+".pickle")
	##import pdb; pdb.set_trace()
	#error = 0

	#for key in ["angles","tb","hgt", "Ze_rr", "Ze_ci", "Ze_sn", "Ze_gr", "Ze_ha", "Att_cw", "Att_rr", "Att_ci", "Att_sn", "Att_gr", "Att_ha", "Att_atmo"]:
		#if np.any(reference.r[key] != t.r[key]):
			#error += 1
			#print key, "max. difference:", np.max(reference.r[key] - t.r[key])
	#if error > 0:
		##import pdb;pdb.set_trace()
		#raise IOError("do not match")
	##else:
		##print key, "OK"

		##for key in self.__dict__.keys():
			##print key,type(self.__dict__[key]),self.__dict__[key]

