import sys
import pyPamtra

inputFile="../test/referenceProfile/testProfiles.dat"

t = pyPamtra.pyPamtra()
t.readProfile(inputFile)

t.settings["data_path"]='/home/mech/models/pamtra/data/'



testNo = sys.argv[1]

if testNo == "1":


	t.runPamtra([24,90,150])
	
elif testNo == "2":


	t.settings["passive"]=False
	t.settings["active"]=True
	t.settings["EM_snow"]='liudb'
	t.settings["isnow_n0"]=0
	t.settings["liu_type"]=4


	t.runPamtra(35)

elif testNo == "3":
	t.settings["write_ascii"]=True
	t.settings["write_nc"]=True

	#t.settings["passive"]=True
	#t.settings["active"]=False
	#t.settings["EM_snow"]='liudb'
	#t.settings["isnow_n0"]=0
	#t.settings["liu_type"]=8
	
	t.runPamtra(35)
else:
	sys.exit("unknown test number "+testNo)
	
	
t.writeToNumpy('../test/tmp/pythontest'+testNo+'.pickle')

reference = pyPamtra.pyPamtra()
