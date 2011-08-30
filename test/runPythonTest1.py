import pyPamtra

inputFile="../test/referenceProfile/testProfiles.dat"

t = pyPamtra.pyPamtra()
t.readProfile(inputFile)

t.settings["output_path"]='../test/tmp'
t.settings["data_path"]='/home/mech/models/pamtra/data/'
t.settings["write_ascii"]=True
t.settings["write_nc"]=False

t.runPamtra([24,90,150])
