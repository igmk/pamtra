import pyPamtra

inputFile="../test/referenceProfile/testProfiles.dat"

t = pyPamtra.pyPamtra()
t.readProfile(inputFile)


t.settings["output_path"]='../test/tmp'
t.settings["data_path"]='/home/mech/models/pamtra/data/'
t.settings["write_ascii"]=True
t.settings["write_nc"]=False

t.settings["passive"]=False
t.settings["active"]=True
t.settings["EM_snow"]='liudb'
t.settings["isnow_n0"]=0
t.settings["liu_type"]=4


t.runPamtra(35)

