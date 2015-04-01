import pyPamtra
import matplotlib.pyplot as plt


fnames = "/data/gop/final/gop9/cols_lmk_jue/cols_jue201304/20130425_gop9_lmk_Selhausen_7x7.nc.gz"
kind = "gop_collumn"

descriptorFile = "../descriptorfiles/descriptor_file_COSMO_1mom.txt"
pam = pyPamtra.importer.readCosmoDe1MomDataset(fnames,kind,descriptorFile,forecastIndex = 1,colIndex=0,tmpDir="/tmp/",fnameInTar="",concatenateAxis=1,debug=False,verbosity=0)
pam.nmlSet["data_path"] = '/net/karif//mmaahn/pamtra_data/'

pam.runPamtra(35)

plt.figure()
pyPamtra.plot.plotTB(pam,levels=[1,400])