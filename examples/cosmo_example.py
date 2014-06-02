import pyPamtra
import pyPamtraLibWrapper
import pyPamtraImport
import pyPamtraLib
import pyPamtraPlotter

reload(pyPamtraPlotter)

reload(pyPamtraImport)
reload(pyPamtra)


fnames = "/data/gop/final/gop9/cols_lmk_jue/cols_jue201304/20130425_gop9_lmk_Selhausen_7x7.nc.gz"
kind = "gop_collumn"

pam = pyPamtraImport.readCosmoDe1MomDataset(fnames,kind,forecastIndex = 1,colIndex=0,tmpDir="/tmp/",fnameInTar="",concatenateAxis=1,debug=False,verbosity=0)
pam.nmlSet["data_path"] = '/net/marin//mmaahn/pamtra_data/'

pam.runPamtra(35)

pyPamtraPlotter.plotTB(pam,levels=[1,400])