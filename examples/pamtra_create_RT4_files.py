"""
Example how to write a files to run the passive model with the pure RT4 code
by Frank Evans. runFile contains filenames of the files in path
"""
import numpy as np
import pyPamtra

pam = pyPamtra.pyPamtra()

pamData = {}
pamData["hgt_lev"] = [[1000,1100]]
pamData["temp_lev"] = [[263.15,263.15]]
pamData["press_lev"] = [[80001, 80000]]
pamData["relhum_lev"] = [[0, 0]]
pamData["hydro_q"] = [[1e-5]]
pam.df.addHydrometeor(('ice', 0.6, -1,  917.,  -99.0, -99.0,  -99.,  -99., 3, 1, 'mono', -99.0, -99.0, -99.0, -99.0, 0.005,  -99., "tmatrix", 'heymsfield10_particles',20.0))
pam.createProfile(**pamData)

pam.nmlSet["save_ssp"] = True
pam.nmlSet["save_psd"] = True

pam.runPamtra(89.0)
pam.dumpForRT4("runFile.txt","scatfiles")