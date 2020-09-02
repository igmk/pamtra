from __future__ import print_function

import pyPamtra
import shutil
import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import imp

imp.reload(pyPamtra)

descriptorFile = "../descriptorfiles/descriptor_file_COSMO_1mom.txt"


plt.figure(1)
plt.clf()

plt.figure(2)
plt.clf()

pam = pyPamtra.pyPamtra()
pam.df.readFile(descriptorFile)
pam.readPamtraProfile("../profile/example_input.lev")
pam.filterProfiles(np.array([[True, False], [False, False]]))
pam.nmlSet['tmatrix_db'] = 'file'
pam.nmlSet['tmatrix_db_path'] = 'example_db/'

#pam.set["pyVerbose"] = 10
pam.runPamtra(35.5)

Ze = np.ma.masked_equal(pam.r["Ze"], -9999).ravel()
plt.figure(1)
plt.title("tb")
plt.plot(pam.r["tb"].ravel(), label="Mie")
plt.figure(2)
plt.title("Ze")
plt.plot(Ze.compressed(), label="Mie")

print("##########################")


pam = pyPamtra.pyPamtra()
pam.df.readFile(descriptorFile)
pam.readPamtraProfile("../profile/example_input.lev")
pam.filterProfiles(np.array([[True, False], [False, False]]))
pam.nmlSet['tmatrix_db'] = 'file'
pam.nmlSet['tmatrix_db_path'] = 'example_db/'


pam.df.data["scat_name"][:] = "tmatrix"
pam.df.data["as_ratio"][:] = 1.0

#pam.set["pyVerbose"] = 10
pam.runPamtra(35.5)

Ze = np.ma.masked_equal(pam.r["Ze"], -9999).ravel()
plt.figure(1)
plt.title("tb")
plt.plot(pam.r["tb"].ravel(), label="tmatrix 1.0")
plt.figure(2)
plt.title("Ze")
plt.plot(Ze.compressed(), label="tmatrix 1.0")


print("##########################")

pam = pyPamtra.pyPamtra()
pam.df.readFile(descriptorFile)
pam.readPamtraProfile("../profile/example_input.lev")
pam.filterProfiles(np.array([[True, False], [False, False]]))
pam.nmlSet['tmatrix_db'] = 'file'
pam.nmlSet['tmatrix_db_path'] = 'example_db/'


pam.df.data["scat_name"][:] = "tmatrix"
pam.df.data["as_ratio"][:] = 1.0

#pam.set["pyVerbose"] = 10
pam.runPamtra(35.5)

Ze = np.ma.masked_equal(pam.r["Ze"], -9999).ravel()
plt.figure(1)
plt.title("tb")
plt.plot(pam.r["tb"].ravel(), label="tmatrix 1.0 from database", marker='x')
plt.figure(2)
plt.title("Ze")
plt.plot(Ze.compressed(), label="tmatrix 1.0 from database", marker='x')

print("##########################")


pam = pyPamtra.pyPamtra()
pam.df.readFile(descriptorFile)
pam.readPamtraProfile("../profile/example_input.lev")
pam.filterProfiles(np.array([[True, False], [False, False]]))
pam.nmlSet['tmatrix_db'] = 'file'
pam.nmlSet['tmatrix_db_path'] = 'example_db/'


pam.df.data["scat_name"][:] = "tmatrix"
pam.df.data["as_ratio"][:] = 0.6


#pam.set["pyVerbose"] = 10
pam.runPamtra(35.5)

Ze = np.ma.masked_equal(pam.r["Ze"], -9999).ravel()
plt.figure(1)
plt.title("tb")
plt.plot(pam.r["tb"].ravel(), label="tmatrix 0.6")
plt.figure(2)
plt.title("Ze")
plt.plot(Ze.compressed(), label="tmatrix 0.6")

print("##########################")

pam = pyPamtra.pyPamtra()
pam.df.readFile(descriptorFile)
pam.readPamtraProfile("../profile/example_input.lev")
pam.filterProfiles(np.array([[True, False], [False, False]]))
pam.nmlSet['tmatrix_db'] = 'file'
pam.nmlSet['tmatrix_db_path'] = 'example_db/'


pam.df.data["scat_name"][:] = "tmatrix"
pam.df.data["as_ratio"][:] = 0.6


#pam.set["pyVerbose"] = 10
pam.runPamtra(35.5)

Ze = np.ma.masked_equal(pam.r["Ze"], -9999).ravel()
plt.figure(1)
plt.title("tb")
plt.plot(pam.r["tb"].ravel(), label="tmatrix 0.6 from database", marker='+')
plt.figure(2)
plt.title("Ze")
plt.plot(Ze.compressed(), label="tmatrix 0.6 from database", marker='+')

print("##########################")


plt.figure(1)
plt.legend()

plt.figure(2)
plt.legend()

print('cleaning up')
shutil.rmtree('example_db')
print('cleaning up done')
