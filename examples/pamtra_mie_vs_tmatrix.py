import pyPamtra
import os
import netCDF4

reload(pyPamtra)

descriptorFile = "../descriptorfiles/descriptor_file_COSMO_1mom.txt"


plt.figure(1)
plt.clf()

plt.figure(2)
plt.clf()

pam = pyPamtra.pyPamtra()
pam.df.readFile(descriptorFile)
pam.readPamtraProfile("../profile/example_input.lev")

print "##########################"

#pam.set["pyVerbose"] = 10
pam.runPamtra(35.5)

Ze = np.ma.masked_equal(pam.r["Ze"],-9999).ravel()
plt.figure(1)
plt.title("tb")
plt.plot(pam.r["tb"].ravel(),label="Mie")
plt.figure(2)
plt.title("Ze")
plt.plot(Ze.compressed(),label="Mie")



pam = pyPamtra.pyPamtra()
pam.df.readFile(descriptorFile)
pam.readPamtraProfile("../profile/example_input.lev")

pam.df.data["scat_name"][:] = "tmatrix"
pam.df.data["as_ratio"][:] = 1.0

print "##########################"

#pam.set["pyVerbose"] = 10
pam.runPamtra(35.5)

Ze = np.ma.masked_equal(pam.r["Ze"],-9999).ravel()
plt.figure(1)
plt.title("tb")
plt.plot(pam.r["tb"].ravel(),label="tmatrix 1.0")
plt.figure(2)
plt.title("Ze")
plt.plot(Ze.compressed(),label="tmatrix 1.0")


pam = pyPamtra.pyPamtra()
pam.df.readFile(descriptorFile)
pam.readPamtraProfile("../profile/example_input.lev")

pam.df.data["scat_name"][:] = "tmatrix"
pam.df.data["as_ratio"][:] = 0.6

print "##########################"

#pam.set["pyVerbose"] = 10
pam.runPamtra(35.5)

Ze = np.ma.masked_equal(pam.r["Ze"],-9999).ravel()
plt.figure(1)
plt.title("tb")
plt.plot(pam.r["tb"].ravel(),label="tmatrix 0.6")
plt.figure(2)
plt.title("Ze")
plt.plot(Ze.compressed(),label="tmatrix 0.6")

plt.figure(1)
plt.legend()

plt.figure(2)
plt.legend()