import pyPamtra
import os
import netCDF4

reload(pyPamtra)





plt.figure(1)
plt.clf()

plt.figure(2)
plt.clf()

pam = pyPamtra.pyPamtra()
pam.df.readFile("descriptor_file_COSMO.txt")
pam.readPamtraProfile("profile/example_input.lev")

print "##########################"

pam.runPamtra(35.5)



Ze = np.ma.masked_equal(pam.r["Ze"],-9999).ravel()

plt.figure(1)
plt.title("tb")
plt.plot(pam.r["tb"].ravel(),label="pyPamtra level")
plt.figure(2)
plt.title("Ze")
plt.plot(Ze.compressed(),label="pyPamtra level")

pam = pyPamtra.pyPamtra()
pam.df.readFile("descriptor_file_COSMO.txt")
pam.readPamtraProfile("profile/example_input.lay")
print "##########################"

pam.runPamtra(35.5)
Ze = np.ma.masked_equal(pam.r["Ze"],-9999).ravel()

plt.figure(1)
plt.plot(pam.r["tb"].ravel(),label="pyPamtra layer")
plt.figure(2)
plt.plot(Ze.compressed(),label="pyPamtra layer")

print "##########################"
os.remove("output/example_input_035.5000.nc")


os.system("./pamtra -f 35.5 -d descriptor_file_COSMO.txt -p example_input.lev")

ncData = netCDF4.Dataset("output/example_input_035.5000.nc")
Ze = np.ma.masked_equal(ncData.variables["Ze"],-9999).ravel()

plt.figure(1)
plt.title("tb")
plt.plot(ncData.variables["tb"][:].ravel(),":",lw=4,label="Pamtra level")
plt.figure(2)
plt.title("Ze")
plt.plot(Ze.compressed(),":",lw=4,label="Pamtra level")

print Ze.compressed()
ncData.close()

print "##########################"

os.system("./pamtra -f 35.5 -d descriptor_file_COSMO.txt -p example_input.lay")
ncData = netCDF4.Dataset("output/example_input_035.5000.nc")
Ze = np.ma.masked_equal(ncData.variables["Ze"],-9999).ravel()

plt.figure(1)
plt.title("tb")
plt.plot(ncData.variables["tb"][:].ravel(),":",lw=4,label="Pamtra layer")
plt.legend()

plt.figure(2)
plt.title("Ze")
plt.plot(Ze.compressed(),":",lw=4,label="Pamtra layer")
print Ze.compressed()
plt.legend()
ncData.close()
