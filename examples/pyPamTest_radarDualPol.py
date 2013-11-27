import pyPamtra
import pyPamtraLibWrapper
import pyPamtraImport

pam = pyPamtraImport.createUsStandardProfile(hgt_lev=np.arange(1000,1300,200))

freqs = [9.5]

pam.set["verbose"] = 0
pam.set["pyVerbose"] =0

pam.nmlSet["data_path"] = "/work/mmaahn/pamtra_data/"
pam.nmlSet["passive"] = False
pam.nmlSet["radar_mode"] = "simple"
pam.nmlSet["radar_dualpol"] = "both"

pam.df.addHydrometeor(('ice', 0.6, -1, 917,917 *  pi / 6., 3, pi/4., 2, 0, 10, 'exp', 3000, 3e8, -99.0, -99.0, 100e-6,  1000e-6, 'tmatrix', 'test'))
pam.p["hydro_q"][:] = 0.001

pam.runPamtra(freqs,checkData=False)

#plot(pam.r["Ze"].ravel())
#plt.plot(pam.r["radar_vel"],pam.r["radar_spectra"][0,0,0,2])

print "ZDR, LDR", pam.r["zdr"], pam.r["ldr"]
print "Ze", pam.r["Ze"]

#print pam.r["Ze"][0,0,1,0], pam.r["radar_moments"][0,0,1,0], pam.r["radar_slopes"][0,0,1,0]

