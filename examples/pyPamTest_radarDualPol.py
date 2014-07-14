import pyPamtra
import pyPamtraLibWrapper
import pyPamtraImport

pam = pyPamtra.pyPamtra()
pam.df.addHydrometeor(('ice', 0.6, -1, 917,917 *  pi / 6., 3, pi/4., 2, 0, 10, 'exp', 3000, 3e8, -99.0, -99.0, 100e-6,  1000e-6, 'tmatrix', 'heymsfield10_particles', 0.0))

pam = pyPamtraImport.createUsStandardProfile(pam,hgt_lev=np.arange(1000,1300,200))

freqs = [9.5]

pam.set["verbose"] = 0
pam.set["pyVerbose"] =0

pam.nmlSet["data_path"] = "/work/mmaahn/pamtra_data/"
pam.nmlSet["passive"] = False
pam.nmlSet["radar_mode"] = "spectrum"
pam.nmlSet["randomseed"] = 10
#pam.nmlSet["radar_dualpol"] = "both"
pam.nmlSet["radar_polarisation"] = "NN,NN,NN,NN,NN"

pam.p["hydro_q"][:] = 0.001

pam.runPamtra(freqs,checkData=False)

#plot(pam.r["Ze"].ravel())
#plt.plot(pam.r["radar_vel"],pam.r["radar_spectra"][0,0,0,2])

plt.figure()
for i_p in range(pam.set["radar_npol"]):
  plt.plot(pam.r["radar_vel"],pam.r["radar_spectra"][0,0,0,0,i_p])


#print "ZDR, LDR", pam.r["zdr"], pam.r["ldr"]
#print "Ze", pam.r["Ze"]

##print pam.r["Ze"][0,0,1,0], pam.r["radar_moments"][0,0,1,0], pam.r["radar_slopes"][0,0,1,0]

