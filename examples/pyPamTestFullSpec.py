import pyPamtra
import pyPamtraLibWrapper
import pyPamtraImport

pam = pyPamtraImport.createUsStandardProfile(hgt_lev=np.arange(1000,1300,200))

freqs = [35.5,80,150]

pam.set["verbose"] = 5
pam.set["pyVerbose"] =5

pam.nmlSet["data_path"] = "/work/mmaahn/pamtra_data/"
#pam.nmlSet["save_psd"] = True
pam.nmlSet["radar_mode"] = "spectrum"


pam.df.addHydrometeor(('ice', -99.0, -1, 917,917 *  pi / 6., 3, pi/4., 2, 0, 10, 'exp', 3000, 3e8, -99.0, -99.0, 100e-6,  1000e-6, 'mie-sphere', 'test'))
pam.p["hydro_q"][:] = 0.001

pam.nmlSet["hydro_fullspec"] = True
pam.df.addFullSpectra()

pam.df.dataFullSpec["d_bound_ds"][0,0,0,0,:] = np.linspace(100e-6,1000e-6,10)
pam.df.dataFullSpec["delta_d_ds"][0,0,0,0,:] = np.diff(pam.df.dataFullSpec["d_bound_ds"][0,0,0,0,:])
pam.df.dataFullSpec["density2scat"][0,0,0,0,:] = 917
pam.df.dataFullSpec["diameter2scat"][0,0,0,0,:] = pam.df.dataFullSpec["d_bound_ds"][0,0,0,0,:]
pam.df.dataFullSpec["f_ds"][0,0,0,0,:] = 3e8 * np.exp(-3000 * pam.df.dataFullSpec["d_bound_ds"][0,0,0,0,:])
pam.df.dataFullSpec["area_ds"][0,0,0,0,:] = pi/4. *  pam.df.dataFullSpec["d_bound_ds"][0,0,0,0,:] ** 2
pam.df.dataFullSpec["mass_ds"][0,0,0,0,:] =pi / 6. *917 *  pam.df.dataFullSpec["d_bound_ds"][0,0,0,0,:] ** 3
pam.df.dataFullSpec["as_ratio"][0,0,0,0,:] = 0.6


#pam.df.data4D["p_1"] = np.ones((5,1,5,5)) * -99

pam.runPamtra(freqs,checkData=False)

#plot(pam.r["Ze"].ravel())
plt.plot(pam.r["radar_vel"],pam.r["radar_spectra"][0,0,0,2])

#results, pyPamtraLib = PamtraFortranWrapper(pam.set,pam.nmlSet,pam.df,1)


#print pam.r["Ze"][0,0,-1,0]
#print pam.r["tb"][0,0,-1,1,-1,-1]

#print "ref", 21.6671557514, 45.4300752495




