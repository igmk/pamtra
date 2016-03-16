import pyPamtra
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np

verbosity = 0
plt.figure()

pam = pyPamtra.pyPamtra()
#!name       as_ratio    liq_ice     rho_ms    a_ms    b_ms    alpha    beta   moment_in   nbin      dist_name        p_1     p_2     p_3     p_4     d_1       d_2           scat_name   vel_size_mod           canting
pam.df.addHydrometeor(('ice', 1.0, -1, 200,-99,-99, np.pi/4., 2, 0, 100, 'exp', 3000, 3e8, -99.0, -99.0, 100e-6,  1000e-6, 'mie-sphere', 'heymsfield10_particles',0.0))

pam = pyPamtra.importer.createUsStandardProfile(pam,hgt_lev=np.arange(1000,1300,200))
pam.p["airturb"][:] = 0.2
freqs = [35.5]#,80,150]

pam.set["verbose"] = 0
pam.set["pyVerbose"] =0
pam.nmlSet["passive"] = False
pam.nmlSet["randomseed"] = 0
pam.nmlSet["radar_mode"] = "spectrum"
pam.nmlSet["radar_aliasing_nyquist_interv"] = 3
pam.nmlSet["hydro_adaptive_grid"] = False
pam.nmlSet["conserve_mass_rescale_dsd"] = False


pam.p["hydro_q"][:] = 0.002
pam.nmlSet["save_psd"] = True

pam.set['verbose'] = verbosity
pam.set['pyVerbose'] = verbosity

# pam.writeNmlFile("test.nml")
# pam.writePamtraProfile("test.lev")

pam.runPamtra(freqs,checkData=False)
plt.plot(pam.r["radar_vel"],pam.r["radar_spectra"][0,0,0,0,0])

# sys.exit()

pamFS = pyPamtra.pyPamtra()
#!name       as_ratio    liq_ice     rho_ms    a_ms    b_ms    alpha    beta   moment_in   nbin      dist_name        p_1     p_2     p_3     p_4     d_1       d_2           scat_name   vel_size_mod           canting
pamFS.df.addHydrometeor(('ice', 1.0, -1, 200,-99,-99, np.pi/4., 2, 0, 100, 'exp', 3000, 3e8, -99.0, -99.0, 100e-6,  1000e-6, 'mie-sphere', 'heymsfield10_particles',0.0))

pamFS = pyPamtra.importer.createUsStandardProfile(pamFS,hgt_lev=np.arange(1000,1300,200))
pamFS.p["airturb"][:] = 0.2
pamFS.set["verbose"] = 5
pamFS.set["pyVerbose"] =0
pamFS.nmlSet["passive"] = False
pamFS.nmlSet["randomseed"] = 0
pamFS.nmlSet["hydro_adaptive_grid"] = False
pamFS.nmlSet["conserve_mass_rescale_dsd"] = False
pamFS.nmlSet["radar_mode"] = "spectrum"
pamFS.nmlSet["radar_use_hildebrand"] = True


pamFS.nmlSet["radar_aliasing_nyquist_interv"] = 3
pamFS.p["hydro_q"][:] = 0.002

pamFS.nmlSet["randomseed"] = 0
pamFS.nmlSet["save_psd"] = True

pamFS.nmlSet["hydro_fullspec"] = True
pamFS.df.addFullSpectra()

pamFS.df.dataFullSpec["d_bound_ds"][0,0,0,0,:] = np.linspace(100e-6,1000e-6,101)
pamFS.df.dataFullSpec["d_ds"][0,0,0,0,:] = pamFS.df.dataFullSpec["d_bound_ds"][0,0,0,0,:-1] + 0.5 * np.diff(pamFS.df.dataFullSpec["d_bound_ds"][0,0,0,0,:])
pamFS.df.dataFullSpec["rho_ds"][0,0,0,0,:] = 200
pamFS.df.dataFullSpec["n_ds"][0,0,0,0,:] = 3e8 * np.exp(-3000 * pamFS.df.dataFullSpec["d_ds"][0,0,0,0,:]) *np.diff(pamFS.df.dataFullSpec["d_bound_ds"][0,0,0,0,:])
pamFS.df.dataFullSpec["area_ds"][0,0,0,0,:] = np.pi/4. *  pamFS.df.dataFullSpec["d_ds"][0,0,0,0,:] ** 2
pamFS.df.dataFullSpec["mass_ds"][0,0,0,0,:] = np.pi/6. *pamFS.df.dataFullSpec["rho_ds"][0,0,0,0,:] *  pamFS.df.dataFullSpec["d_ds"][0,0,0,0,:] ** 3
pamFS.df.dataFullSpec["as_ratio"][0,0,0,0,:] = 1.0


pamFS.set['verbose'] = verbosity
pamFS.set['pyVerbose'] = verbosity

pamFS.runPamtra(freqs,checkData=False)

#plt.plot(pam.r["radar_vel"],pam.r["radar_spectra"][0,0,0,0])
plt.plot(pamFS.r["radar_vel"],pamFS.r["radar_spectra"][0,0,0,0,0])


print pamFS.r["Ze"], pam.r["Ze"]
