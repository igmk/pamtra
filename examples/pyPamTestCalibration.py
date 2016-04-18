import pyPamtra
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np


fig1 = plt.figure()



for miscalibration in [-10,-5,0,5,10]:

  pam = pyPamtra.pyPamtra()
  #!name       as_ratio    liq_ice     rho_ms    a_ms    b_ms    alpha    beta   moment_in   nbin      dist_name        p_1     p_2     p_3     p_4     d_1       d_2           scat_name   vel_size_mod           canting
  pam.df.addHydrometeor(('ice', 1.0, -1, 200,-99,-99, np.pi/4., 2, 0, 100, 'exp', 3000, 3e8, -99.0, -99.0, 100e-6,  1000e-6, 'mie-sphere', 'heymsfield10_particles',0.0))

  pam = pyPamtra.importer.createUsStandardProfile(pam,hgt_lev=np.arange(1000,1300,200))
  pam.p["airturb"][:] = 0.2
  freqs = [35.5]#,80,150]

  pam.set["verbose"] = 0
  pam.set["pyVerbose"] =0
  pam.nmlSet["randomseed"] = 10
  pam.nmlSet["radar_mode"] = "spectrum"
  pam.nmlSet["radar_aliasing_nyquist_interv"] = 3
  pam.nmlSet["hydro_adaptive_grid"] = False
  pam.nmlSet["radar_receiver_miscalibration"] = miscalibration



  pam.p["hydro_q"][:] = 0.002
  pam.nmlSet["save_psd"] = True

  pam.runPamtra(freqs,checkData=False)
  plt.plot(pam.r["radar_vel"],pam.r["radar_spectra"][0,0,0,0,0])

  print miscalibration, "Ze", pam.r["Ze"].ravel()[0], "W", pam.r["radar_moments"].ravel()[0], "SW", pam.r["radar_moments"].ravel()[1], "Sk", pam.r["radar_moments"].ravel()[2], "Ku", pam.r["radar_moments"].ravel()[3], "lS", pam.r["radar_slopes"].ravel()[0], "rS", pam.r["radar_slopes"].ravel()[1]


#same but with additional miscalibrated noise

for miscalibration in [-10,-5,0,5,10]:

  pam = pyPamtra.pyPamtra()
  #!name       as_ratio    liq_ice     rho_ms    a_ms    b_ms    alpha    beta   moment_in   nbin      dist_name        p_1     p_2     p_3     p_4     d_1       d_2           scat_name   vel_size_mod           canting
  pam.df.addHydrometeor(('ice', 1.0, -1, 200,-99,-99, np.pi/4., 2, 0, 100, 'exp', 3000, 3e8, -99.0, -99.0, 100e-6,  1000e-6, 'mie-sphere', 'heymsfield10_particles',0.0))

  pam = pyPamtra.importer.createUsStandardProfile(pam,hgt_lev=np.arange(1000,1300,200))
  pam.p["airturb"][:] = 0.2
  freqs = [35.5]#,80,150]

  pam.set["verbose"] = 0
  pam.set["pyVerbose"] =0
  pam.nmlSet["randomseed"] = 10
  pam.nmlSet["radar_mode"] = "spectrum"
  pam.nmlSet["radar_aliasing_nyquist_interv"] = 3
  pam.nmlSet["hydro_adaptive_grid"] = False
  pam.nmlSet["radar_receiver_miscalibration"] = miscalibration
  pam.nmlSet["radar_pnoise0"] -= miscalibration



  pam.p["hydro_q"][:] = 0.002
  pam.nmlSet["save_psd"] = True

  pam.runPamtra(freqs,checkData=False)
  plt.plot(pam.r["radar_vel"],pam.r["radar_spectra"][0,0,0,0,0])

  print miscalibration, "Ze", pam.r["Ze"].ravel()[0], "W", pam.r["radar_moments"].ravel()[0], "SW", pam.r["radar_moments"].ravel()[1], "Sk", pam.r["radar_moments"].ravel()[2], "Ku", pam.r["radar_moments"].ravel()[3], "lS", pam.r["radar_slopes"].ravel()[0], "rS", pam.r["radar_slopes"].ravel()[1]


