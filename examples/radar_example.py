from __future__ import print_function

import pyPamtra
import numpy as np
import matplotlib.pyplot as plt
import sys


plt.figure()
for turb in [1,0.5,0.1,0.01,0]:
  Ws = []

  pam = pyPamtra.pyPamtra()
  # pam.df.addHydrometeor(('ice', -99.0, -1, 917,917 *  np.pi / 6., 3, np.pi/4., 2, 0, 10, 'exp', 3000, 3e8, -99.0, -99.0, 100e-6,  200e-6, 'mie-sphere', 'heymsfield10_particles',0.0))
  # pam.df.addHydrometeor(("ice", -99., -1, 917., 130., 3.0, 0.684, 2., 3, 1, "mono_cosmo_ice", -99., -99., -99., -99., -99., -99., "mie-sphere", "heymsfield10_particles",0.0))
  # pam.df.addHydrometeor(('ice', 0.5, -1, 917,917 *  np.pi / 6., 3, np.pi/4., 2, 0, 10, 'exp', 3000, 3e8, -99.0, -99.0, 100e-6,  1000e-6, 'mie-sphere', 'heymsfield10_particles', 90.0))
  pam.df.addHydrometeor(('ice', -99.0, -1, 200,-99,-99, np.pi/4., 2, 0, 100, 'exp', 4000, 3e6, -99.0, -99.0, 100e-6,  1.2e-3, 'mie-sphere', 'heymsfield10_particles',0.0))

  pam = pyPamtra.importer.createUsStandardProfile(pam,hgt_lev=[100,200,300])

  pam.set["verbose"] = 0
  pam.p["hydro_q"][0,0,1] = 1e-2
  pam.p["airturb"][:] = turb
  pam.p["wind_w"][:] =0.0
  pam.nmlSet["radar_mode"] = "spectrum"
  pam.nmlSet["radar_noise_distance_factor"] = 2.0
  pam.nmlSet["radar_save_noise_corrected_spectra"]=  False
  pam.nmlSet["randomseed"]=  10
  pam.nmlSet['radar_airmotion'] = True
  pam.nmlSet['radar_airmotion_model'] = 'constant'
  pam.nmlSet['radar_airmotion_vmax'] = 0.0
  pam.nmlSet['radar_airmotion_vmin'] = 0.0
  pam.nmlSet['radar_aliasing_nyquist_interv'] = 1
  pam.nmlSet['radar_no_Ave'] = 5

  pam.runPamtra(35)
  plt.plot(pam.r["radar_vel"][0],pam.r["radar_spectra"][0,0,1,0,0],label=turb)
  print(turb, pam.r["Ze"][0,0,1,0,0], pam.r["radar_moments"][0,0,1,0,0], pam.r["radar_slopes"][0,0,1,0,0])

plt.legend()
plt.xlabel('Doppler velocity')
plt.ylabel('spectral reflectivity')
plt.show()
plt.legend()
# print pam.fortObject.vars_output.out_debug_radarback_wturb

