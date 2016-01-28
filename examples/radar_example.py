import pyPamtra
import numpy as np
import matplotlib.pyplot as plt

pam = pyPamtra.pyPamtra()
pam.df.addHydrometeor(('ice', -99.0, -1, 917,917 *  np.pi / 6., 3, np.pi/4., 2, 0, 10, 'exp', 3000, 3e8, -99.0, -99.0, 100e-6,  200e-6, 'mie-sphere', 'heymsfield10_particles',0.0))
#pam.df.addHydrometeor(("ice", -99., -1, 917., 130., 3.0, 0.684, 2., 3, 1, "mono_cosmo_ice", -99., -99., -99., -99., -99., -99., "mie-sphere", "heymsfield10_particles",0.0))

pam = pyPamtra.importer.createUsStandardProfile(pam,hgt_lev=[100,200,300])

pam.set["verbose"] = 0
pam.p["hydro_q"][0,0,1] = 1e-2
pam.p["airturb"][:] = 1e-2
pam.nmlSet["data_path"] = '/net/marin//mmaahn/pamtra_data/'
pam.nmlSet["radar_mode"] = "spectrum"
pam.nmlSet["radar_noise_distance_factor"] = 2.0
pam.nmlSet["radar_save_noise_corrected_spectra"]=  False
pam.nmlSet["randomseed"]=  10


pam.runPamtra(35)
plt.figure()
plt.plot(pam.r["radar_vel"],pam.r["radar_spectra"][0,0,1,0,0])
print pam.r["Ze"][0,0,1,0,0], pam.r["radar_moments"][0,0,1,0,0], pam.r["radar_slopes"][0,0,1,0,0]
plt.show()
#pam.writeResultsToNetCDF("/tmp/test.nc")

