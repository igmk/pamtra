import pyPamtra
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
import sys

freqs = [35.5]#,80,150]
verbosity =0

pam = pyPamtra.pyPamtra()
pam.df.addHydrometeor(('ice', -99.0, -1, 200,-99,-99, np.pi/4., 2, 0, 100, 'exp', 4000, 3e6, -99.0, -99.0, 100e-6,  1.2e-3, 'mie-sphere', 'heymsfield10_particles',0.0))

pam = pyPamtra.importer.createUsStandardProfile(pam,hgt_lev=np.arange(1000,1300,200))

EDR = 1e-0
wind_uv = 30
beamwidth_deg = 0.3
integration_time = 1.4
frequency = freqs[0]
pam.addSpectralBroadening(EDR, wind_uv, beamwidth_deg, integration_time, frequency, kolmogorov=0.5)

pam.set["verbose"] = verbosity
pam.set["pyVerbose"] =verbosity
pam.nmlSet["randomseed"] = 1
pam.nmlSet["radar_mode"] = "spectrum"
pam.nmlSet["radar_aliasing_nyquist_interv"] = 1
pam.nmlSet["save_psd"] = True

pam.p["hydro_q"][:] = 0.00001

pam.runPamtra(freqs,checkData=False)
#plt.figure()
plt.plot(pam.r['radar_vel'].ravel(),pam.r['radar_spectra'].ravel())

pam2 = pyPamtra.pyPamtra()
pam2.df.addHydrometeor(('ice', -99.0, -1, 200,-99,-99, np.pi/4., 2, 0, 100, 'exp', 4000, 3e6, -99.0, -99.0, 100e-6,  1.2e-3, 'mie-sphere', 'heymsfield10_particles',0.0))

pam2 = pyPamtra.importer.createUsStandardProfile(pam2,hgt_lev=np.arange(1000,1300,200))
pam2.p['wind_uv'][:] = wind_uv
pam2.p['turb_edr'][:] = EDR


pam2.set["verbose"] = verbosity
pam2.set["pyVerbose"] =verbosity
pam2.nmlSet["radar_fwhr_beamwidth_deg"] = beamwidth_deg
pam2.nmlSet["radar_integration_time"] = integration_time
pam2.nmlSet["randomseed"] = 1
pam2.nmlSet["radar_mode"] = "spectrum"
pam2.nmlSet["radar_aliasing_nyquist_interv"] = 1
pam2.nmlSet["save_psd"] = True

pam2.p["hydro_q"][:] = 0.00001

pam2.runPamtra(freqs,checkData=False)


plt.plot(pam2.r['radar_vel'].ravel(),pam2.r['radar_spectra'].ravel())
plt.show()