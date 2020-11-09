from __future__ import print_function

import pyPamtra
import numpy as np
import matplotlib.pyplot as plt


pam = pyPamtra.pyPamtra()
#pam.df.addHydrometeor(('ice', 1.0, -1, 917,917 *  np.pi / 6., 3, np.pi/4., 2, 0, 100, 'exp', 3000, 3e8, -99.0, -99.0, 100e-6,  1000e-6, 'mie-sphere', 'heymsfield10_particles', 90.0))
pam.df.addHydrometeor(('ice', 0.3, -1, 917,917 *  np.pi / 6., 3, np.pi/4., 2, 0, 100, 'exp', 3000, 3e8, -99.0, -99.0, 100e-6,  1000e-6, 'tmatrix', 'heymsfield10_particles', 90.0))
#pam.df.addHydrometeor(('ice', 1.0, -1, 917,917 *  np.pi / 6., 3, np.pi/4., 2, 0, 100, 'exp', 3000, 3e8, -99.0, -99.0, 100e-6,  1000e-6, 'ssrg-rt3', 'heymsfield10_particles', 90.0))

pam = pyPamtra.importer.createUsStandardProfile(pam,hgt_lev=np.arange(1000,1300,200))

freqs = [9.5]

pam.set["verbose"] = 0
pam.set["pyVerbose"] = 0

pam.nmlSet["passive"] = False
pam.nmlSet["radar_mode"] = "spectrum"
#pam.nmlSet["radar_mode"] = "simple"
pam.nmlSet["randomseed"] = 10
#pam.nmlSet["radar_dualpol"] = "both"
pam.nmlSet["radar_polarisation"] = "NN,HH,VV,HV,VH"

pam.p["hydro_q"][:] = 0.001

pam.runPamtra(freqs,checkData=False)

#plot(pam.r["Ze"].ravel())
#plt.plot(pam.r["radar_vel"],pam.r["radar_spectra"][0,0,0,2])

for i_p, pol in enumerate(pam.set['radar_pol']):
    print(pol, pam.r["Ze"][...,i_p,0])
print('ave HH VV ', np.log10(0.5*(10**pam.r["Ze"][...,1,0]+10**pam.r["Ze"][...,2,0])))
plt.figure(11)
plt.clf()
for i_p in range(pam.set["radar_npol"]):
  plt.plot(pam.r["radar_vel"][0,:],pam.r["radar_spectra"][0,0,0,0,i_p,:],label=pam.set["radar_pol"][i_p] + " " + str(pam.r["Ze"][...,i_p,0]))
  
plt.legend(loc="upper left")  
plt.title(pam.df.data["canting"])

plt.show()

#print "ZDR, LDR", pam.r["zdr"], pam.r["ldr"]
#print "Ze", pam.r["Ze"]

##print pam.r["Ze"][0,0,1,0], pam.r["radar_moments"][0,0,1,0], pam.r["radar_slopes"][0,0,1,0]