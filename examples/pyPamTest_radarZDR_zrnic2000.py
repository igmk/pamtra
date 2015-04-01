import pyPamtra
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt

freqs = [3e8/0.0545/1e9]
nD = 31
temps = [0,10,20,30]
#temps = [0,300 + pyPamtra.meteoSI.Tnull,301 + pyPamtra.meteoSI.Tnull]
styles = ["k:","k-.","k-","k--","r"]

nT = len(temps)
pam = pyPamtra.pyPamtra()
#!name       as_ratio    liq_ice     rho_ms    a_ms    b_ms    alpha    beta   moment_in   nbin      dist_name        p_1     p_2     p_3     p_4     d_1       d_2           scat_name   vel_size_mod           canting
#pam.df.addHydrometeor(('ice', 1.0, -1, 917,917 *  pi / 6., 3, pi/4., 2, 0, 10, 'exp', 3000, 3e8, -99.0, -99.0, 100e-6,  1000e-6, 'tmatrix', 'heymsfield10_particles',0.0))
pam.df.addHydrometeor(("rain", -99., 1, -99., -99., -99., -99., -99., 0, 1, "fullSpec", -99., -99., -99., -99., -99., -99., "tmatrix", "khvorostyanov01_drops",90.0))
#pam.df.addHydrometeor(("rain", -99., 1, -99., -99., -99., -99., -99., 0, 1, "fullSpec", -99., -99., -99., -99., -99., -99., "mie-sphere", "khvorostyanov01_drops",90.0))

pam = pyPamtra.importer.createUsStandardProfile(pam,hgt_lev=np.array([[np.arange(1000,1300,200).tolist()]*nT]*nD))

for tt in xrange(nT):
  pam.p["temp_lev"][:,tt,:] = temps[tt] - pyPamtra.meteoSI.Tnull

pam.p["relhum_lev"][:] = 10

pam.p["airturb"][:] = 0
pam.set["verbose"] = 0
pam.set["pyVerbose"] =0
pam.nmlSet["data_path"] = "/work/mmaahn/pamtra_data/"
pam.nmlSet["randomseed"] = 0
pam.nmlSet["radar_mode"] = "simple"
pam.nmlSet["tmatrix_db"] = "none"
#pam.nmlSet["tmatrix_db"] = "file"
#pam.nmlSet["tmatrix_db_path"] = "/ssdwork/mmaahn/tmatrix_db/"
pam.nmlSet["radar_polarisation"] = "NN,HH,VV"#,HV"
pam.nmlSet["liq_mod"] = "Ray"


pam.p["hydro_q"][:] = 0.002

pam.nmlSet["randomseed"] = 0

pam.nmlSet["hydro_fullspec"] = True
pam.df.addFullSpectra()

rho_water = 1000.

D_es = np.linspace(0.5,8,nD)*1e-3
for dd,d_e in enumerate(D_es):

  AR = 1.03 - 0.062*d_e*1e3
  #AR = 1
  DMax = d_e/(AR)**(1./3.)
  
  print d_e, DMax, AR
  
  pam.df.dataFullSpec["d_ds"][dd,:,0,0,:] = DMax
  pam.df.dataFullSpec["d_bound_ds"][dd,:,0,0,:] = [0.99 ,1.01 ] * pam.df.dataFullSpec["d_ds"][dd,:,0,0,:] 

  pam.df.dataFullSpec["rho_ds"][dd,:,0,0,:] =rho_water
  pam.df.dataFullSpec["n_ds"][dd,:,0,0,:] = 1.
  pam.df.dataFullSpec["area_ds"][dd,:,0,0,:] = np.pi/4. *  pam.df.dataFullSpec["d_ds"][dd,:,0,0,:] ** 2
  pam.df.dataFullSpec["mass_ds"][dd,:,0,0,:] =np.pi / 6. *rho_water *  pam.df.dataFullSpec["d_ds"][dd,:,0,0,:] ** 3 * AR

  pam.df.dataFullSpec["as_ratio"][dd,:,0,0,:] = AR

  pam.df.dataFullSpec["canting"][:] = 90


pam.runPamtra(freqs,checkData=False)



#plt.figure()
#plt.title("NN")
#for tt in xrange(nT):
  #plt.plot(D_es,pam.r["Ze"][:,tt,0,0,0],styles[tt],label=str(temps[tt])+"C")
#plt.plot(D_es, np.log10((D_es*1000)**6)*10,"r")
#plt.ylim(-20,70)  
#plt.ylabel("Ze [dBz]")

plt.figure()
plt.title("HH")
for tt in xrange(nT):
  plt.plot(D_es*1000,pam.r["Ze"][:,tt,0,0,1],styles[tt],label=str(temps[tt])+"C")  
plt.ylim(-20,70)  
plt.ylabel("Ze [dBz]")

plt.figure()
plt.title("ZDR")
for tt in xrange(nT):
  plt.plot(D_es*1000,pam.r["Ze"][:,tt,0,0,1]-pam.r["Ze"][:,tt,0,0,2],styles[tt],label=str(temps[tt])+"C")
plt.ylim(0,10)  
plt.ylabel("ZDR [dB]")

#plt.figure()
#plt.title("LDR")
#for tt in xrange(nT):
  #plt.plot(D_es*1000,pam.r["Ze"][:,tt,0,0,3]-pam.r["Ze"][:,tt,0,0,1],styles[tt],label=str(temps[tt])+"C")
##plt.ylim(0,10)  
#plt.ylabel("LDR [dB]")

#plt.plot(pam.r["radar_vel"],pam.r["radar_spectra"][0,0,0,0])


