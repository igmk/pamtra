import pyPamtra
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
import sys

pam = pyPamtra.pyPamtra()
pam.df.addHydrometeor(('ice', -99.0, -1, 917,917 *  np.pi / 6., 3, np.pi/4., 2, 0, 100, 'exp', 3000, 3e9, -99.0, -99.0, 100e-6,  1.2e-3, 'mie-sphere', 'heymsfield10_particles',0.0))

pam = pyPamtra.importer.createUsStandardProfile(pam,hgt_lev=np.arange(1000,1300,200))
pam.p["airturb"][:] = 0.2
freqs = [35.5]#,80,150]

pam.set["verbose"] = -666
pam.set["pyVerbose"] =0
pam.nmlSet["data_path"] = "/work/mmaahn/pamtra_data/"
pam.nmlSet["randomseed"] = 100
pam.nmlSet["radar_mode"] = "spectrum"
pam.nmlSet["radar_aliasing_nyquist_interv"] = 0
pam.nmlSet["save_psd"] = True

pam.p["hydro_q"][:] = 0.002

pam.nmlSet["hydro_fullspec"] = False

pam.runPamtra(freqs,checkData=False)
#sys.exit()
    ##allocate(out_debug_radarvel(radar_nfft))
    ##allocate(out_debug_radarback(radar_nfft))
    ##allocate(out_debug_radarback_wturb(radar_nfft))
    ##allocate(out_debug_radarback_wturb_wnoise(radar_nfft))
    ##out_debug_diameter(:) = 0.d0
    ##out_debug_back_of_d(:) = 0.d0

diameter = pam.fortObject.vars_output.out_debug_diameter
sigma_D = pam.fortObject.vars_output.out_debug_back_of_d
eta_wo_turb = pam.fortObject.vars_output.out_debug_radarback
eta_turb = pam.fortObject.vars_output.out_debug_radarback_wturb
eta_turb_noisy = pam.fortObject.vars_output.out_debug_radarback_wturb_wnoise
velocity = pam.fortObject.vars_output.out_debug_radarvel


eta_wo_turb[eta_wo_turb<= 0.0] = 1e-8

matplotlib.rcParams.update({'font.size': 18})


plt.figure()
plt.plot(pam.r["psd_d"].ravel(),pam.r["psd_n"].ravel())
plt.xlabel("Diameter [m]")
plt.ylabel("number concentration [1/m^4]")
plt.title("Particle size distribution")
#plt.xscale("log")
plt.xlim(100e-6,1200e-6)
plt.savefig("radar_howto_0.png")

plt.figure()
plt.plot(diameter, np.log10(sigma_D))
plt.xlabel("Diameter [m]")
plt.ylabel("Backscattering [dB]")
plt.title("Backscattering per diameter")
#plt.xscale("log")
plt.xlim(100e-6,1200e-6)
plt.savefig("radar_howto_1.png")


plt.figure()
plt.plot(velocity, np.log10(eta_wo_turb)*10)
plt.xlabel("Velocity [m/s]]")
plt.ylabel("Backscattering [dB]")
plt.ylim(-60,0)
plt.xlim(-7.5,7.5)
plt.title("diameter to velocity")
plt.savefig("radar_howto_2.png")

plt.figure()
plt.plot(velocity,np.log10(eta_turb)*10)


plt.xlabel("Velocity [m/s]]")
plt.ylabel("Backscattering [dB]")
plt.ylim(-60,0)
plt.xlim(-7.5,7.5)
plt.title("applying turbulence")
plt.savefig("radar_howto_3.png")

plt.figure()
plt.plot(velocity, 10*np.log10(eta_turb_noisy))
plt.xlabel("Velocity [m/s]]")
plt.ylabel("Backscattering [dB]")
plt.ylim(-60,0)
plt.xlim(-7.5,7.5)
plt.title("applying noise")
plt.savefig("radar_howto_4.png")
