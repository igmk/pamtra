import pyPamtra
import pyPamtraImport

pam = pyPamtraImport.createUsStandardProfile(hgt_lev=[100,200,300])

pam.df.addHydrometeor(("ice", -99., -1, 917., 130., 3.0, 0.684, 2., 3, 1, "mono_cosmo_ice", -99., -99., -99., -99., -99., -99., "mie-sphere", "heymsfield10_particles",0.0))

pam.p["hydro_q"][0,0,1] = 1e-2

pam.nmlSet["data_path"] = '/work/mmaahn/pamtra_data/'
pam.nmlSet["radar_mode"] = "spectrum"
pam.nmlSet["radar_noise_distance_factor"] = 2.0
pam.runPamtra(35)

plt.plot(pam.r["radar_vel"],pam.r["radar_spectra"][0,0,1,0])
print pam.r["Ze"][0,0,1,0], pam.r["radar_moments"][0,0,1,0], pam.r["radar_slopes"][0,0,1,0]

