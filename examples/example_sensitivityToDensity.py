import pyPamtra
import numpy as np

pam = pyPamtra.pyPamtra()
pam.df.addHydrometeor(("iwc_q", 0.6, -1, 917., -99., -99., 0.684, 2., 3, 1, "mono_cosmo_ice", -99., -99., -99., -99., -99., -99., "mie-sphere", "heymsfield10_particles",0.0))
pam = pyPamtra.importer.createUsStandardProfile(pam,hgt_lev=[[1000,1100]]*10)
pam.nmlSet['data_path'] =  '/work/mmaahn/pamtra_data/'

pam.df.data4D["rho_ms"] = np.linspace(100,900,10).reshape(pam._shape4D)

pam.p["hydro_q"][:] = 1e-4
pam.runPamtra(35.5)

plt.figure()
plt.plot(pam.df.data4D["rho_ms"].ravel(),pam.r["Ze"].ravel())
plt.xlabel("rho_ms [kg/m^3]")
plt.ylabel("Ze [dBz]")