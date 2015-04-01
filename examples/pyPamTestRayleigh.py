import pyPamtra
from copy import deepcopy

pam = pyPamtra.pyPamtra()
#!name       as_ratio    liq_ice     rho_ms    a_ms    b_ms    alpha    beta   moment_in   nbin      dist_name        p_1     p_2     p_3     p_4     d_1       d_2           scat_name   vel_size_mod           canting
pam.df.addHydrometeor(('rain1', -99., 1,  -99., -99., -99.,  -99.,  -99., 1, 1, 'mono', -99.0, -99.0, -99.0, -99.0, 0.001,  -99.0, 'mie-sphere', 'heymsfield10_particles',0.0))
pam.df.addHydrometeor(('rain2', -99., 1,  -99., -99., -99.,  -99.,  -99., 1, 1, 'mono', -99.0, -99.0, -99.0, -99.0, 0.0001,  -99.0, 'mie-sphere', 'heymsfield10_particles',0.0))
pam.df.addHydrometeor(('ice1', -99., -1,  -99., 0.0121, 1.9,  -99.,  -99., 1, 1, 'mono', -99.0, -99.0, -99.0, -99.0, 0.001,  -99.0, 'mie-sphere', 'heymsfield10_particles',0.0))
pam.df.addHydrometeor(('ice2', -99., -1,  -99., 0.0121, 1.9,  -99.,  -99., 1, 1, 'mono', -99.0, -99.0, -99.0, -99.0, 0.0001,  -99.0, 'mie-sphere', 'heymsfield10_particles',0.0))

pam = pyPamtra.importer.createUsStandardProfile(pam,hgt_lev=np.array([np.arange(1000,1300,200)]*4))

pam.p["hydro_n"][0,0,0,0] = 1
pam.p["hydro_n"][1,0,0,1] = 1e6
pam.p["hydro_n"][2,0,0,2] = 1
pam.p["hydro_n"][3,0,0,3] = 1e4
freqs = [10]#,80,150]

pam.set["verbose"] = 0
pam.set["pyVerbose"] =0
pam.nmlSet["data_path"] = "/work/mmaahn/pamtra_data/"


pam.runPamtra(freqs,checkData=False)
plt.figure()
plt.plot(pam.r["Ze"].ravel())
