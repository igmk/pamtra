import pyPamtra
import numpy as np


reload(pyPamtra)

pam = pyPamtra.pyPamtra()

pam.df.addHydrometeor(("ice", -99., -1, 917., 130., 3.0, 0.684, 2., 3, 1, "mono_cosmo_ice", -99., -99., -99., -99., -99., -99., "mie-sphere", "heymsfield10_particles",0.0))

pam = pyPamtra.importer.createUsStandardProfile(pam,hgt_lev=[100,200,300])

pam.set["verbose"] = 3
#pam.set["pyVerbose"] = 3
#pam.nmlSet['ground_type'] = 'F'
pam.p['sfc_type'][:] = 0.
pam.p['groundtemp'][:] = 277.
#pam.p['wind10u'][:] = np.sqrt(25.)
#pam.p['wind10v'][:] = np.sqrt(24.)

pam.runPamtra(35)

import matplotlib.pyplot as plt
plt.figure()
pyPamtra.plot.plotTB(pam,levels=[1,400])
