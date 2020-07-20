import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec
import matplotlib.ticker

from rosen98_gasabs import gasabs_module

freq = np.arange(10, 600.1, 0.2) # GHz
tempK = 290
pres = 900 # milibar or hPa
vapden = 10 # G/M^3   WATER VAPOR DENSITY


fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(1, 1, 1)

o2abs = np.vectorize(gasabs_module.o2abs)(tempK, pres, vapden, freq)
abh2o = np.vectorize(gasabs_module.abh2o)(tempK, pres, vapden, freq)
absn2 = np.vectorize(gasabs_module.absn2)(tempK, pres, freq)

ax.plot(freq, abh2o + o2abs + absn2, zorder=2,
    color='black', label='Total (T = %.0f K, p = %.0f hPa)' % (tempK, pres))

ax.plot(freq, abh2o, zorder=1, linestyle='--',
    color='tab:blue', label=r'Water vapor ($\rho_{water}$ = %.1f g m$^{-3}$)' % vapden,)

ax.plot(freq, o2abs, zorder=1, linestyle=':',
    color='tab:red', label='Oxygen')

ax.plot(freq, absn2, zorder=1, linestyle='-.',
    color='tab:green', label='Nitrogen',)


ax.set_xlim(freq.min(), freq.max())
ax.set_yscale('log')

ax.legend(loc='lower left', ncol=3)

ax.set_ylabel('Absorption coefficient (km$^{-1}$)')
ax.set_xlabel('Frequency (GHz)')
fig.subplots_adjust(bottom=0.2)

plt.show()
