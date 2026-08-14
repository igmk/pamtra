from __future__ import print_function

import pyPamtra
import os
import shutil
import subprocess
import netCDF4
from matplotlib import pylab as plt
import numpy as np

if os.getenv('PAMTRA_DATADIR') is None:
  raise SystemError('Set $PAMTRA_DATADIR environment variable first')


def find_pamtra_binary():
  """Locate the standalone pamtra CLI binary (built via meson-python).

  A regular `pip install .` puts it on PATH; an editable dev install
  doesn't (meson-python's editable installs only cover the Python
  extension), so also check the `build/` dir made by `meson compile -C
  build` / `pixi run build-cli`.
  """
  found = shutil.which("pamtra")
  if found is not None:
    return found
  repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
  candidate = os.path.join(repo_root, "build", "pamtra")
  if os.path.isfile(candidate):
    return candidate
  raise SystemError(
      "pamtra binary not found on PATH or in ../build/pamtra -- build it "
      "with `pip install .` (puts it on PATH) or `pixi run build-cli` "
      "(-> build/pamtra)")


pamtraBinary = find_pamtra_binary()

descriptorFile = "../descriptorfiles/descriptor_file_COSMO_1mom.txt"

freqs = [35.5]
plt.figure(1)
plt.clf()

plt.figure(2)
plt.clf()

pam = pyPamtra.pyPamtra()
pam.df.readFile(descriptorFile)
pam.readPamtraProfile("../profile/example_input.lev")
pam.readNmlFile("pamtra_vs_pyPamtra.nml")

print("##########################")



pam.set["pyVerbose"] = 10
pam.set["verbose"] = 10

pam.runPamtra(freqs)



Ze = np.ma.masked_equal(pam.r["Ze"],-9999).ravel()

plt.figure(1)
plt.title("tb")
plt.plot(pam.r["tb"].ravel(),label="pyPamtra level")
plt.figure(2)
plt.title("Ze")
plt.plot(Ze.compressed(),label="pyPamtra level")
print(Ze.compressed())

pam = pyPamtra.pyPamtra()
pam.df.readFile(descriptorFile)
pam.readPamtraProfile("../profile/example_input.lay")
pam.readNmlFile("pamtra_vs_pyPamtra.nml")
print("##########################")


pam.runPamtra(freqs)
Ze = np.ma.masked_equal(pam.r["Ze"],-9999).ravel()

plt.figure(1)
plt.plot(pam.r["tb"].ravel(),label="pyPamtra layer")
plt.figure(2)
plt.plot(Ze.compressed(),label="pyPamtra layer")
print(Ze.compressed())

print("##########################")
try:os.remove("../output/example_input_035.5000.nc")
except: pass

subprocess.run([pamtraBinary, "-f", ",".join(map(str,freqs)),
                "-d", descriptorFile,
                "-p", "../profile/example_input.lev",
                "-n", "pamtra_vs_pyPamtra.nml",
                "-o", "../output"], check=True)
ncData = netCDF4.Dataset("../output/example_input_035.5000.nc")
Ze = np.ma.masked_equal(ncData.variables["Ze"],-9999).ravel()

plt.figure(1)
plt.title("tb")
plt.plot(ncData.variables["tb"][:].ravel(),":",lw=4,label="Pamtra level")
plt.figure(2)
plt.title("Ze")
plt.plot(Ze.compressed(),":",lw=4,label="Pamtra level")

print(Ze.compressed())
ncData.close()

print("##########################")

subprocess.run([pamtraBinary, "-f", ",".join(map(str,freqs)),
                "-d", descriptorFile,
                "-p", "../profile/example_input.lay",
                "-n", "pamtra_vs_pyPamtra.nml",
                "-o", "../output"], check=True)
ncData = netCDF4.Dataset("../output/example_input_035.5000.nc")
Ze = np.ma.masked_equal(ncData.variables["Ze"],-9999).ravel()

plt.figure(1)
plt.title("tb")
plt.plot(ncData.variables["tb"][:].ravel(),":",lw=4,label="Pamtra layer")
plt.legend()

plt.figure(2)
plt.title("Ze")
plt.plot(Ze.compressed(),":",lw=4,label="Pamtra layer")
print(Ze.compressed())
plt.legend()
ncData.close()

#for nx in [0,1]:
  #for ny in [0,1]:
    #for nz in range(pam.p["max_nlyrs"]):
      #print nx, ny, nz, pam.p["nlyrs"][nx,ny], pam.p["hgt_lev"][nx,ny,nz]
plt.show()

