import os
import glob


for fname in np.sort(glob.glob("*.py")):
  if fname == "run_all_examples.py" : continue
  print "#"*60
  print fname
  print "#"*60
  execfile(fname)