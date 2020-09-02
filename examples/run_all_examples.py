from __future__ import print_function

import os
import glob
import numpy as np

for fname in np.sort(glob.glob("*.py")):
  if fname == "run_all_examples.py" : continue
  print("#"*60)
  print(fname)
  print("#"*60)
  exec(compile(open(fname).read(), fname, 'exec'))