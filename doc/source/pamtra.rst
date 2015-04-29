..  _pamtra:

Pamtra
======

.. toctree::
   :maxdepth: 1

   pamtraDev    

Usage: pamtra [options]

Available options:

* -n|--namelist     namelist file (default None)
* -p|--profile      profile file  (default profile/standard.dat)
* -o|--output       profile directory  (default output)
* -d|--descriptor   descriptor file  (default descriptor_file.txt)
* -f|--freqs        comma seperated list of frequencies (no blanks) (default 89.0)
* -v|--verbose      integer specifying verbose level between -1 and 4 (default 0)
* -h|--help         print this help

Example: ./pamtra -v 1 -p rt_comp_single.dat -n run_params.nml -f 35,94,
See namelist file for further pamtra options