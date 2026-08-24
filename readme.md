# PAMTRA package - Passive and Active Microwave TRANsfer 

[![Documentation Status](https://readthedocs.org/projects/pamtra/badge/?version=latest)](https://pamtra.readthedocs.io/en/latest/?badge=latest)
[![CI](https://github.com/igmk/pamtra/actions/workflows/ci.yml/badge.svg)](https://github.com/igmk/pamtra/actions/workflows/ci.yml)
[![Coverage Status](https://coveralls.io/repos/github/igmk/pamtra/badge.svg?branch=master)](https://coveralls.io/github/igmk/pamtra?branch=master)
[![PyPI Downloads](https://static.pepy.tech/badge/pamtra)](https://pepy.tech/projects/pamtra)


Python/Fortran 90 package to solve the passive and active microwave radiative transfer in a plan parallel horizontally homogeneous atmosphere with hydrometeors

## Manual and Installation

Recommended, for Linux (x86_64) and macOS (Apple Silicon / arm64):

```
pip install pamtra
```

(If you previously used the old `make pyinstall` workflow, see the
installation docs below for a warning about removing its leftover
`~/lib/python/pyPamtra` by hand -- it can silently shadow this.)

Not available as a wheel for Windows or Intel macOS -- see
https://pamtra.readthedocs.io/en/latest/installation.html for those and
other install options (conda-forge/pixi, building from source, HPC
clusters).

## Mailing list

If you want to join the mailing list you can find it here https://lists.uni-koeln.de/mailman/listinfo/meteo-pamtra. There you can ask any question related to the usage of PAMTRA.

## Reference

If you use PAMTRA, please cite our paper:

Mech, M., Maahn, M., Kneifel, S., Ori, D., Orlandi, E., Kollias, P., Schemann, V., and Crewell, S.: PAMTRA 1.0: the Passive and Active Microwave radiative TRAnsfer tool for simulating radiometer and radar measurements of the cloudy atmosphere, Geosci. Model Dev., 13, 4229–4251, https://doi.org/10.5194/gmd-13-4229-2020, 2020.

