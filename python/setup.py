import os
import sys

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

print('#'*75)
print('WARNING '*10)
print('This setup file is currently only for readthedocs! '
  'It does not compile the Fortran module')
print('WARNING '*10)
print('#'*75)


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def configuration(parent_package='', top_path=None):

    config = Configuration(
        'Pamtra', parent_package, top_path,
        version='0.1',
        author="IGMK",
        license="GPL v3",
        python_requires='2.7',
        # url='https://github.com/maahn/pyPamtraRadarSimulator',
        # download_url='https://github.com/maahn/pyPamtraRadarSimulator/releases/download/0.1/pyPamtraRadarSimulator-0.1.zip',
        long_description=read('../readme.md'),
        classifiers=[
           "Development Status :: 3 - Alpha",
           "License :: OSI Approved :: MIT License",
           "Operating System :: OS Independent",
           "Programming Language :: Fortran",
           "Programming Language :: Python",
           "Intended Audience :: Science/Research",
           "Topic :: Scientific/Engineering :: Atmospheric Science",
        ]
        )

    kw = {}
  #   if sys.platform == 'darwin':
  #       kw['extra_link_args'] = ['-undefined dynamic_lookup', '-bundle']
  #   config.add_extension('pyPamtraLib',
  #                        sources=[
  # 'src/kinds.f90',
  # 'src/nan.f90',
  # 'src/vars_index.f90',
  # 'src/report_module.f90',
  # 'src/rt_utilities.f90',
  # 'src/settings.f90',
  # 'src/constants.f90',
  # 'src/zlib_stuff.f90',
  # 'src/mod_fastem4_coef.f90',
  # 'src/gasabs_module.f90',
  # 'src/conversions.f90',
  # 'src/descriptor_file.f90',
  # 'src/vars_atmosphere.f90',
  # 'src/vars_rt.f90',
  # 'src/vars_hydroFullSpec.f90',
  # 'src/mod_io_strings.f90',
  # 'src/getopt.f90',
  # 'src/parse_options.f90',
  # 'src/radmat.f90',
  # 'src/convolution.f90',
  # 'src/get_gasabs.f90',
  # 'src/vars_output.f90',
  # 'src/azimuth_emissivity_module.f90',
  # 'src/hyperbolic_step.f90',
  # 'src/slope_variance.f90',
  # 'src/reflection_correction_module.f90',
  # 'src/large_scale_correction_module.f90',
  # 'src/small_scale_correction_module.f90',
  # 'src/foam_utility_module.f90',
  # 'src/liu.f90',
  # 'src/fresnel.f90',
  # 'src/fastemx.f90',
  # 'src/tessem2.f90',
  # 'src/ocean_sfc_optics.f90',
  # 'src/mod_mwatlas_nt_bin.f90',
  # 'src/telsem2.f90',
  # 'src/land_sfc_optics.f90',
  # 'src/sfc_optics.f90',
  # 'src/sfc_matrices.f90',
  # 'src/run_rt.f90',
  # 'src/scat_utilities.f90',
  # 'src/mpm93.f90',
  # 'src/eps_water.f90',
  # 'src/mie_scat_utilities.f90',
  # 'src/mie_spheres.f90',
  # 'src/scatdb.f90',
  # 'src/dda_db_liu.f90',
  # 'src/dda_db_hong.f90',
  # 'src/hongdb.f90',
  # 'src/dia2vel.f90',
  # 'src/rescale_spectra.f90',
  # 'src/radar_moments.f90',
  # 'src/radar_spectrum.f90',
  # 'src/radar_spectral_broadening.f90',
  # 'src/radar_simulator.f90',
  # 'src/rosen98_gasabs.f90',
  # 'src/surface.f90',
  # 'src/eps_ice.f90',
  # 'src/eps_mix.f90',
  # 'src/equcom.f90',
  # 'src/land_emis_ssmi.f90',
  # 'src/equare.f90',
  # 'src/ref_water.f90',
  # 'src/ref_ice.f90',
  # 'src/e_sat_gg_water.f90',
  # 'src/interpolation.f90',
  # 'src/collect_output.f90',
  # 'src/save_active.f90',
  # 'src/random.f90',
  # 'src/rt4.f90',
  # 'src/radtran4.f90',
  # 'src/radintg4.f90',
  # 'src/radscat4.f90',
  # 'src/drop_size_dist.f90',
  # 'src/make_dist.f90',
  # 'src/make_dist_param.f90',
  # 'src/make_mass_size.f90',
  # 'src/calc_moment.f90',
  # 'src/make_soft_spheroid.f90',
  # 'src/check_print.f90',
  # 'src/tmatrix.f90',
  # 'src/rayleigh_gans.f90',
  # 'src/scatProperties.f90',
  # 'src/hydrometeor_extinction.f90',
  # 'src/scatcnv.f90',
  # 'src/tmatrix_lpq.f90',
  # 'src/get_scat_mat.f90',
  # 'src/refractive_index.f90',
  # 'src/dsort.f90',
  # 'src/rho_air.f90',
  # 'src/viscosity_air.f90',
  # 'src/versionNumber.auto.f90',
  # 'src/smooth_savitzky_golay.f90',
  # 'src/radar_hildebrand_sekhon.f90',
  # 'src/tmatrix_amplq.lp.f90',
  # 'src/deallocate_everything.f90',
  #                        ],
  #                        library_dirs=['~/.local/lib/', '/usr/local/lib/'],
  #                        libraries=['fftw3', 'lapack'],
  #                        **kw)

    return config


if __name__ == "__main__":

    setup(configuration=configuration,
          package_dir = {'pyPamtra': 'pyPamtra'},
          packages=['pyPamtra',
                    # 'pyPamtra.pyPamtraLib'
                    ],
          platforms=['any'],
          requires=['numpy', 'scipy' ,'matplotlib', 'netCDF4'])
