import os
import shlex

import sys
from glob import glob

import numpy
import setuptools
from Cython.Build import cythonize
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension
from numpy.distutils.misc_util import Configuration

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


# try:
#     os.environ['LDFLAGS'] += " -shared"
# except KeyError:
#     os.environ['LDFLAGS'] = " -shared"

kw = {}
if sys.platform == 'darwin':
    kw['extra_link_args'] = ['-undefined dynamic_lookup', '-bundle']
    os.environ["CC"] = "/usr/local/bin/gcc-11"

library_dirs = ['~/.local/lib/', '/usr/local/lib/', '/opt/local/lib/']

def configuration(parent_package='', top_path=None):

    config = Configuration(
        'pyPamtra',
        parent_package,
        top_path,
        version='1.1',
        author="Pamtra Team",
        author_email="meteo-pamtra@uni-koeln.de",
        description="atmospheric microwace passive"
                    " and active instrument simulator",
        license="GPL v3",
        python_requires='>=3.5',
        url='https://github.com/igmk/pamtra',
        long_description=read('readme.md'),
        classifiers=[
            "Development Status :: 2 - Beta",
            "License :: OSI Approved :: GPL v3 License",
            "Operating System :: OS Independent",
            "Programming Language :: Fortran",
            "Programming Language :: Python",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: Atmospheric Science",
        ]
    )
    return config


fortPath = 'src'
pyPamtra = Extension(
    'pyPamtra.pyPamtraLib',
    sources=[
        '%s/pyPamtraLib.pyf' %fortPath,
        "%s/kinds.f90"%fortPath,
        '%s/pyPamtraLib.f90' %fortPath,
        "%s/nan.f90"%fortPath,
        "%s/vars_index.f90"%fortPath,
        "%s/report_module.f90"%fortPath,
        "%s/rt_utilities.f90"%fortPath,
        "%s/settings.f90"%fortPath,
        "%s/constants.f90"%fortPath,
        "%s/mod_fastem4_coef.f90"%fortPath,
        "%s/gasabs_module.f90"%fortPath,
        "%s/conversions.f90"%fortPath,
        "%s/descriptor_file.f90"%fortPath,
        "%s/vars_atmosphere.f90"%fortPath,
        "%s/vars_rt.f90"%fortPath,
        "%s/vars_hydroFullSpec.f90"%fortPath,
        "%s/mod_io_strings.f90"%fortPath,
        "%s/getopt.f90"%fortPath,
        "%s/parse_options.f90"%fortPath,
        "%s/radmat.f90"%fortPath,
        "%s/convolution.f90"%fortPath,
        "%s/get_gasabs.f90"%fortPath,
        "%s/vars_output.f90"%fortPath,
        "%s/azimuth_emissivity_module.f90"%fortPath,
        "%s/hyperbolic_step.f90"%fortPath,
        "%s/slope_variance.f90"%fortPath,
        "%s/reflection_correction_module.f90"%fortPath,
        "%s/large_scale_correction_module.f90"%fortPath,
        "%s/small_scale_correction_module.f90"%fortPath,
        "%s/foam_utility_module.f90"%fortPath,
        "%s/liu.f90"%fortPath,
        "%s/fresnel.f90"%fortPath,
        "%s/fastemx.f90"%fortPath,
        "%s/tessem2.f90"%fortPath,
        "%s/ocean_sfc_optics.f90"%fortPath,
        "%s/mod_mwatlas_nt_bin.f90"%fortPath,
        "%s/telsem2.f90"%fortPath,
        "%s/land_sfc_optics.f90"%fortPath,
        "%s/sfc_optics.f90"%fortPath,
        "%s/sfc_matrices.f90"%fortPath,
        "%s/run_rt.f90"%fortPath,
        "%s/scat_utilities.f90"%fortPath,
        "%s/mpm93.f90"%fortPath,
        "%s/eps_water.f90"%fortPath,
        "%s/mie_scat_utilities.f90"%fortPath,
        "%s/mie_spheres.f90"%fortPath,
        "%s/scatdb.c"%fortPath,
        "%s/dda_db_liu.f90"%fortPath,
        "%s/dda_db_hong.f90"%fortPath,
        "%s/hongdb.f90"%fortPath,
        "%s/dia2vel.f90"%fortPath,
        "%s/rescale_spectra.f90"%fortPath,
        "%s/radar_moments.f90"%fortPath,
        "%s/radar_spectrum.f90"%fortPath,
        "%s/radar_spectral_broadening.f90"%fortPath,
        "%s/radar_simulator.f90"%fortPath,
        "%s/rosen98_gasabs.f90"%fortPath,
        "%s/surface.f90"%fortPath,
        "%s/eps_ice.f90"%fortPath,
        "%s/eps_mix.f90"%fortPath,
        "%s/equcom.f90"%fortPath,
        "%s/land_emis_ssmi.f90"%fortPath,
        "%s/equare.f90"%fortPath,
        "%s/ref_water.f90"%fortPath,
        "%s/ref_ice.f90"%fortPath,
        "%s/e_sat_gg_water.f90"%fortPath,
        "%s/interpolation.f90"%fortPath,
        "%s/collect_output.f90"%fortPath,
        "%s/save_active.f90"%fortPath,
        "%s/random.f90"%fortPath,
        "%s/rt4.f90"%fortPath,
        "%s/radtran4.f90"%fortPath,
        "%s/radintg4.f90"%fortPath,
        "%s/radscat4.f90"%fortPath,
        "%s/drop_size_dist.f90"%fortPath,
        "%s/make_dist.f90"%fortPath,
        "%s/make_dist_param.f90"%fortPath,
        "%s/make_mass_size.f90"%fortPath,
        "%s/calc_moment.f90"%fortPath,
        "%s/make_soft_spheroid.f90"%fortPath,
        "%s/check_print.f90"%fortPath,
        "%s/tmatrix.f90"%fortPath,
        "%s/rayleigh_gans.f90"%fortPath,
        "%s/scatProperties.f90"%fortPath,
        "%s/hydrometeor_extinction.f90"%fortPath,
        "%s/scatcnv.f90"%fortPath,
        "%s/tmatrix_lpq.f"%fortPath,
        "%s/get_scat_mat.f90"%fortPath,
        "%s/refractive_index.f90"%fortPath,
        "%s/dsort.f90"%fortPath,
        "%s/rho_air.f90"%fortPath,
        "%s/viscosity_air.f90"%fortPath,
        "%s/versionNumber.auto.f90"%fortPath,
        "%s/smooth_savitzky_golay.f90"%fortPath,
        "%s/radar_hildebrand_sekhon.f90"%fortPath,
        "%s/tmatrix_amplq.lp.f"%fortPath,
        "%s/deallocate_everything.f90"%fortPath,
        "%s/zlib_stuff.f90"%fortPath, 
    ],
    library_dirs=library_dirs,
    libraries=['fftw3', 'lapack', "z", "blas"],
    extra_compile_args = [ "-fPIC", "-cpp", "-c"],
    **kw)


usStandardPath = "tools/py_usStandard"
usStandard = Extension(
    'usStandard.usStandardAtmosphere',
    sources=[
        '%s/src/usStandardAtmosphere.pyf' %usStandardPath,
        "%s/src/usStandard.f90"%usStandardPath,
    ],
    library_dirs=library_dirs,
    extra_compile_args = [ "-fPIC", "-cpp", "-c"],
    **kw)

# cMie = Extension(
#     name = "pamtra2.libs.singleScattering.cMie",
#     sources = ["%s/Mie/cython/cMie.pyx" % singleScattering_path,
#              "%s/Mie/src/cMie.c" % singleScattering_path],
#     include_dirs = [numpy.get_include()],
#     extra_compile_args = ["-O3", "-ffast-math",
#                           "-Wall", "-lm", "-fPIC", "-std=c99"],
#     language='c'
# )


if __name__ == "__main__":

    setup(
        configuration=configuration,
        packages=[
        'pyPamtra',
        'pyPamtra.fortranNamelist',
        'usStandard'
        ],
        package_dir={
            'pyPamtra': 'python/pyPamtra' ,
            'usStandard': '%s/usStandard'%usStandardPath ,
        },
        package_data={
            # Don't aks me why, but I need an extra random character in front
            # of the *dat files...
            # 'pamtra2.libs.refractiveIndex': ['*.dat'],
            # 'pamtra2.libs.singleScattering': ['*.dat'],
        },
        platforms=['any'],
        install_requires=["pandas", "numpy", "scipy", "matplotlib", "netcdf4",
                  'cython'],
        build_requires=['numpy', 'cython'],
        setup_requires=["pytest-runner"],
        tests_require=["pytest"],
        ext_modules=cythonize(
            [pyPamtra,usStandard]),

    )
