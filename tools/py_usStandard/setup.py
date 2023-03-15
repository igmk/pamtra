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
        'usStandard',
        parent_package,
        top_path,
        version='2.0',
        author="Pamtra Team",
        author_email="meteo-pamtra@uni-koeln.de",
        description="US standard Atmosphere",
        license="GPL v3",
        python_requires='>=3.5',
        url='https://github.com/maahn/py_usStandard',
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

usStandard = Extension(
    'usStandard.usStandardAtmosphere',
    sources=[
        '%s/src/usStandardAtmosphere.pyf' %rootPath,
        "%s/src/usStandard.f90"%rootPath,
    ],
    library_dirs=library_dirs,
    extra_compile_args = [ "-fPIC", "-cpp", "-c"],
    **kw)

rootPath = './'
if __name__ == "__main__":

    setup(
        configuration=configuration,
        packages=[
        'usStandard',
                  ],
        package_dir={
            'usStandard': '%s/usStandard'%rootPath ,
        },
        package_data={
        },
        platforms=['any'],
        install_requires=["numpy", 'cython'],
        build_requires=['numpy', 'cython'],
        ext_modules=cythonize(
            [usStandard]),

    )
