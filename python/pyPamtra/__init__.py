# -*- coding: utf-8 -*-

from __future__ import division

from .core import pyPamtra


# import subpackages
from . import meteoSI
from . import importer
from . import plot
from . import tools
from . import descriptorFile

# Version is pulled from git tag!!
from importlib.metadata import version, PackageNotFoundError
try:
    __version__ = ".".join(version("pyPamtra").split(".")[:3])
    __versionFull__ = version("pyPamtra")
except PackageNotFoundError:
    # package is not installed
    __version__ = "NotAvailable"
    __versionFull__ = "NotAvailable"
    pass