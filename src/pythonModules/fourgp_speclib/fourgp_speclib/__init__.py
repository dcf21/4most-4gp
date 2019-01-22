#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This python module defines classes representing arrays of spectra, and libraries to keep them in.
"""

import logging
from numpy import RankWarning
from warnings import simplefilter

from .spectrum_library_sqlite import SpectrumLibrarySqlite
from .spectrum_library import SpectrumLibrary
from .spectrum_array import SpectrumArray
from .spectrum import Spectrum, spectrum_splice
from .spectrum_smooth import SpectrumSmoothFactory, SpectrumSmooth, SpectrumPolynomial

# Allow MySQL binding to silently fail if system doesn't have MySQLdb package installed
try:
    from .spectrum_library_mysql import SpectrumLibraryMySql
except ImportError:
    pass

__version__ = "20190201.1"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # TODO: Remove this when stable.

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(handler)

simplefilter("ignore", RankWarning)
simplefilter("ignore", RuntimeWarning)
