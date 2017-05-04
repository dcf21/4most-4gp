#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This python module defines classes representing arrays of spectra, and libraries to keep them in.
"""

__version__ = "0.1.0"

import logging
from numpy import RankWarning
from warnings import simplefilter

from . import spectrum_library_sqlite, spectrum, spectrum_polynomial

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # TODO: Remove this when stable.

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    "%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(handler)

simplefilter("ignore", RankWarning)
simplefilter("ignore", RuntimeWarning)


