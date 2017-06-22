#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This python module defines a wrapper for synthesising spectra using TurboSpectrum.
"""

import logging
from numpy import RankWarning
from warnings import simplefilter

from .turbospectrum import TurboSpectrum
from .solar_abundances import solar_abundances

__version__ = "0.1.0"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # TODO: Remove this when stable.

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    "%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(handler)

simplefilter("ignore", RankWarning)
simplefilter("ignore", RuntimeWarning)


