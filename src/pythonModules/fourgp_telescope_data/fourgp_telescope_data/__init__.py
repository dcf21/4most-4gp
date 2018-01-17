#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
This python package is a small collection of data about 4MOST -- e.g. its wavelength coverage and spectral resolution.
"""

import logging
from numpy import RankWarning
from warnings import simplefilter

from .fourmost import FourMost

__version__ = "0.1.0"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # TODO: Remove this when stable.

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(handler)

simplefilter("ignore", RankWarning)
simplefilter("ignore", RuntimeWarning)


