#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

import logging
from numpy import RankWarning
from warnings import simplefilter

from .convolve import SpectrumConvolver
from .interpolate import SpectrumInterpolator
from .resample import SpectrumResampler
from .gaussian_noise import GaussianNoise
from .redden import SpectrumReddener

__version__ = "0.1.0"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # TODO: Remove this when stable.

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    "%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(handler)

simplefilter("ignore", RankWarning)
simplefilter("ignore", RuntimeWarning)


