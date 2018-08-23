#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

import logging
from numpy import RankWarning
from warnings import simplefilter

from .convolve import SpectrumConvolver
from .interpolate import SpectrumInterpolator
from .gaussian_noise import GaussianNoise
from .resample import SpectrumResampler
from .redden import SpectrumReddener
from .snr_conversion import SNRConverter, SNRValue
from .spectrum_properties import SpectrumProperties

__version__ = "release-2019-09-01-01"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # TODO: Remove this when stable.

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    "%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(handler)

simplefilter("ignore", RankWarning)
simplefilter("ignore", RuntimeWarning)


