#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
from numpy import RankWarning
from warnings import simplefilter

from .rv_random import random_radial_velocity
from .brani_code import RvInstanceBrani
from .cross_correlation import RvInstanceCrossCorrelation

__version__ = "20190201.1"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # TODO: Remove this when stable.

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(handler)

simplefilter("ignore", RankWarning)
simplefilter("ignore", RuntimeWarning)
