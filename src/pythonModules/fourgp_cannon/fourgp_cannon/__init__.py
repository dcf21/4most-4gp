#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
from numpy import RankWarning
from warnings import simplefilter

# from .cannon_instance import \
#    CannonInstance, \
#    CannonInstanceWithContinuumNormalisation, \
#    CannonInstanceWithRunningMeanNormalisation

from .cannon_instance_release_2018_01_09_1 import \
    CannonInstance_2018_01_09, \
    CannonInstanceWithRunningMeanNormalisation_2018_01_09, \
    CannonInstanceWithContinuumNormalisation_2018_01_09

CannonInstance = CannonInstanceWithContinuumNormalisation = CannonInstanceWithRunningMeanNormalisation = None

__version__ = "release-2019-09-01-01"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # TODO: Remove this when stable.

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(handler)

simplefilter("ignore", RankWarning)
simplefilter("ignore", RuntimeWarning)
