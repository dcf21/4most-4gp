#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

import logging
from numpy import RankWarning
from warnings import simplefilter

from .rv_instance import RvInstance, random_radial_velocity

__version__ = "release-2019-09-01-01"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # TODO: Remove this when stable.

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(handler)

simplefilter("ignore", RankWarning)
simplefilter("ignore", RuntimeWarning)


