#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "0.1.0"

import logging
from numpy import RankWarning
from warnings import simplefilter

from .rv_instance import RvInstance

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # TODO: Remove this when stable.

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    "%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(handler)

simplefilter("ignore", RankWarning)
simplefilter("ignore", RuntimeWarning)


