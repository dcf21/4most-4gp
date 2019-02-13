#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
from numpy import RankWarning
from warnings import simplefilter

__version__ = "20190301.1"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # TODO: Remove this when stable.

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(handler)

simplefilter("ignore", RankWarning)
simplefilter("ignore", RuntimeWarning)
