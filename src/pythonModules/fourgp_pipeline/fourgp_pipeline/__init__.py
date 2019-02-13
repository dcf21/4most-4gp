#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This python module defines the Pipeline implemented by 4GP.
"""

import logging
from numpy import RankWarning
from warnings import simplefilter

from .pipeline_manager import PipelineManager
from .pipeline import Pipeline
from .spectrum_analysis import SpectrumAnalysis

__version__ = "20190301.1"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # TODO: Remove this when stable.

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)-8s] %(message)s"))
logger.addHandler(handler)

simplefilter("ignore", RankWarning)
simplefilter("ignore", RuntimeWarning)
