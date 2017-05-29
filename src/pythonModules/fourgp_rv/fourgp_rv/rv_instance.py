#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from multiprocessing import cpu_count
import logging
from astropy.table import Table

import fourgp_speclib

logger = logging.getLogger(__name__)


class RvInstance(object):
    """
    A class which is adapted from Brani's code for forward-modelling the radial velocities of stars.
    """

    def __init__(self, spectrum_library, threads=None):
        """
        Instantiate the RV code, and read from disk the library of template spectra.
        
        :param spectrum_library:
            A SpectrumLibrary containing the template spectra we use for modelling.
        """

        assert isinstance(training_set, fourgp_speclib.SpectrumLibrary), \
            "Argument to RvInstance should be a SpectrumLibrary."
        self._spectrum_library = spectrum_library

        # Work out how many CPUs we should allow the Cannon to use
        if threads is None:
            threads = cpu_count()

