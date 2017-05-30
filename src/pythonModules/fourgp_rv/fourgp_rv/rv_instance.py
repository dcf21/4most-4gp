#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from multiprocessing import cpu_count
import logging
import itertools

import fourgp_speclib

logger = logging.getLogger(__name__)


class RvInstance(object):
    """
    A class which is adapted from Brani's code for forward-modelling the radial velocities of stars.
    """

    # Define the grid of template spectra supplied by Brani for use by RV code
    # [Axis variable, (min,max,step)]
    grid_axes = [["Teff", (4000, 8250, 250)],
                 ["Fe/H", (0.5, 3.0, 0.5)],
                 ["logg", (1.5, 5.5, 0.5)]
                 ]

    grid_axis_values = [np.arange(axis[1][0], axis[1][1], axis[1][2])
                        for axis in grid_axes]

    expected_number_spectra = np.prod([i.shape[0] for i in grid_axis_values])

    def __init__(self, spectrum_library, threads=None):
        """
        Instantiate the RV code, and read from disk the library of template spectra.
        
        :param spectrum_library:
            A SpectrumLibrary containing the template spectra we use for modelling.

        :type spectrum_library:
            SpectrumLibrary

        :param threads:
            The number of concurrent CPU threads to use.

        :type threads:
            int
        """

        assert isinstance(spectrum_library, fourgp_speclib.SpectrumLibrary), \
            "Argument to RvInstance should be a SpectrumLibrary."
        self._spectrum_library = spectrum_library

        # Work out how many CPUs we should allow the Cannon to use
        if threads is None:
            threads = cpu_count()
        self._threads = threads

        # Load template spectra
        spectrum_list = self._spectrum_library.search()

        # Check that we have the expected number of template spectra
        assert len(spectrum_list) == self.expected_number_spectra, \
            "Supplied SpectrumLibrary did not contain the expected number of spectra. " \
            "Was expecting {} spectra; got {} spectra.". \
            format(self.expected_number_spectra, len(spectrum_list))

        # Check that the template spectra are sampled on the right grid
        grid_axis_indices = [range(int((axis[1][1] - axis[1][0]) / axis[1][2])) for axis in self.grid_axes]
        grid_axis_index_combinations = itertools.product(*grid_axis_indices)

        # Check that each grid point in turn exists
        for axis_indices in grid_axis_index_combinations:
            search_criteria = {}
            for axis_counter, index in enumerate(axis_indices):
                metadata_key = self.grid_axes[axis_counter][0] + "_index"
                search_criteria[metadata_key] = index
            matches = self._spectrum_library.search(**search_criteria)
            assert len(matches) == 1, "Could not find spectrum matching {}".format(search_criteria)

    @classmethod
    def from_spectrum_library_sqlite(cls, library_path, *args, **kwargs):
        """
        Open a Spectrum Library on disk at supplied file path, an instantiate an RvInstance using the template
        spectra in it.

        :param library_path:
            The file path of the library of template spectra to use.

        :type library_path:
            str

        :return:
            RvInstance
        """
        library = fourgp_speclib.SpectrumLibrarySqlite(path=library_path)
        return cls(spectrum_library=library, *args, **kwargs)
