#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for all SQL implementations of spectrum libraries.
"""

from os import path as os_path
import uuid
import unittest
import numpy as np
import fourgp_speclib


class TestSpectrumLibrarySQL(object):
    """
    This class is a mixin which adds lots of standard tests to any SQL-based SpectrumLibrary unittest class.
    """

    def test_refresh(self):
        """
        Check that we can refresh database connection.
        """
        self._lib.refresh_database()

    def test_spectrum_retrieval(self):
        """
        Check that we can store a single spectra into the SpectrumLibrary and retrieve it again.
        """

        # Create a random spectrum to insert into the spectrum library
        size = 50
        raster = np.arange(size)
        values = np.random.random(size)
        value_errors = np.random.random(size)
        input_spectrum = fourgp_speclib.Spectrum(wavelengths=raster,
                                                 values=values,
                                                 value_errors=value_errors,
                                                 metadata={"origin": "unit-test"})

        # Insert it into the spectrum library
        self._lib.insert(input_spectrum, "dummy_filename")

        # Load it back as a SpectrumArray
        my_spectrum_array = self._lib.open(filenames=["dummy_filename"])

        # Pick spectrum out of SpectrumArray
        my_spectrum = my_spectrum_array.extract_item(0)

        # Check that we got back the same spectrum we put in
        self.assertEqual(my_spectrum, input_spectrum)
