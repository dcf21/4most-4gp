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

    def test_search_illegal_metadata(self):
        """
        Check that we can search for spectra on a simple metadata constraint.
        """

        # Insert ten random spectra into SpectrumLibrary
        size = 50
        input_spectrum = fourgp_speclib.Spectrum(wavelengths=np.arange(size),
                                                 values=np.random.random(size),
                                                 value_errors=np.random.random(size),
                                                 metadata={"origin": "unit-test"})
        self._lib.insert(input_spectrum, "dummy_filename")

        # Search on an item of metadata which doesn't exist
        with self.assertRaises(AssertionError):
            my_spectra = self._lib.search(x_value=23)

    def test_search_1d_numerical_range(self):
        """
        Check that we can search for spectra on a simple metadata numerical range constraint.
        """

        # Insert ten random spectra into SpectrumLibrary
        size = 50
        x_values = range(10)
        for x in x_values:
            input_spectrum = fourgp_speclib.Spectrum(wavelengths=np.arange(size),
                                                     values=np.random.random(size),
                                                     value_errors=np.random.random(size),
                                                     metadata={"origin": "unit-test",
                                                               "x_value": x})
            self._lib.insert(input_spectrum, "x_{}".format(x))

        # Search for spectra with x in a defined range
        x_range = [4.5, 8.5]
        filenames_expected = ["x_{}".format(x) for x in x_values if (x > x_range[0] and x < x_range[1])]
        my_spectra = self._lib.search(x_value=x_range)
        filenames_got = [str(item["filename"]) for item in my_spectra]
        filenames_got.sort()

        # Check that we got back the same spectrum we put in
        self.assertEqual(filenames_expected, filenames_got)

    def test_search_1d_numerical_value(self):
        """
        Check that we can search for spectra on a simple metadata numerical point-value constraint.
        """

        # Insert ten random spectra into SpectrumLibrary
        size = 50
        x_values = range(10)
        for x in x_values:
            input_spectrum = fourgp_speclib.Spectrum(wavelengths=np.arange(size),
                                                     values=np.random.random(size),
                                                     value_errors=np.random.random(size),
                                                     metadata={"origin": "unit-test",
                                                               "x_value": x})
            self._lib.insert(input_spectrum, "x_{}".format(x))

        # Search for spectra with matching x_value
        my_spectra = self._lib.search(x_value=5)
        filenames_got = [str(item["filename"]) for item in my_spectra]

        # Check that we got back the same spectrum we put in
        self.assertEqual(filenames_got, ["x_5"])

    def test_search_1d_string_range(self):
        """
        Check that we can search for spectra on a simple metadata string range constraint.
        """

        # Insert random spectra into SpectrumLibrary
        alphabet = "abcdefghijklmnopqrstuvwxyz"
        size = 50
        x_values = range(12)
        for x in x_values:
            input_spectrum = fourgp_speclib.Spectrum(wavelengths=np.arange(size),
                                                     values=np.random.random(size),
                                                     value_errors=np.random.random(size),
                                                     metadata={"origin": "unit-test",
                                                               "x_value": alphabet[x:x + 3]})
            self._lib.insert(input_spectrum, "x_{}".format(x))

        # Search for spectra with x in a defined range
        my_spectra = self._lib.search(x_value=["dxx", "h"])
        x_values_expected = ["efg", "fgh", "ghi"]
        filenames_got = [str(item["filename"]) for item in my_spectra]
        x_values_got = [str(i["x_value"]) for i in self._lib.get_metadata(filenames=filenames_got)]
        x_values_got.sort()

        # Check that we got back the same spectrum we put in
        self.assertEqual(x_values_expected, x_values_got)

    def test_search_1d_string_value(self):
        """
        Check that we can search for spectra on a simple metadata string point-value constraint.
        """

        # Insert random spectra into SpectrumLibrary
        alphabet = "abcdefghijklmnopqrstuvwxyz"
        size = 50
        x_values = range(10)
        for x in x_values:
            input_spectrum = fourgp_speclib.Spectrum(wavelengths=np.arange(size),
                                                     values=np.random.random(size),
                                                     value_errors=np.random.random(size),
                                                     metadata={"origin": "unit-test",
                                                               "x_value": alphabet[x:x + 3]})
            self._lib.insert(input_spectrum, "x_{}".format(x))

        # Search for spectra with matching x_value
        my_spectra = self._lib.search(x_value="def")
        filenames_got = [str(item["filename"]) for item in my_spectra]
        x_values_got = [str(i["x_value"]) for i in self._lib.get_metadata(filenames=filenames_got)]
        x_values_got.sort()

        # Check that we got back the same spectrum we put in
        self.assertEqual(x_values_got, ["def"])
