#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
Unit tests for the Spectrum class
"""

import os
from os import path as os_path
import uuid
import unittest
import numpy as np
import fourgp_speclib


class TestSpectrum(unittest.TestCase):
    def setUp(self):
        """
        Create a Spectrum object.
        """

        self._size = 50
        self._raster = np.arange(self._size)
        self._values = np.arange(100, self._size + 100)
        self._value_errors = np.random.random(self._size)
        self._spectrum = fourgp_speclib.Spectrum(wavelengths=self._raster,
                                                 values=self._values,
                                                 value_errors=self._value_errors,
                                                 metadata={"origin": "unit-test"})

    def test_metadata(self):
        self.assertEqual(self._spectrum.metadata["origin"], "unit-test")

    def test_length(self):
        self.assertEqual(len(self._spectrum), self._size)

    def test_data_sizes_must_match_1(self):
        with self.assertRaises(AssertionError):
            fourgp_speclib.Spectrum(wavelengths=np.arange(self._size + 1),
                                    values=self._values,
                                    value_errors=self._value_errors)

    def test_data_sizes_must_match_2(self):
        with self.assertRaises(AssertionError):
            fourgp_speclib.Spectrum(wavelengths=self._raster,
                                    values=np.arange(self._size + 1),
                                    value_errors=self._value_errors)

    def test_data_sizes_must_match_3(self):
        with self.assertRaises(AssertionError):
            fourgp_speclib.Spectrum(wavelengths=self._raster,
                                    values=self._values,
                                    value_errors=np.arange(self._size + 1))

    def test_spectrum_retrieval_binary(self):
        unique_filename = str(uuid.uuid4())
        unique_path = os_path.join("/tmp", unique_filename+".npy")
        self._spectrum.to_file(unique_path, binary=True)
        new_spectrum = self._spectrum.from_file(unique_path, binary=True)
        os.unlink(unique_path)

        # Check that we got back the same spectrum we put in
        self.assertEqual(self._spectrum, new_spectrum)

    def test_spectrum_retrieval_text(self):
        unique_filename = str(uuid.uuid4())
        unique_path = os_path.join("/tmp", unique_filename+".txt")
        self._spectrum.to_file(unique_path, binary=False)
        new_spectrum = self._spectrum.from_file(unique_path, binary=False)
        os.unlink(unique_path)

        # Check that we got back the same spectrum we put in
        self.assertEqual(self._spectrum, new_spectrum)

    def test_spectrum_retrieval_unspecified_format(self):
        unique_filename = str(uuid.uuid4())
        unique_path = os_path.join("/tmp", unique_filename+".npy")
        self._spectrum.to_file(unique_path)
        new_spectrum = self._spectrum.from_file(unique_path)
        os.unlink(unique_path)

        # Check that we got back the same spectrum we put in
        self.assertEqual(self._spectrum, new_spectrum)

    def test_addition_multiplication(self):
        """
        Try adding spectra together repeatedly using the __sum__ and __isum__ methods. Check that this is the same
        as multiplying the spectrum by a fixed integer.
        """

        # Create an empty numpy array to insert raster of multipliers into
        multiplier = np.empty(self._size)
        failures = 0
        for i in range(2, 5):
            sum_1 = fourgp_speclib.Spectrum(wavelengths=self._raster,
                                            values=np.zeros(self._size),
                                            value_errors=self._value_errors)
            sum_2 = fourgp_speclib.Spectrum(wavelengths=self._raster,
                                            values=np.zeros(self._size),
                                            value_errors=self._value_errors)

            # Test __add__ method
            for j in range(i):
                sum_1 = sum_1 + self._spectrum

            # Test __iadd__ method
            for j in range(i):
                sum_2 += self._spectrum

            # Test __mul__ method
            multiplier.fill(i)
            b = fourgp_speclib.Spectrum(wavelengths=self._raster,
                                        values=multiplier,
                                        value_errors=self._value_errors)
            sum_3 = b * self._spectrum

            # Test __imul__ method
            b *= self._spectrum

            # Check that all three calculations reached the same result
            if sum_1 != sum_2:
                failures += 1
            if sum_1 != sum_3:
                failures += 1
            if sum_1 != b:
                failures += 1

            del sum_1, sum_2, sum_3, b

        # Check that none of the calculations failed
        self.assertEqual(failures, 0)

    def test_subtraction_division(self):
        """
        Try subtracting a spectrum from zero N times. Then divide by minus N times, and ensure we get back to where we
        started.
        """

        # Create an empty numpy array to insert raster of multipliers into
        multiplier = np.empty(self._size)
        failures = 0
        for i in range(2, 5):
            sum_1 = fourgp_speclib.Spectrum(wavelengths=self._raster,
                                            values=np.zeros(self._size),
                                            value_errors=self._value_errors)
            sum_2 = fourgp_speclib.Spectrum(wavelengths=self._raster,
                                            values=np.zeros(self._size),
                                            value_errors=self._value_errors)

            # Test __sub__ method
            for j in range(i):
                sum_1 = sum_1 - self._spectrum

            # Test __isub__ method
            for j in range(i):
                sum_2 -= self._spectrum

            # Test __div__ method
            multiplier.fill(-i)
            b = fourgp_speclib.Spectrum(wavelengths=self._raster,
                                        values=multiplier,
                                        value_errors=self._value_errors)
            sum_3 = sum_1 / b
            sum_4 = sum_2 / b

            # Test __idiv__ method
            sum_1 /= b
            sum_2 /= b

            # Check that all three calculations reached the same result
            for item in [sum_1, sum_2, sum_3, sum_4]:
                if item != self._spectrum:
                    failures += 1

            del sum_1, sum_2, sum_3, sum_4, b

        # Check that none of the calculations failed
        self.assertEqual(failures, 0)

    def tearDown(self):
        """
        Tear down Spectrum object.
        """
        del self._spectrum


# Run tests if we are run from command line
if __name__ == '__main__':
    unittest.main()
