#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
Unit tests for the SpectrumPolynomial class
"""

import unittest
import numpy as np
import fourgp_speclib


class TestSpectrumPolynomial(unittest.TestCase):
    def setUp(self):
        """
        Create test Spectrum objects.
        """

        self._size = 50
        self._observed_raster = np.arange(self._size, dtype=np.float64)
        self._observed_values = np.arange(self._size, dtype=np.float64)
        self._observed_value_errors = np.ones(self._size, dtype=np.float64)
        self._observed = fourgp_speclib.Spectrum(wavelengths=self._observed_raster,
                                                 values=self._observed_values,
                                                 value_errors=self._observed_value_errors,
                                                 metadata={"origin": "unit-test"})

        self._absorption_values = np.ones(self._size, dtype=np.float64)
        self._absorption_value_errors = np.ones(self._size, dtype=np.float64)
        self._absorption = fourgp_speclib.Spectrum(wavelengths=self._observed_raster,
                                                   values=self._absorption_values,
                                                   value_errors=self._absorption_value_errors,
                                                   metadata={"origin": "unit-test"})

        self._polynomial = fourgp_speclib.SpectrumSmoothFactory(function_family=fourgp_speclib.SpectrumPolynomial,
                                                                wavelengths=self._observed_raster,
                                                                terms=3)

    def test_data_sizes_must_match_1(self):
        with self.assertRaises(AssertionError):
            raster = np.arange(self._size + 10, dtype=np.float64)
            absorption = fourgp_speclib.Spectrum(wavelengths=raster,
                                                 values=raster,
                                                 value_errors=raster)
            self._polynomial.fit_to_continuum_via_template(other=self._observed, template=absorption)

    def test_data_sizes_must_match_2(self):
        with self.assertRaises(AssertionError):
            raster = np.arange(self._size + 10, dtype=np.float64)
            other = fourgp_speclib.Spectrum(wavelengths=raster,
                                            values=raster,
                                            value_errors=raster)
            self._polynomial.fit_to_continuum_via_template(other=other, template=self._absorption)

    def test_fitting(self):
        polynomial = self._polynomial.fit_to_continuum_via_template(other=self._observed, template=self._absorption)
        coefficients_expected = np.asarray([0, 1, 0], dtype=np.float64)
        self.assertLess(np.max(np.abs(np.asarray(polynomial.coefficients) - coefficients_expected)), 0.1)

    def tearDown(self):
        """
        Tear down Spectrum objects.
        """
        del self._observed, self._absorption


# Run tests if we are run from command line
if __name__ == '__main__':
    unittest.main()
