#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module provides a class implementing a polynomial spectrum.
"""

import numpy as np
import logging

from spectrum import hash_numpy_array, Spectrum

logger = logging.getLogger(__name__)


class SpectrumPolynomial(Spectrum):
    """
    A class implementing a polynomial spectrum.
    """

    def __init__(self, wavelengths, terms=2, metadata=None):
        """
        
        :param wavelengths: 
            A 1D array listing the wavelengths at which this array of spectra are sampled.
            
        :type wavelengths:
            np.ndarray
            
        :param terms:
            The number of terms that this polynomial spectrum should have. A spectrum with zero terms still has a flat
            component; a polynomial with two terms is quadratic.
            
        :type terms:
            int
            
        :param metadata: 
            Dictionary of metadata about this spectrum.
            
        :type metadata:
            dict
        """
        self.terms = terms
        self.wavelengths = wavelengths

        super(SpectrumPolynomial, self).__init__(wavelengths=wavelengths,
                                                 values=self.values,
                                                 value_errors=self.value_errors,
                                                 metadata=metadata)

    @property
    def terms(self):
        return self._terms

    @terms.setter
    def terms(self, value):
        assert isinstance(value, int), \
            "Number of terms in polynomial spectrum must be an integer."

        assert value >= 0, \
            "A polynomial spectrum must have a positive number of terms."

        self._terms = value

        # Initialise polynomial coefficients to zero
        self.coefficients = [0] * self._terms

    @property
    def coefficients(self):
        return self._coefficients

    @coefficients.setter
    def coefficients(self, value):
        assert isinstance(value, (list, tuple)), \
            "Coefficients of polynomial spectrum must be a list of numerical terms."

        assert len(value) == self._terms + 1, \
            "Number of coefficients of polynomial spectrum must equal {}.".format(self._terms + 1)

        assert all(isinstance(i, (float, int)) for i in value), \
            "Polynominal spectrum coefficients must all be floats or ints."

        self._coefficients = value

    @property
    def wavelengths(self):
        return self._wavelengths

    @wavelengths.setter
    def wavelengths(self, value):
        assert isinstance(value, np.ndarray), \
            "The wavelength raster of a polynomial spectrum should be a numpy ndarray."

        self._wavelengths = value
        self._value_errors = np.zeros_like(self._wavelengths)
        self._update_raster_hash()
        self._update_values_array()

    @property
    def value_errors(self):
        return self._value_errors

    @value_errors.setter
    def value_errors(self, value):
        pass

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, value):
        pass

    def _update_values_array(self):
        """
        Evaluate polynomial at every point on wavelength raster.
        :return: 
            None
        """
        self._values = np.zeros_like(self._wavelengths)

        # Express polynomial as ((a_2 * x) + a_1)*x) + a_0
        for order in range(1, self._terms + 2):
            self._values *= self._wavelengths
            self._values += self._coefficients[-order]
