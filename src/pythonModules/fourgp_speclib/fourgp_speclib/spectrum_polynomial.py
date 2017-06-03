#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module provides a class implementing a polynomial spectrum.

This is useful for crudely modelling continuum. To produce a continuum-normalised version of an observed spectrum,
 first fit a polynomial to it, which is multiplied by a continuum-normalised absorption profile in the fitting process.
 The observed spectrum is then divided by the polynomial.
"""

import numpy as np
import logging
from scipy.optimize import least_squares

from spectrum import Spectrum

logger = logging.getLogger(__name__)


class SpectrumPolynomial(Spectrum):
    """
    A class implementing a polynomial spectrum.
    """

    def __init__(self, wavelengths, terms=2, coefficients=None, metadata=None):
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

        assert isinstance(terms, int), \
            "Number of terms in polynomial spectrum must be an integer."

        assert terms >= 0, \
            "A polynomial spectrum must have a positive number of terms."
        self._terms = terms

        if coefficients is None:
            coefficients = (0,) * (self._terms + 1)
        self._coefficients = coefficients

        assert isinstance(wavelengths, np.ndarray), \
            "The wavelength raster of a polynomial spectrum should be a numpy ndarray."
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
        self.coefficients = (0,) * (self._terms + 1)
        self._update_values_array()

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
        self._update_values_array()

    @property
    def wavelengths(self):
        return self._wavelengths

    @wavelengths.setter
    def wavelengths(self, value):
        assert isinstance(value, np.ndarray), \
            "The wavelength raster of a polynomial spectrum should be a numpy ndarray."

        self._wavelengths = np.asarray(value, dtype=np.float64)
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
        self._values = self._evaluate_polynomial(self._wavelengths, self.coefficients)

    def _evaluate_polynomial(self, raster, coefficients):
        """
        Evaluate polynomial at every point on wavelength raster.
        :return: 
            None
        """
        output = np.zeros_like(raster, dtype=np.float64)
        coefficients = np.asarray(coefficients, dtype=np.float64)

        # Express polynomial as ((a_2 * x) + a_1)*x) + a_0
        for order in range(1, self._terms + 2):
            output *= raster
            output += coefficients[-order]

        return output

    def fit_to_continuum(self, other, template, lambda_min_norm=None, lambda_max_norm=None):
        """
        Fit a polynomial to the continuum in the spectrum <other>. The spectrum <template> is an absorption line
        profile, such that we expect:
        
        other = polynomial * template
        
        :param other:
            The spectrum to fit this polynomial to.
            
        :type other:
            Spectrum
        
        :param template:
            The absorption line profile superimposed on the continuum in <other>. A value of 1 implies full transmission
            at a particular wavelength; a value of 0 implies complete absorption at this wavelength.
            
        :type template:
            Spectrum
            
        :param lambda_min_norm:
            Used to specify a wavelength range for an initial estimate of the continuum level. This is the short-ward
            end of the wavelength range.
            
        :type lambda_min_norm:
            float
            
        :param lambda_max_norm:
            Used to specify a wavelength range for an initial estimate of the continuum level. This is the long-ward
            end of the wavelength range.
            
        :type lambda_max_norm:
            float
            
        :return:
            None
        """

        assert isinstance(other, Spectrum), \
            "The fit_to method requires a Spectrum object to fit a polynomial to. Supplied object has type {}.". \
                format(type(other))

        assert isinstance(template, Spectrum), \
            "The fit_to method requires a Spectrum object to fit a polynomial to. Supplied object has type {}.". \
                format(type(other))

        assert other.raster_hash == template.raster_hash, \
            "The two spectra passed to the fit_to_continuum() method must be sampled on the same raster."

        if lambda_min_norm is None:
            lambda_min_norm = other.wavelengths[0]

        if lambda_max_norm is None:
            lambda_max_norm = other.wavelengths[-1]

        # Exclude any wavelengths that are masked out or not finite
        mask = (other.mask * template.mask *
                (other.value_errors > 0) * np.isfinite(other.values) * np.isfinite(template.values))

        other_values_masked = other.values[mask]
        other_value_errors_masked = other.value_errors[mask]
        other_wavelengths_masked = other.wavelengths[mask]
        template_values_masked = template.values[mask]

        assert len(other) > 0, "No good data in spectrum <other>."
        assert len(template) > 0, "No good data in spectrum <template>."

        def error_func(coefficients, other_wavelengths_masked_, other_values_masked_,
                       other_value_errors_masked_, template_values_masked_):
            return ((other_values_masked_ -
                     template_values_masked_ * self._evaluate_polynomial(other_wavelengths_masked_, coefficients)) /
                    other_value_errors_masked_)

        # roughly fit the template to the observed spectrum
        norm_factor = np.median(
            other.values[(other.wavelengths > lambda_min_norm) & (other.wavelengths < lambda_max_norm)])
        coefficients_initial = np.asarray((norm_factor,) + (0,) * self._terms)

        result = least_squares(error_func, coefficients_initial, args=(other_wavelengths_masked,
                                                                       other_values_masked,
                                                                       other_value_errors_masked,
                                                                       template_values_masked))

        self.coefficients = tuple(result.x)
