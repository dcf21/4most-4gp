#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import logging

from spectrum import hash_numpy_array, Spectrum

logger = logging.getLogger(__name__)


class SpectrumArray(object):
    """
    An object representing an array of spectra.
    
    :ivar np.ndarray wavelengths:
        A 1D array listing the wavelengths at which this array of spectra are sampled.
        
    :ivar np.ndarray fluxes:
        A 2D array listing the flux measurements for each spectrum in this SpectrumArray.
        
    :ivar np.ndarray flux_errors:
        A 2D array listing the standard errors in the flux measurements for each spectrum in this SpectrumArray.
        
    :ivar str raster_hash:
        A string hash of the wavelength raster, used to quickly check whether spectra are sampled on a common raster.
    """

    def __init__(self, wavelengths, fluxes, flux_errors):
        """
        Instantiate new SpectrumArray object.
        
        :param wavelengths: 
            A 1D array listing the wavelengths at which this array of spectra are sampled.
        
        :param fluxes: 
            A 2D array listing the flux measurements for each spectrum in this SpectrumArray.
            
        :param flux_errors: 
            A 2D array listing the standard errors in the flux measurements for each spectrum in this SpectrumArray.
        """
        self.wavelengths = wavelengths
        self.fluxes = fluxes
        self.flux_errors = flux_errors
        self._update_raster_hash()

    def __str__(self):
        return "<{module}.{name} instance".format(module=self.__module__,
                                                  name=type(self).__name__)

    def __repr__(self):
        return "<{0}.{1} object at {2}>".format(self.__module__,
                                                type(self).__name__, hex(id(self)))

    def _update_raster_hash(self):
        """
        Update the internal string hash of the wavelength raster that this spectrum array is sampled on.
        
        This hash is used to quickly check whether two spectra are sampled on the same raster before doing arithmetic
        operations on them.
        
        :return:
            None
        """
        self.raster_hash = hash_numpy_array(self.wavelengths)

    def extract_item(self, index):
        """
        Extract a single spectrum from a SpectrumArray. This creates a numpy view of the spectrum, without copying the
        data.
        
        :param index:
            Index of the spectrum to extract
            
        :type index:
            int
            
        :return:
            Spectrum object
        """

        assert 0 <= index < self.fluxes.shape[0], "Index of SpectrumArray out of range."

        return Spectrum(wavelengths=self.wavelengths,
                        fluxes=self.fluxes[index, :],
                        flux_errors=self.flux_errors[index, :])