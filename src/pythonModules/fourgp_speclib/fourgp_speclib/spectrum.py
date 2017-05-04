#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from os import path as os_path
import numpy as np
import hashlib
import logging

logger = logging.getLogger(__name__)


def hash_numpy_array(item):
    """
    Efficiently produce a string hash of a numpy array, for quick checking of whether arrays match.
    
    :param item:
        Any numpy array
        
    :type item:
        np.ndarray
        
    :return:
        String hash
    """
    raw = item.view(np.uint8)
    hashlib.sha1(raw).hexdigest()


def requires_common_raster(method):
    """
    A decorator for spectrum methods that require an other spectrum to be sampled on the same wavelength raster as us.

    :param method:
        A method belonging to a sub-class of Spectrum.
    """

    def wrapper(spectrum, other, *args, **kwargs):

        assert spectrum.raster_hash == other.raster_hash, \
            "Cannot do arithmetic on spectra sampled on a different wavelength rasters"

        return method(spectrum, other=other, *args, **kwargs)

    return wrapper


class Spectrum(object):
    """
    An object representing a single spectrum.
    
    :ivar np.ndarray wavelengths:
        A 1D array listing the wavelengths at which this array of spectra are sampled.
        
    :ivar np.ndarray fluxes:
        A 1D array listing the flux measurements in this spectrum.
        
    :ivar np.ndarray flux_errors:
        A 1D array listing the standard errors in the flux measurements.
        
    :ivar np.ndarray mask:
        A 1D array listing which wavelength samples we've currently selected to use.
        
    :ivar bool mask_set:
        Boolean selecting whether any wavelengths are currently masked out.
        
    :ivar str raster_hash:
        A string hash of the wavelength raster, used to quickly check whether spectra are sampled on a common raster.
    """

    def __init__(self, wavelengths, fluxes, flux_errors):
        """
        Instantiate a new Spectrum object.
        
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
        self.mask = np.ones_like(self.wavelengths)
        self.mask_set = False

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

    @requires_common_raster
    def copy_mask_from(self, other):
        """
        Duplicate the wavelength mask set up on another spectrum array, and apply it to this array.
        
        :param other:
            The spectrum array to copy the mask from
            
        :type other:
            SpectrumArray
        
        :return:
            None
        """
        self.mask = other.mask.copy()

    def mask_include(self, wavelength_min=0, wavelength_max=np.inf):
        """
        Update the wavelength mask to include wavelengths in the specified range.
        
        :param wavelength_min: 
            The shortest wavelength in the range to be added into the mask.
            
        :type wavelength_min:
            float
        
        :param wavelength_max:
            The longest wavelength in the range to be added into the mask.
            
        :type wavelength_max:
            float
        
        :return: 
            None
        """

        window = (self.wavelengths >= wavelength_min) * (self.wavelengths <= wavelength_max)
        self.mask[window] = True
        self.mask_set = not np.all(self.mask)

    def mask_exclude(self, wavelength_min=0, wavelength_max=np.inf):
        """
        Update the wavelength mask to exclude wavelengths in the specified range.
        
        :param wavelength_min: 
            The shortest wavelength in the range to be removed from the mask.
            
        :type wavelength_min:
            float
        
        :param wavelength_max:
            The longest wavelength in the range to be removed from the mask.
            
        :type wavelength_max:
            float
        
        :return: 
            None
        """

        window = (self.wavelengths >= wavelength_min) * (self.wavelengths <= wavelength_max)
        self.mask[window] = False
        self.mask_set = not np.all(self.mask)

    @requires_common_raster
    def __add__(self, other):
        """
        Add the fluxes in another spectrum to the fluxes in this one, and return a new Spectrum object.
        
        :param other:
            The Spectrum object to add to this one.
            
        :type other:
            Spectrum
        
        :return:
            Spectrum object containing the sum of the two spectra.
        """

        new_flux_errors = np.hypot(self.flux_errors, other.flux_errors)

        if not (self.mask_set or other.mask_set):
            new_flux = self.fluxes + other.fluxes
        else:
            new_flux = self.fluxes * self.mask + other.fluxes * other.mask
            new_flux_errors[~(self.mask * other.mask)] = np.inf

        return Spectrum(wavelengths=self.wavelengths, fluxes=new_flux, flux_errors=new_flux_errors)

    @requires_common_raster
    def __sub__(self, other):
        """
        Subtract the fluxes in another spectrum from the fluxes in this one, and return a new Spectrum object.
        
        :param other:
            The Spectrum object to subtract from this one.
            
        :type other:
            Spectrum
        
        :return:
            Spectrum object containing the difference of the two spectra.
        """

        new_flux_errors = np.hypot(self.flux_errors, other.flux_errors)

        if not (self.mask_set or other.mask_set):
            new_flux = self.fluxes - other.fluxes
        else:
            new_flux = self.fluxes * self.mask - other.fluxes * other.mask
            new_flux_errors[~(self.mask * other.mask)] = np.inf

        return Spectrum(wavelengths=self.wavelengths, fluxes=new_flux, flux_errors=new_flux_errors)

    @requires_common_raster
    def __iadd__(self, other):
        """
        Add the fluxes in another spectrum to the fluxes in this one.
        
        :param other:
            The Spectrum object to add to this one.
            
        :type other:
            Spectrum
        
        :return:
            self
        """

        self.flux_errors = np.hypot(self.flux_errors, other.flux_errors)

        if not (self.mask_set or other.mask_set):
            self.fluxes = self.fluxes + other.fluxes
        else:
            self.fluxes = self.fluxes * self.mask + other.fluxes * other.mask
            self.flux_errors[~(self.mask * other.mask)] = np.inf

        return self

    @requires_common_raster
    def __isub__(self, other):
        """
        Subtract the fluxes in another spectrum from the fluxes in this one.
        
        :param other:
            The Spectrum object to subtract from this one.
            
        :type other:
            Spectrum
        
        :return:
            self
        """

        self.flux_errors = np.hypot(self.flux_errors, other.flux_errors)

        if not (self.mask_set or other.mask_set):
            self.flux = self.fluxes - other.fluxes
        else:
            self.flux = self.fluxes * self.mask - other.fluxes * other.mask
            self.flux_errors[~(self.mask * other.mask)] = np.inf

        return self

    @requires_common_raster
    def __mul__(self, other):
        """
        Multiply the fluxes in another spectrum by the fluxes in this one, and return a new Spectrum object.
        
        :param other:
            The Spectrum object to multiply by this one.
            
        :type other:
            Spectrum
        
        :return:
            Spectrum object containing the sum of the two spectra.
        """

        if not (self.mask_set or other.mask_set):
            new_flux = self.fluxes * other.fluxes
        else:
            new_flux = self.fluxes * self.mask * other.fluxes * other.mask

        new_flux_errors = np.hypot(self.flux_errors / self.fluxes, other.flux_errors / other.fluxes) * np.abs(new_flux)

        if self.mask_set or other.mask_set:
            new_flux_errors[~(self.mask * other.mask)] = np.inf

        return Spectrum(wavelengths=self.wavelengths, fluxes=new_flux, flux_errors=new_flux_errors)

    @requires_common_raster
    def __div__(self, other):
        """
        Divide the fluxes in this spectrum by the fluxes in another one, and return a new Spectrum object.
        
        :param other:
            The Spectrum object to divide this one by.
            
        :type other:
            Spectrum
        
        :return:
            Spectrum object containing the quotient of the two spectra.
        """

        if not (self.mask_set or other.mask_set):
            new_flux = self.fluxes / other.fluxes
        else:
            new_flux = (self.fluxes * self.mask) / (other.fluxes * other.mask)

        new_flux_errors = np.hypot(self.flux_errors / self.fluxes, other.flux_errors / other.fluxes) * np.abs(new_flux)

        if self.mask_set or other.mask_set:
            new_flux_errors[~(self.mask * other.mask)] = np.inf

        return Spectrum(wavelengths=self.wavelengths, fluxes=new_flux, flux_errors=new_flux_errors)

    @requires_common_raster
    def __iadd__(self, other):
        """
        Add the fluxes in another spectrum to the fluxes in this one.
        
        :param other:
            The Spectrum object to add to this one.
            
        :type other:
            Spectrum
        
        :return:
            self
        """

        self.flux_errors = np.hypot(self.flux_errors, other.flux_errors)

        if not (self.mask_set or other.mask_set):
            self.fluxes = self.fluxes + other.fluxes
        else:
            self.fluxes = self.fluxes * self.mask + other.fluxes * other.mask
            self.flux_errors[~(self.mask * other.mask)] = np.inf

        return self

    @requires_common_raster
    def __isub__(self, other):
        """
        Subtract the fluxes in another spectrum from the fluxes in this one.
        
        :param other:
            The Spectrum object to subtract from this one.
            
        :type other:
            Spectrum
        
        :return:
            self
        """

        self.flux_errors = np.hypot(self.flux_errors, other.flux_errors)

        if not (self.mask_set or other.mask_set):
            self.flux = self.fluxes - other.fluxes
        else:
            self.flux = self.fluxes * self.mask - other.fluxes * other.mask
            self.flux_errors[~(self.mask * other.mask)] = np.inf

        return self
