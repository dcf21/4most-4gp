#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from fourgp_speclib import Spectrum


class SpectrumInterpolator(object):
    """
    A class containing utility functions for interpolating spectra onto different wavelength rasters
    """

    def __init__(self, input_spectrum):
        assert isinstance(input_spectrum, Spectrum), "The SpectrumInterpolate class can only operate on Spectrum objects."
        self._input = input_spectrum

    def match_to_other_spectrum(self, other, interpolate_errors=True, interpolate_mask = True):
        """
        Resample this spectrum onto the wavelength raster of another Spectrum object.

        :param other:
            The other Spectrum object whose raster we should resample this Spectrum onto.

        :param interpolate_errors:
            Should we bother interpolating the errors as well as the data itself? If not, the errors will be meaningless
            after the interpolation operation, but the function will return 30% quicker.

        :param interpolate_mask:
            Should we bother interpolating the spectrum's mask as the data itself? If not, the mask will be cleared, but
            the function will return 30% quicker.

        :return:
            New Spectrum object.
        """

        new_values = np.interp(x=other.wavelengths, xp=self._input.wavelengths, fp=self._input.values)

        if interpolate_errors:
            new_value_errors = np.interp(x=other.wavelengths, xp=self._input.wavelengths, fp=self._input.value_errors)
        else:
            new_value_errors = self._input.value_errors

        output = Spectrum(wavelengths=other.wavelengths,
                          values=new_values,
                          value_errors=new_value_errors)

        if interpolate_mask and self._input.mask_set:
            output.mask = np.interp(x=other.wavelengths, xp=self._input.wavelengths, fp=self._input.mask) > 0.5
            output.mask_set = not np.all(output.mask)

        return output
