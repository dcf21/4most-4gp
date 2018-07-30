#!/usr/bin/env python2.7
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

    def onto_raster(self, output_raster, interpolate_errors=True, interpolate_mask=True):
        """
        Resample this spectrum onto a user-specified wavelength raster.

        :param output_raster:
            The raster we should resample this Spectrum onto.

        :param interpolate_errors:
            Should we bother interpolating the errors as well as the data itself? If not, the errors will be meaningless
            after the interpolation operation, but the function will return 30% quicker.

        :param interpolate_mask:
            Should we bother interpolating the spectrum's mask as the data itself? If not, the mask will be cleared, but
            the function will return 30% quicker.

        :return:
            New Spectrum object.
        """

        new_values = np.interp(x=output_raster,
                               xp=self._input.wavelengths,
                               fp=self._input.values)

        if interpolate_errors:
            new_value_errors = np.interp(x=output_raster,
                                         xp=self._input.wavelengths,
                                         fp=self._input.value_errors)
        else:
            new_value_errors = np.zeros_like(new_values)

        output = Spectrum(wavelengths=output_raster,
                          values=new_values,
                          value_errors=new_value_errors,
                          metadata=self._input.metadata.copy()
                          )

        if interpolate_mask and self._input.mask_set:
            output.mask = np.interp(x=output_raster,
                                    xp=self._input.wavelengths,
                                    fp=self._input.mask) > 0.5
            output.mask_set = not np.all(output.mask)

        return output

    def match_to_other_spectrum(self, other, interpolate_errors=True, interpolate_mask=True):
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

        return self.onto_raster(output_raster=other.wavelengths,
                                interpolate_errors=interpolate_errors,
                                interpolate_mask=interpolate_mask)
