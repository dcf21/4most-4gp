# -*- coding: utf-8 -*-

from scipy.ndimage.filters import gaussian_filter1d

from fourgp_speclib import Spectrum


class SpectrumConvolver(object):
    """
    A class containing utility functions for convolving Spectrum objects with point spread functions
    """

    def __init__(self, input_spectrum):
        assert isinstance(input_spectrum, Spectrum), "The SpectrumConvolver class can only operate on Spectrum objects."
        self._input = input_spectrum

    def gaussian_convolve(self, sigma):
        """
        Convolve this spectrum with a Gaussian PSF.

        :param sigma:
            Standard deviation of point spread function in pixels.

        :return:
            New Spectrum object.
        """

        new_values = gaussian_filter1d(input=self._input.values, sigma=sigma)

        output = Spectrum(wavelengths=self._input.wavelengths,
                          values=new_values,
                          value_errors=self._input.value_errors,
                          metadata=self._input.metadata.copy()
                          )

        if self._input.mask_set:
            output.copy_mask_from(self._input)

        return output
