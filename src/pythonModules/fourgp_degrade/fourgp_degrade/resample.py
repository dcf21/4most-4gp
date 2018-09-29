# -*- coding: utf-8 -*-

import numpy as np

from fourgp_speclib import Spectrum


class SpectrumResampler(object):
    """
    A class containing utility functions for resampling spectra onto different wavelength rasters
    """

    def __init__(self, input_spectrum):
        assert isinstance(input_spectrum, Spectrum), "The SpectrumResampler class can only operate on Spectrum objects."
        self._input = input_spectrum

    @staticmethod
    def _pixel_left_edges(raster):
        """
        Turn a raster of the central wavelengths of pixels, into the left-hand edge of each pixel.

        For example:

        >>> raster
        array([0., 1., 2., 3., 4.])
        >>> left_edges
        array([-0.5,  0.5,  1.5,  2.5,  3.5,  4.5])

        :param raster:
            Numpy array containing the input raster of the central wavelengths of pixels
        :return:
            Numpy array containing the left-hand edge of each pixel
        """

        left_edges = (raster[1:] + raster[:-1]) / 2
        left_edges = np.insert(arr=left_edges,
                               obj=0,
                               values=raster[0] * 1.5 - raster[1] * 0.5
                               )
        left_edges = np.append(arr=left_edges,
                               values=raster[-1] * 1.5 - raster[-2] * 0.5
                               )
        return left_edges

    @staticmethod
    def _pixel_widths(left_edges):
        """
        Turn a raster of the central wavelengths of pixels, into the widths of each pixel.

        For example:

        >>> left_edges
        array([0., 1., 2., 3., 4.])
        >>> pixel_widths
        array([1, 1, 1, 1])

        :param left_edges:
            Numpy array containing the left-edges of the pixels (N+1 entries)
        :return:
            Numpy array containing the widths of each pixel (N entries)
        """

        pixel_widths = left_edges[1:] - left_edges[:-1]
        return pixel_widths

    @staticmethod
    def _resample(x_new, x_in, y_in):
        """
        Function which actually resamples a 1D spectrum onto a new raster of wavelengths.

        :param x_new:
            The new raster of wavelengths

        :type x_new:
            np.ndarray

        :param x_in:
            The existing raster of wavelengths

        :type x_in:
            np.ndarray

        :param y_in:
            The existing spectrum, on the existing raster

        :type y_in:
            np.ndarray

        :return:
            A numpy array containing the resampled spectrum on the new raster
        """

        # Make sure that input arrays are numpy objects
        x_in = np.asarray(x_in)
        y_in = np.asarray(y_in)
        x_new = np.asarray(x_new)

        # Make sure that input data is sensible
        assert x_in.ndim == 1, \
            "Input x array should have exactly one dimension. Passed array has {} dimensions".format(x_in.ndim)
        assert y_in.ndim == 1, \
            "Input y array should have exactly one dimension. Passed array has {} dimensions".format(y_in.ndim)
        assert x_in.shape[0] > 3, \
            "Input spectrum must have at least three pixels for resampling to produce sensible output"
        assert x_in.shape == y_in.shape, \
            "Input x and y vectors have different lengths"
        assert x_new.ndim == 1, \
            "New x array should have exactly one dimension. Passed array has {} dimensions".format(x_new.ndim)

        # Make an array of the left edge of each pixel in the input raster. The final entry is the right-edge of the
        # last pixel, so if we have N input pixels, we output N+1 left edges (length N+1)
        x_in_pixel_left_edges = SpectrumResampler._pixel_left_edges(x_in)

        # Compute the width of each pixel in the input raster (length N)
        x_in_pixel_width = SpectrumResampler._pixel_widths(x_in_pixel_left_edges)

        # Do the same for the output raster
        x_new_pixel_left_edges = SpectrumResampler._pixel_left_edges(x_new)
        x_new_pixel_width = SpectrumResampler._pixel_widths(x_new_pixel_left_edges)

        # Create an array of the integrated flux per A leftwards of wavelength i in <x_in_pixel_left_edges> (length N+1)
        x_in_integrated = np.cumsum(np.insert(y_in * x_in_pixel_width, 0, 0))

        # Create an array of the integrated flux within each pixel of the output array, divided by the pixel's width
        y_new = (np.interp(xp=x_in_pixel_left_edges,
                           fp=x_in_integrated,
                           x=x_new_pixel_left_edges[1:]) -
                 np.interp(xp=x_in_pixel_left_edges,
                           fp=x_in_integrated,
                           x=x_new_pixel_left_edges[:-1])
                 ) / x_new_pixel_width

        # Return output
        return y_new

    def onto_raster(self, output_raster, resample_errors=True, resample_mask=True):
        """
        Resample this spectrum onto a user-specified wavelength raster.

        :param output_raster:
            The raster we should resample this Spectrum onto.

        :param resample_errors:
            Should we bother resampling the errors as well as the data itself? If not, the errors will be meaningless
            after the resampling operation, but the function will return 30% quicker.

        :param resample_mask:
            Should we bother resampling the spectrum's mask as the data itself? If not, the mask will be cleared, but
            the function will return 30% quicker.

        :return:
            New Spectrum object.
        """

        new_values = self._resample(x_new=output_raster,
                                    x_in=self._input.wavelengths,
                                    y_in=self._input.values)

        if resample_errors:
            new_value_errors = self._resample(x_new=output_raster,
                                              x_in=self._input.wavelengths,
                                              y_in=self._input.value_errors)
        else:
            new_value_errors = np.zeros_like(new_values)

        output = Spectrum(wavelengths=output_raster,
                          values=new_values,
                          value_errors=new_value_errors,
                          metadata=self._input.metadata.copy()
                          )

        if resample_mask and self._input.mask_set:
            output.mask = self._resample(x_new=output_raster,
                                         x_in=self._input.wavelengths,
                                         y_in=self._input.mask) > 0.5
            output.mask_set = not np.all(output.mask)

        return output

    def match_to_other_spectrum(self, other, resample_errors=True, resample_mask=True):
        """
        Resample this spectrum onto the wavelength raster of another Spectrum object.

        :param other:
            The other Spectrum object whose raster we should resample this Spectrum onto.

        :param resample_errors:
            Should we bother resampling the errors as well as the data itself? If not, the errors will be meaningless
            after the resampling operation, but the function will return 30% quicker.

        :param resample_mask:
            Should we bother resampling the spectrum's mask as the data itself? If not, the mask will be cleared, but
            the function will return 30% quicker.

        :return:
            New Spectrum object.
        """

        return self.onto_raster(output_raster=other.wavelengths,
                                resample_errors=resample_errors,
                                resample_mask=resample_mask)
