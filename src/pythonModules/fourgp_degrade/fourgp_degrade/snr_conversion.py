# -*- coding: utf-8 -*-

"""
A class which converts between SNR/pixel and SNR/A.
"""

from math import sqrt
import numpy as np


class SNRValue:
    """
    Class to hold an SNR value, both per pixel and per A
    """

    def __init__(self, value_per_a, value_per_pixel):
        """
        Class to hold an SNR value, both per pixel and per A

        :param value_per_a:
            Value of this SNR per A
        :param value_per_pixel:
            Value of this SNR per pixel
        """
        self.value_per_a = value_per_a
        self.value_per_pixel = value_per_pixel

    def per_pixel(self):
        """
        Return the SNR per pixel
        :return:
            float
        """
        return self.value_per_pixel

    def per_a(self):
        """
        Return the SNR per A
        :return:
            float
        """
        return self.value_per_a


class SNRConverter:
    """
    Class to convert SNRs between SNR/pixel and SNR/A
    """

    def __init__(self, raster_from_file=None, raster=None, snr_at_wavelength=6100):
        """
        Class to convert SNRs between SNR/pixel and SNR/A

        :param raster_from_file:
            The filename of a file from which we should read the wavelength raster being used. Optional.
        :param raster:
            A numpy array containing the wavelength raster being used. Optional.
        :param snr_at_wavelength:
            The wavelength at which we define SNR
        """
        self.snr_at_wavelength = snr_at_wavelength

        # Load wavelength raster
        if raster is not None:
            self.raster = raster
        elif raster_from_file is not None:
            self.raster = np.loadtxt(raster_from_file).transpose()[0]
        else:
            raise ValueError("No wavelength raster supplied")

        # Differentiate raster to get pixels per A
        raster_diff = np.diff(raster[raster > float(snr_at_wavelength)])
        self.pixels_per_angstrom = 1.0 / raster_diff[0]

    def per_pixel(self, value):
        """
        Create an SNR value with defined value per pixel.

        :param value:
            SNR/pixel
        :type value:
            float
        :return:
            SNRValue
        """
        return SNRValue(value_per_pixel=value,
                        value_per_a=value * sqrt(self.pixels_per_angstrom))

    def per_a(self, value):
        """
        Create an SNR value with defined value per A.

        :param value:
            SNR/A
        :type value:
            float
        :return:
            SNRValue
        """
        return SNRValue(value_per_pixel=value / sqrt(self.pixels_per_angstrom),
                        value_per_a=value)
