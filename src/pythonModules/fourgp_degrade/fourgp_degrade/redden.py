# -*- coding: utf-8 -*-

import numpy as np

from fourgp_speclib import Spectrum


class SpectrumReddener(object):
    """
    A class containing utility functions for reddening Spectrum objects

    Cardelli et al., ApJ 345, 245 (1989)
    Howarth, MNRAS 203, 301-304 (1983) -- for x < 0.3 (see below)
    """

    def __init__(self, input_spectrum):
        """
        :param input_spectrum:
            Spectrum object. Wavelength raster must be specified in A. Flux and flux errors can be in arbitrary units.

        :type input_spectrum:
            Spectrum
        """
        assert isinstance(input_spectrum, Spectrum), "The SpectrumReddener class can only operate on Spectrum objects."
        self._input = input_spectrum

    def redden(self, e_bv, r=3.1):
        return self.deredden(-e_bv, r)

    def deredden(self, e_bv, r=3.1):
        """
        Redden this spectrum.

        :param e_bv:
            E(B-V). Positive values deredden a spectrum. Negative values increase the reddening of a spectrum.
        :type e_bv:
            float
        :param r:
            A(V)/E(B-V). Typically assumed to be 3.1 for a standard dust grain size spectrum.
        :type r:
            float

        :return:
            New Spectrum object.
        """

        x = 1e4 / self._input.wavelengths
        extinction = np.zeros_like(self._input.wavelengths)

        # Regime 1
        mask = (x >= 8)
        h = x - 8
        a = -1.073 - 0.628 * h + 0.137 * h * h - 0.070 * h * h * h
        b = 13.670 + 4.257 * h - 0.420 * h * h + 0.374 * h * h * h
        extinction[mask] = r * a[mask] + b[mask]

        # Regime 2
        mask = (x >= 5.9) * (x < 8)
        h = x - 5.9
        fa = -0.04473 * h * h - 0.009779 * h * h * h
        fb = 0.2130 * h * h + 0.1207 * h * h * h
        a = 1.752 - 0.316 * x - 0.104 / ((x - 4.67) ** 2 + 0.341) + fa
        b = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) ** 2 + 0.263) + fb
        extinction[mask] = r * a[mask] + b[mask]

        # Regime 3
        mask = (x >= 3.3) * (x < 5.9)
        a = 1.752 - 0.316 * x - 0.104 / ((x - 4.67) ** 2 + 0.341)
        b = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) ** 2 + 0.263)
        extinction[mask] = r * a[mask] + b[mask]

        # Regime 4
        mask = (x >= 1.1) * (x < 3.3)
        y = x - 1.82
        a = (1
             + 0.17699 * y
             - 0.50447 * y * y
             - 0.02427 * y ** 3
             + 0.72085 * y ** 4
             + 0.01979 * y ** 5
             - 0.77530 * y ** 6
             + 0.32999 * y ** 7)
        b = (1.41338 * y
             + 2.28305 * y * y
             + 1.07233 * y ** 3
             - 5.38434 * y ** 4
             - 0.62251 * y ** 5
             + 5.30260 * y ** 6
             - 2.09002 * y ** 7)
        extinction[mask] = r * a[mask] + b[mask]

        # Regime 5
        mask = (x >= 0.2) * (x < 1.1)
        a = 0.574 * x ** 1.61
        b = -0.527 * x ** 1.61
        extinction[mask] = r * a[mask] + b[mask]

        # Regime 6
        mask = (x < 0.2)

        # This is the far-IR range. We use Howarth (1983), but normalized
        # to Cardelli et al., 1989 at x=0.2
        # 0.05056 is A_lambda/E(B-V) of Howarth, 1983 at x=0.2
        xx = 0.2 ** 1.61
        hilfy = (r * 0.574 - 0.527) * xx / 0.05056
        extinction_6 = hilfy * x * ((1.86 - 0.48 * x) * x - 0.1)
        extinction[mask] = extinction_6[mask]

        multiplier = 10 ** (0.4 * (extinction * e_bv))

        output = Spectrum(wavelengths=self._input.wavelengths,
                          values=self._input.values * multiplier,
                          value_errors=self._input.value_errors * multiplier,
                          metadata=self._input.metadata.copy()
                          )

        if self._input.mask_set:
            output.copy_mask_from(self._input)

        return output
