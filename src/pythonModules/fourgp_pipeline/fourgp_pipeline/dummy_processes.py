# -*- coding: utf-8 -*-

"""
This file contains dummy implementations for parts of the 4GP pipeline which aren't implemented yet.
"""

from fourgp_speclib import Spectrum


class ContinuumNormalisationDummy:
    """
    This is a dummy continuum normalisation routine, which simply returns the spectrum it is given.
    """

    def __init__(self):
        pass

    @staticmethod
    def continuum_normalise(spectrum_flux_normalised):
        """
        Dummy continuum normalisation routine, which simply returns the spectrum it is given.

        :param spectrum_flux_normalised:
            The flux-normalised spectrum which we are to continuum normalise.
        :type spectrum_flux_normalised:
            Spectrum
        :return:
            The continuum-normalised spectrum.
        """

        assert isinstance(spectrum_flux_normalised, Spectrum), \
            "The continuum normalisation routine has been requested to continuum normalise an object which isn't a " \
            "Spectrum object."

        return spectrum_flux_normalised
