# -*- coding: utf-8 -*-

"""
This file contains dummy implementations for parts of the 4GP pipeline which aren't implemented yet.
"""

from fourgp_speclib import Spectrum
from .spectrum_analysis import SpectrumAnalysis


class ContinuumNormalisationDummy:
    """
    This is a dummy continuum normalisation routine, which simply returns the spectrum it is given.
    """

    def __init__(self):
        pass

    @staticmethod
    def continuum_normalise(spectrum_flux_normalised, spectrum_analysis):
        """
        Dummy continuum normalisation routine, which simply returns the spectrum it is given.

        :param spectrum_flux_normalised:
            The flux-normalised spectrum which we are to continuum normalise.
        :type spectrum_flux_normalised:
            Spectrum
        :param spectrum_analysis:
            Structure containing the results of the analysis of the spectrum so far.
        :type spectrum_analysis:
            SpectrumAnalysis
        :return:
            The continuum-normalised spectrum.
        """

        assert isinstance(spectrum_flux_normalised, Spectrum), \
            "The continuum normalisation routine has been requested to continuum normalise an object which isn't a " \
            "Spectrum object."

        assert isinstance(spectrum_analysis, SpectrumAnalysis), \
            "The spectrum analysis data structure we have been passed is of the incorrect type."

        return spectrum_flux_normalised


class DecisionMakerDummy:
    """
    This is a dummy decision maker, which simply checks whether there was an error.
    """

    def __init__(self):
        pass

    @staticmethod
    def consult(spectrum_analysis):
        """
        Dummy dummy decision maker, which simply checks whether there was an error.

        :param spectrum_analysis:
            Structure containing the results of the analysis of the spectrum so far.
        :type spectrum_analysis:
            SpectrumAnalysis
        :return:
            Boolean flag indicating whether analysis was good.
        """

        return not spectrum_analysis.failure
