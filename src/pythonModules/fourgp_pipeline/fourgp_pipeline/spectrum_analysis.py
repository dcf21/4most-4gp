# -*- coding: utf-8 -*-

"""

The `SpectrumAnalysis` class is created by the `Pipeline` class, and represents
the analysis of a particular spectrum.

"""


class SpectrumAnalysis:
    """
    Class which is created by the `Pipeline` class, and represents
    the analysis of a particular spectrum.
    """

    def __init__(self, input_spectrum):
        """
        Class which represents the analysis of a particular spectrum.

        :param input_spectrum:
            The Spectrum object we are to analyse.
        :type input_spectrum:
            Spectrum
        """

        self.input_spectrum = input_spectrum
        self.failure = False
        self.intermediate_results = []
