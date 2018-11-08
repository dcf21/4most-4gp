# -*- coding: utf-8 -*-

"""
This module implements a simple cross-correlation code for determining RVs.

It is a heavily cleaned up version of Jane Lin's GUESS code, as used by GALAH. Original code taken from
<https://github.com/jlin0504/GUESS>.
"""

from math import sqrt

import numpy as np
from scipy.optimize import leastsq
import logging

import fourgp_speclib
from fourgp_degrade.resample import SpectrumResampler

from .templates_resample import logarithmic_raster

logger = logging.getLogger(__name__)


class RvInstanceCrossCorrelation(object):
    """
    A class which is adapted from Jane Lin's GUESS code, as used by GALAH.
    """

    def __init__(self, spectrum_library):
        """
        Instantiate the RV code, and read from disk the library of template spectra used for cross correlation.

        :param spectrum_library:
            A SpectrumLibrary containing the template spectra we use for modelling.

        :type spectrum_library:
            SpectrumLibrary
        """

        assert isinstance(spectrum_library, fourgp_speclib.SpectrumLibrary), \
            "Argument to RvInstanceCrossCorrelation should be a SpectrumLibrary."
        self._spectrum_library = spectrum_library

        # Load template spectra
        spectrum_list = self._spectrum_library.search(continuum_normalised=1)

        # Make a list of templates by 4MOST wavelength arm
        self.templates_by_arm = {}
        self.arm_rasters = {}
        self.arm_properties = {}
        self.window_functions = {}

        # Sort template spectra by the arm they are sampled on
        for template_id in [i["specId"] for i in spectrum_list]:
            template_metadata = self._spectrum_library.get_metadata(ids=(template_id,))[0]
            mode = template_metadata['mode']
            arm_name = template_metadata['arm_name']

            if mode not in self.templates_by_arm:
                self.templates_by_arm[mode] = {}
            if arm_name not in self.templates_by_arm[mode]:
                self.templates_by_arm[mode][arm_name] = []

            # Add this template to the list of available templates for this wavelength arm
            self.templates_by_arm[mode][arm_name].append(template_id)

            # If we haven't already recreated the fixed-log-step wavelength raster for this arm, do it now
            if arm_name not in self.arm_rasters:
                self.arm_rasters[arm_name] = logarithmic_raster(lambda_min=template_metadata['lambda_min'],
                                                                lambda_max=template_metadata['lambda_max'],
                                                                lambda_step=template_metadata['lambda_step'])

                self.arm_properties[arm_name] = {
                    'lambda_min': template_metadata['lambda_min'],
                    'lambda_max': template_metadata['lambda_max'],
                    'lambda_step': template_metadata['lambda_step'],
                    'multiplicative_step': self.arm_rasters[arm_name][1] / self.arm_rasters[arm_name][0]
                }

                self.window_functions[arm_name] = self.window_function(template_length=len(self.arm_rasters[arm_name]))
                print("Arm {} requested length {}; got length {}".format(arm_name, len(self.arm_rasters[arm_name]),
                                                                         self.window_functions[arm_name].shape))

        # Load library of template spectra for each arm
        self.template_spectra = {}

        for mode in self.templates_by_arm:
            for arm_name in self.templates_by_arm[mode]:
                self.template_spectra[arm_name] = self._spectrum_library.open(
                    ids=self.templates_by_arm[mode][arm_name],
                    shared_memory=True
                )

        # Multiply template spectra by window function
        for arm_name in self.template_spectra:
            for index in range(len(self.template_spectra[arm_name])):
                template_spectrum = self.template_spectra[arm_name].extract_item(index=index)
                print("Shape of template: {}".format(template_spectrum.values.shape))
                print("Shape of window: {}".format(self.window_functions[arm_name].shape))
                template_spectrum.values *= self.window_functions[arm_name]

    @staticmethod
    def window_function(template_length):
        """
        Create a window function that we multiply each template by to avoid edge-effects. This goes to zero at either
        end.

        :param template_length:
            The number of pixels in the template spectrum we want to create a window function for.
        :return:
            np.array
        """

        # Create a linear ramp covering 10% of the length of the spectrum
        left_ramp = np.arange(int(template_length / 10))
        left_ramp = left_ramp / max(left_ramp)

        # Ramp to use at right-hand end of spectrum is simply a back-to-front version of the left ramp
        right_ramp = left_ramp[::-1]

        middle_length = template_length - len(left_ramp) - len(right_ramp)

        window_function = np.hstack([left_ramp,
                                     np.ones(middle_length),
                                     right_ramp
                                     ])

        return window_function

    def resample_single_arm(self, input_spectrum, arm_name):
        """
        Resample an input spectrum onto a raster with fixed logarithmic stride, representing a single 4MOST arm. We
        use the same wavelength rasters that the template spectra were sampled onto.

        :param input_spectrum:
            A Spectrum object, containing an observed spectrum
        :param arm_name:
            The name of the arm within this 4MOST mode
        :return:
            A Spectrum object containing a single arm of 4MOST data, resampled with a fixed logarithmic stride.
        """

        new_raster = self.arm_rasters[arm_name]

        resampler = SpectrumResampler(input_spectrum=input_spectrum)

        resampled_spectrum = resampler.onto_raster(output_raster=new_raster, resample_errors=True, resample_mask=False)

        return resampled_spectrum

    def estimate_rv_from_single_arm(self, input_spectrum, mode, arm_name):
        """
        Estimate the RV of a spectrum on the basis of data from a single arm. We return a list of RV estimates from
        cross correlation with each of the template spectra, and the chi-squared mismatch of each template spectrum
        which should be used to weight the RV estimates.

        :param input_spectrum:
            A Spectrum object, containing an observed spectrum. This should represent one arm of data only, and should
            be sampled onto the same fixed logarithmic wavelength stride as the template spectra.
        :param mode:
            The name of the 4MOST mode this arm is part of -- either LRS or HRS
        :param arm_name:
            The name of the arm within this 4MOST mode
        :return:
            List of [RV value, weight]
        """

        rv_fits = []

        # Loop over all template spectra
        for template_index, template_id in enumerate(self.templates_by_arm[mode][arm_name]):
            template_spectrum = self.template_spectra[arm_name].extract_item(template_index)

            input_array = input_spectrum.values * self.window_functions[arm_name]

            cross_correlation = np.correlate(a=template_spectrum.values,
                                             v=input_array,
                                             mode='same')

            # Find the index of the maximum of cross correlation function
            max_position = np.where(cross_correlation == max(cross_correlation))[0][0]

            # Make sure we don't go off the end of the array
            if max_position < 1:
                max_position = 1
            if max_position > len(input_array) - 2:
                max_position = len(input_array) - 2

            # Now make three points which straddle the maximum
            x_vals = np.array([max_position - 1, max_position, max_position + 1])
            y_vals = cross_correlation[x_vals]

            # Fit a quadratic curve through three points
            def quadratic_curve(p, x):
                return p[0] * x ** 2 + p[1] * x + p[2]

            def quadratic_peak_x(p):
                return -p[1] / (2 * p[0])

            def quadratic_mismatch(p, x, y):
                return quadratic_curve(p, x) - y

            p0, p1, p2 = leastsq(func=quadratic_mismatch,
                                 x0=np.array([0, 0, 0]),
                                 args=(x_vals, y_vals)
                                 )[0]

            peak_x = quadratic_peak_x(p=(p0, p1, p2))

            pixel_shift = peak_x - len(input_array) // 2

            multiplicative_shift = pixel_shift * self.arm_properties[arm_name]['multiplicative_step']

            c = 299792458.0
            velocity = c * (pow(multiplicative_shift, 2) - 1) / (pow(multiplicative_shift, 2) + 1)
            weight = max(cross_correlation)

            rv_fits.append((velocity, weight))

        return rv_fits

    def estimate_rv(self, input_spectrum, mode):
        """
        Estimate the RV of a spectrum on the basis of all of the 4MOST arms of either HRS or LRS.

        :param input_spectrum:
            A Spectrum object, containing an observed spectrum
        :param mode:
            The name of the 4MOST mode this arm is part of -- either LRS or HRS
        :return:
            [RV, error in RV]
        """

        rv_estimates = []

        arm_names = self.templates_by_arm[mode].keys()

        for arm_name in arm_names:
            new_rv_estimates = self.estimate_rv_from_single_arm(
                input_spectrum=self.resample_single_arm(input_spectrum=input_spectrum, arm_name=arm_name),
                mode=mode,
                arm_name=arm_name
            )
            rv_estimates.extend(new_rv_estimates)

        rv_mean = (sum([rv * weight for rv, weight in rv_estimates]) /
                   sum([weight for rv, weight in rv_estimates]))

        rv_variance = (sum([pow(rv - rv_mean, 2) * weight for rv, weight in rv_estimates]) /
                       sum([weight for rv, weight in rv_estimates]))

        rv_std_dev = sqrt(rv_variance)

        return rv_mean, rv_std_dev
