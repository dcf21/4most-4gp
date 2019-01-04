# -*- coding: utf-8 -*-

"""
This module implements a simple cross-correlation code for determining RVs.

It is a heavily cleaned up version of Jane Lin's GUESS code, as used by GALAH. Original code taken from
<https://github.com/jlin0504/GUESS>.
"""

import logging
from math import sqrt
from operator import itemgetter

import fourgp_speclib
import numpy as np
from fourgp_degrade.resample import SpectrumResampler
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import leastsq

from .templates_resample import logarithmic_raster

logger = logging.getLogger(__name__)


class RvInstanceCrossCorrelation(object):
    """
    A class which is adapted from Jane Lin's GUESS code, as used by GALAH.
    """

    def __init__(self, spectrum_library, upsampling=1):
        """
        Instantiate the RV code, and read from disk the library of template spectra used for cross correlation.

        :param spectrum_library:
            A SpectrumLibrary containing the template spectra we use for modelling.
        :param upsampling:
            The factor by which to up-sample the spectrum before doing cross-correlation. A value of 1 means we don't
            up sample.
        :type spectrum_library:
            SpectrumLibrary
        """

        assert isinstance(spectrum_library, fourgp_speclib.SpectrumLibrary), \
            "Argument to RvInstanceCrossCorrelation should be a SpectrumLibrary."

        self._spectrum_library = spectrum_library
        self.upsampling = upsampling

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
                self.arm_properties[arm_name] = {
                    'lambda_min': template_metadata['lambda_min'],
                    'lambda_max': template_metadata['lambda_max'],
                    'lambda_step': template_metadata['lambda_step'],
                }

        # Load library of template spectra for each arm
        self.template_spectra_raw = {}

        for mode in self.templates_by_arm:
            for arm_name in self.templates_by_arm[mode]:
                self.template_spectra_raw[arm_name] = self._spectrum_library.open(
                    ids=self.templates_by_arm[mode][arm_name],
                    shared_memory=True
                )

                self.arm_rasters[arm_name] = self.template_spectra_raw[arm_name].wavelengths

                self.arm_properties[arm_name]['multiplicative_step'] = (self.arm_rasters[arm_name][1] /
                                                                        self.arm_rasters[arm_name][0])

                window_function_length = len(self.arm_rasters[arm_name])
                if self.upsampling > 1:
                    window_function_length = (len(self.arm_rasters[arm_name]) - 1) * upsampling

                self.window_functions[arm_name] = self.window_function(
                    template_length=window_function_length
                )

        # Multiply template spectra by window function and normalise
        self.template_spectra_tapered = {}

        for arm_name in self.template_spectra_raw:
            new_template_list = []

            for index in range(len(self.template_spectra_raw[arm_name])):
                template_spectrum = self.template_spectra_raw[arm_name].extract_item(index=index)

                if self.upsampling > 1:
                    template_spectrum = self.upsample_spectrum(
                        input=template_spectrum,
                        upsampling_factor=self.upsampling)

                # Multiply template by window function
                template_spectrum.values *= self.window_functions[arm_name]

                template_sum = np.sum(template_spectrum.values)

                # Renormalise template
                template_spectrum.values *= 1. / template_sum

                # Add modified template to list of updated templates for this arm
                new_template_list.append(template_spectrum)

            # Create new spectrum array of modified template spectra
            self.template_spectra_tapered[arm_name] = fourgp_speclib.SpectrumArray.from_spectra(
                spectra=new_template_list,
                shared_memory=True
            )

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
        use the same wavelength rasters that the template spectra were sampled onto. We also renormalise the resampled
        spectrum so that all the pixels sum to one.

        :param input_spectrum:
            A Spectrum object, containing an observed spectrum
        :param arm_name:
            The name of the arm within this 4MOST mode
        :return:
            A Spectrum object containing a single arm of 4MOST data, resampled with a fixed logarithmic stride.
        """

        # Look up the raster we are to resample onto
        new_raster = self.arm_rasters[arm_name]

        # Resample the input spectrum onto new raster
        resampler = SpectrumResampler(input_spectrum=input_spectrum)

        resampled_spectrum = resampler.onto_raster(output_raster=new_raster, resample_errors=True, resample_mask=False)

        # Renormalise the output spectrum
        resampled_spectrum_sum = np.sum(resampled_spectrum.values)

        resampled_spectrum.values *= 1. / resampled_spectrum_sum

        return resampled_spectrum

    def estimate_rv_from_single_arm(self, input_spectrum, mode, arm_name, interpolation_scheme="quadratic",
                                    interpolation_pixels=3):
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
        :param interpolation_scheme:
            The type of function to use to interpolate the CCF to measure sub-pixel RVs.
        :param interpolation_pixels:
            The number of pixels around the peak of the CCF to use when interpolating to measure sub-pixel RVs.
        :return:
            List of [RV value, weight]
        """

        assert interpolation_scheme in self.supported_interpolation_schemes()

        rv_fits = []

        if self.upsampling > 1:
            input_spectrum = self.upsample_spectrum(input=input_spectrum, upsampling_factor=self.upsampling)

        # Loop over all template spectra
        for template_index, template_id in enumerate(self.templates_by_arm[mode][arm_name]):
            template_spectrum = self.template_spectra_tapered[arm_name].extract_item(index=template_index)

            # Look up stellar parameters of this template
            template_metadata = self.template_spectra_tapered[arm_name].get_metadata(index=template_index)
            teff = template_metadata['Teff']
            logg = template_metadata['logg']
            fe_h = template_metadata['[Fe/H]']

            input_array = input_spectrum.values * self.window_functions[arm_name]

            cross_correlation = np.correlate(a=template_spectrum.values,
                                             v=input_array,
                                             mode='same')

            # Find the index of the maximum of cross correlation function
            try:
                max_position = np.where(cross_correlation == max(cross_correlation))[0][0]
            except IndexError:
                print("Warning: Cross-correlation failed")
                max_position = np.nan

            if np.isfinite(max_position):
                # Make sure we don't go off the end of the array
                if max_position < (interpolation_pixels // 2):
                    max_position = (interpolation_pixels // 2)
                if max_position > len(input_array) - 1 - (interpolation_pixels // 2):
                    max_position = len(input_array) - 1 - (interpolation_pixels // 2)

                # Now make three points which straddle the maximum
                x_min = max_position - interpolation_pixels // 2
                x_max = x_min + interpolation_pixels - 1
                x_vals = np.array(range(x_min, x_max + 1))
                y_vals = cross_correlation[x_vals]

                # Put peak close to zero for numerical stability
                x_vals = x_vals - max_position
                y_peak = y_vals[1]
                y_vals = y_vals - y_peak

                # Do interpolation
                if interpolation_scheme == "quadratic":
                    if len(x_vals) == 3:
                        if template_index == 0:
                            logger.info("Using analytic quadratic interpolation")
                        peak_x = self.interpolation_quadratic_dcf(x_vals=x_vals, y_vals=y_vals)
                    else:
                        if template_index == 0:
                            logger.info("Using numerical quadratic interpolation")
                        peak_x = self.interpolation_quadratic_galah(x_vals=x_vals, y_vals=y_vals)
                else:
                    if template_index == 0:
                        logger.info("Using spline interpolation")
                    peak_x = self.interpolation_spline(x_vals=x_vals, y_vals=y_vals)
            else:
                peak_x = np.nan

            # Shift peak back from zero to original position
            peak_x = peak_x + max_position

            # Convert position of peak of correlation function into a shift in pixels
            pixel_shift = peak_x - len(input_array) // 2

            # Convert pixel shift into multiplicative change in wavelength, based on fixed logarithmic stride
            multiplicative_shift = pow(self.arm_properties[arm_name]['multiplicative_step'],
                                       pixel_shift / self.upsampling)

            # Convert multiplicative wavelength shift into a radial velocity
            c = 299792458.0
            velocity = -c * (pow(multiplicative_shift, 2) - 1) / (pow(multiplicative_shift, 2) + 1)
            weight = max(cross_correlation)

            rv_fits.append(
                (velocity, weight, (teff, logg, fe_h))
            )

        return rv_fits

    @staticmethod
    def interpolation_quadratic_galah(x_vals, y_vals):
        """
        Use quadratic interpolation to find the sub-pixel position of the peak of the CCF. This is the GALAH way
        of doing it, which uses numerical fitting rather than an analytic solution

        :param x_vals:
            If the CCF is y(x), this is the array of the x values of the data points supplied to interpolate.
        :param y_vals:
            If the CCF is y(x), this is the array of the y values of the data points supplied to interpolate.
        :return:
            Our best estimate of the position of the peak.
        """

        # Fit a quadratic curve through three points
        def quadratic_curve(p, x):
            return p[0] * x ** 2 + p[1] * x + p[2]

        def quadratic_peak_x(p):
            return -p[1] / (2 * p[0])

        def quadratic_mismatch(p, x, y):
            return quadratic_curve(p, x) - y

        p0, p1, p2 = leastsq(func=quadratic_mismatch,
                             x0=np.array([0, 0, 0]),
                             args=(x_vals, y_vals),
                             maxfev=10000
                             )[0]

        peak_x = quadratic_peak_x(p=(p0, p1, p2))
        return peak_x

    @staticmethod
    def interpolation_quadratic_dcf(x_vals, y_vals):
        """
        Use quadratic interpolation to find the sub-pixel position of the peak of the CCF. This is Dominic's version
        which uses analytic solution

        :param x_vals:
            If the CCF is y(x), this is the array of the x values of the data points supplied to interpolate.
        :param y_vals:
            If the CCF is y(x), this is the array of the y values of the data points supplied to interpolate.
        :return:
            Our best estimate of the position of the peak.
        """

        assert len(x_vals) == 3  # This analytic solution only works with three input points
        assert x_vals[1] == 0  # Three points must be centred around x==0
        assert x_vals[2] == - x_vals[0]

        def quadratic_peak_x(p):
            return -p[1] / (2 * p[0])

        p2 = y_vals[1]

        p0 = (y_vals[0] + y_vals[2]) / (2 * x_vals[0] ** 2)

        p1 = (y_vals[0] - y_vals[2]) / (2 * x_vals[0])

        peak_x = quadratic_peak_x(p=(p0, p1, p2))

        return peak_x

    @staticmethod
    def interpolation_spline(x_vals, y_vals):
        """
        Use cubic spline interpolation to find the sub-pixel position of the peak of the CCF

        :param x_vals:
            If the CCF is y(x), this is the array of the x values of the data points supplied to interpolate.
        :param y_vals:
            If the CCF is y(x), this is the array of the y values of the data points supplied to interpolate.
        :return:
            Our best estimate of the position of the peak.
        """

        def quadratic_spline_roots(spl):
            roots = []
            knots = spl.get_knots()
            for a, b in zip(knots[:-1], knots[1:]):
                u, v, w = spl(a), spl((a + b) / 2), spl(b)
                t = np.roots([u + w - 2 * v, w - u, 2 * v])
                t = t[np.isreal(t) & (np.abs(t) <= 1)]
                roots.extend(t * (b - a) / 2 + (b + a) / 2)
            return np.array(roots)

        f = InterpolatedUnivariateSpline(x=x_vals, y=y_vals)
        cr_pts = quadratic_spline_roots(f.derivative())
        cr_vals = f(cr_pts)
        if len(cr_vals) == 0:
            return np.nan
        max_index = np.argmax(cr_vals)
        return cr_pts[max_index]

    @staticmethod
    def supported_interpolation_schemes():
        return "quadratic", "spline"

    def upsample_spectrum(self, input, upsampling_factor):
        """
        Upsample a spectrum object using cubic spline interpolation
        :param input:
            The Spectrum object we should up sample
        :param upsampling_factor:
            The integer factor by which to up-sample the spectrum
        :return:
            An up-sampled Spectrum object
        """

        multiplicative_spacing_in = input.wavelengths[1] / input.wavelengths[0]
        multiplicative_spacing_out = pow(multiplicative_spacing_in, 1. / upsampling_factor)

        # We impose an explicit length on the output, because the arange() here is numerically unstable about whether
        # it includes the final point or not
        raster_in_length = len(input.wavelengths)
        raster_out_length = (raster_in_length - 1) * upsampling_factor

        raster_out = logarithmic_raster(lambda_min=input.wavelengths[0],
                                        lambda_max=input.wavelengths[-1],
                                        lambda_step=input.wavelengths[0] * (multiplicative_spacing_out - 1)
                                        )[:raster_out_length]

        f = InterpolatedUnivariateSpline(x=input.wavelengths, y=input.values)

        return fourgp_speclib.Spectrum(wavelengths=raster_out,
                                       values=f(raster_out),
                                       value_errors=np.zeros_like(raster_out),
                                       metadata=input.metadata
                                       )

    def estimate_rv(self, input_spectrum, mode, arm_names=None, interpolation_scheme="quadratic",
                    interpolation_pixels=3):
        """
        Estimate the RV of a spectrum on the basis of all of the 4MOST arms of either HRS or LRS.

        :param input_spectrum:
            A Spectrum object, containing an observed spectrum
        :param mode:
            The name of the 4MOST mode this arm is part of -- either LRS or HRS
        :param arm_names:
            A list of the 4MOST arms to use, or None to use all possible arms.
        :param interpolation_scheme:
            The type of function to use to interpolate the CCF to measure sub-pixel RVs.
        :param interpolation_pixels:
            The number of pixels around the peak of the CCF to use when interpolating to measure sub-pixel RVs.
        :return:
            [RV, error in RV]
        """

        assert interpolation_scheme in self.supported_interpolation_schemes()

        rv_estimates = []

        if arm_names is None:
            arm_names = self.templates_by_arm[mode].keys()

        # Compile a list of all the RV estimates, from all the arms and all the templates
        for arm_name in arm_names:
            new_rv_estimates = self.estimate_rv_from_single_arm(
                input_spectrum=self.resample_single_arm(input_spectrum=input_spectrum, arm_name=arm_name),
                mode=mode,
                arm_name=arm_name,
                interpolation_scheme=interpolation_scheme,
                interpolation_pixels=interpolation_pixels
            )
            rv_estimates.extend(new_rv_estimates)

        # Sort all the RV estimates into order of weight, and extract the stellar parameters of the best fitting
        # template
        rv_estimates.sort(key=itemgetter(1))
        rv_estimates_by_weight = rv_estimates.copy()
        stellar_parameters = rv_estimates[-1][2]

        # Sort all the RV estimates into order of RV, and chuck out the bottom and top quartiles. This
        # excludes estimates close to the speed of light from the subsequent statistics.
        rv_estimates.sort(key=itemgetter(0))
        if len(rv_estimates) >= 4:
            rv_estimates = rv_estimates[len(rv_estimates) // 4: len(rv_estimates) * 3 // 4]

        # Now form a weighted mean of all the RV estimates
        rv_mean = (sum([rv * weight for rv, weight, metadata in rv_estimates]) /
                   sum([weight for rv, weight, metadata in rv_estimates]))

        # Work out the standard deviation of the mean
        rv_variance = (sum([pow(rv - rv_mean, 2) * weight for rv, weight, metadata in rv_estimates]) /
                       sum([weight for rv, weight, metadata in rv_estimates]))

        rv_std_dev = sqrt(rv_variance)

        return rv_mean, rv_std_dev, stellar_parameters, rv_estimates_by_weight
