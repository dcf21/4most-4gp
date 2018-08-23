#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
This class adds Gaussian noise to spectra. The noise level is specified as an SNR/pixel, which is specified within
a wavelength window. We assume a Gaussian instrumental profile with which equal to the mean pixel spacing. We
automatically detect if multiple arms with different spectral resolution have been spliced together into a single raster
and we treat each arm separately.
"""

import numpy as np
import logging

from fourgp_speclib import Spectrum
from .convolve import SpectrumConvolver
from .resample import SpectrumResampler
from .spectrum_properties import SpectrumProperties

logger = logging.getLogger(__name__)


class GaussianNoise:
    def __init__(self,
                 wavelength_raster,
                 snr_list=None,
                 snr_definitions=None,
                 use_snr_definitions=None
                 ):
        """
        Instantiate a class for adding Gaussian noise to spectra.

        :param wavelength_raster:
            List of the wavelength raster for the output spectrum.

        :param snr_list:
            List of the SNRs that we want 4FS to degrade input spectra to.

        :param snr_definitions:
            List of ways we define SNR. Each should take the form of a tuple (name,min,max), where we take the median
            SNR per pixel between the specified minimum and maximum wavelengths in Angstrom.

        :param use_snr_definitions:
            List of the SNR definitions to use for each wavelength arm.
        """

        # Divide wavelength raster into spectral arms, which have distinct pixel spacing
        # We convolve each wavelength arm separately
        self.wavelength_arms = SpectrumProperties(wavelength_raster).wavelength_arms()['wavelength_arms']
        self.wavelength_raster = wavelength_raster

        logger.info("Detected {} wavelength arms".format(len(self.wavelength_arms)))

        # Read the list of SNRs that we have been requested to degrade spectra to
        if snr_list is None:
            # Default list
            snr_list = (10, 12, 14, 16, 18, 20, 23, 26, 30, 35, 40, 45, 50, 80, 100, 130, 180, 250)
        self.snr_list = snr_list

        # Read the list of SNR definitions supplied, or default to using a window between 6180 and 6680 A
        if snr_definitions is None:
            snr_definitions = [("MEDIANSNR", 6180, 6680)]

        if use_snr_definitions is None:
            use_snr_definitions = ("MEDIANSNR",) * len(self.wavelength_arms)

        self.snr_definitions = snr_definitions
        self.use_snr_definitions = use_snr_definitions

        assert len(self.use_snr_definitions) == len(self.wavelength_arms), \
            "Need an SNR definition for each wavelength arm. " \
            "Received {} definitions, but autodetected {} arms.". \
                format(len(self.use_snr_definitions), len(self.wavelength_arms))

    def close(self):
        # Do cleanup...
        pass

    def process_spectra(self, spectra_list):
        """
        Add Gaussian noise to a list of 4GP Spectrum objects.

        :param spectra_list:
            A list of the spectra we should pass through 4FS. Each entry in the list should be a list of tuple with
            two entries: (input_spectrum, input_spectrum_continuum_normalised). These reflect the contents of the
            third and second columns of Turbospectrum's ASCII output respectively.

        :type spectra_list:
            (list, tuple) of (list, tuple) of Spectrum objects
        """

        output = []  # output[ spectrum_number ][ snr ] = [ full_spectrum, continuum normalised ]

        for spectrum in spectra_list:

            # Convolve and resample onto new wavelength raster.
            # Each wavelength arm is separately convolved by mean pixel spacing.
            resampled_spectrum = []  # resampled_spectrum[ 0=full spectrum ; 1=continuum normalised ][ wavelength_arm ]
            for index, item in enumerate(spectrum):
                resampled_spectrum.append([])
                for (raster, pixel_spacing) in self.wavelength_arms:
                    convolver = SpectrumConvolver(item)
                    convolved = convolver.gaussian_convolve(pixel_spacing)
                    resampler = SpectrumResampler(convolved)
                    resampled = resampler.onto_raster(raster)
                    resampled_spectrum[-1].append(resampled.values)

            # Calculate continuum spectrum by dividing the flux normalised spectrum by continuum normalised spectrum
            continuum_per_arm = []
            for index_arm in range(len(self.wavelength_arms)):
                continuum_per_arm.append(resampled_spectrum[0][index_arm] / resampled_spectrum[1][index_arm])
            continuum = np.concatenate(continuum_per_arm)

            # Measure the integrated signal within the range of each SNR definition
            mean_signal_per_pixel = {}
            for snr_definition_name, wavelength_min, wavelength_max in self.snr_definitions:
                indices = (wavelength_min <= self.wavelength_raster) * (self.wavelength_raster <= wavelength_max)
                pixel_count = np.sum(indices)
                pixel_sum = np.sum(continuum[indices])
                mean_signal_per_pixel[snr_definition_name] = pixel_sum / pixel_count

            # Add noise to spectra at each SNR in turn
            output_item = {}
            output.append(output_item)
            for snr_value in self.snr_list:
                output_values = np.zeros(0)
                output_value_errors = np.zeros(0)
                output_values_cn = np.zeros(0)
                output_value_errors_cn = np.zeros(0)

                # Add noise to each wavelength arm individually
                for index_arm, (raster, pixel_spacing) in enumerate(self.wavelength_arms):
                    snr_definition = self.use_snr_definitions[index_arm]
                    signal_level = mean_signal_per_pixel[snr_definition]
                    noise_level = signal_level / snr_value
                    arm_length = len(raster)

                    # Synthesize some Gaussian noise
                    noise = np.random.normal(loc=0,
                                             scale=noise_level,
                                             size=arm_length)

                    # Add noise into flux-normalised spectrum
                    noised_signal = resampled_spectrum[0][index_arm] + noise
                    noised_signal_errors = np.ones_like(noised_signal) * noise_level

                    # Compute the new continuum-normalised spectrum by dividing by the pure continuum we computed
                    noised_signal_cn = noised_signal / continuum_per_arm[index_arm]
                    noised_signal_errors_cn = noised_signal_errors / continuum_per_arm[index_arm]

                    # Concatenate the various wavelength arms together
                    output_values = np.append(output_values, noised_signal)
                    output_value_errors = np.append(output_value_errors, noised_signal_errors)
                    output_values_cn = np.append(output_values_cn, noised_signal_cn)
                    output_value_errors_cn = np.append(output_value_errors_cn, noised_signal_errors_cn)

                # Convert output spectra into Spectrum objects
                metadata = spectrum[0].metadata.copy()
                metadata['continuum_normalised'] = 0
                metadata['SNR'] = float(snr_value)
                output_spectrum = Spectrum(wavelengths=self.wavelength_raster,
                                           values=output_values,
                                           value_errors=output_value_errors,
                                           metadata=metadata.copy())

                metadata['continuum_normalised'] = 1
                output_spectrum_cn = Spectrum(wavelengths=self.wavelength_raster,
                                              values=output_values_cn,
                                              value_errors=output_value_errors_cn,
                                              metadata=metadata.copy())

                # Add to output data structure
                output_item[snr_value] = (output_spectrum, output_spectrum_cn)

        # Return output spectra to user
        # output[ spectrum_number ][ snr ] = [ flux normalised , continuum normalised ]
        return output
