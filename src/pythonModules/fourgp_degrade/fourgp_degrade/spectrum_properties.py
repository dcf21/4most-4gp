# -*- coding: utf-8 -*-

"""
A class which converts a wavelength raster into separate arms with different pixel spacing, which we should
degrade separately because they are observed separately.
"""

import numpy as np


class SpectrumProperties:
    """
    A class which converts a wavelength raster into separate arms with different pixel spacing, which we should
    degrade separately because they are observed separately.
    """

    def __init__(self, wavelength_raster):
        self.wavelength_raster = wavelength_raster

    def wavelength_arms(self):
        wavelength_raster = self.wavelength_raster

        # Process wavelength raster into spectral arms
        wavelength_arms = []
        break_points = []

        wavelength_old = wavelength_raster[0]
        wavelength_delta = None
        arm_raster = [wavelength_raster[0]]
        arm_pixel_gaps = []
        for i in range(1, len(wavelength_raster)):
            wavelength_delta_new = wavelength_raster[i] - wavelength_old
            wavelength_old = wavelength_raster[i]
            if wavelength_delta is not None:
                ratio = wavelength_delta_new / wavelength_delta
                if ratio < 1:
                    ratio = 1 / ratio
                if ratio > 1.02:
                    mean_interval = np.mean(np.asarray(arm_pixel_gaps))
                    wavelength_arms.append([np.asarray(arm_raster), mean_interval])
                    wavelength_delta = None
                    break_points.append((wavelength_raster[i] + wavelength_raster[i - 1]) / 2)
                    arm_raster = [wavelength_raster[i]]
                    arm_pixel_gaps = []
                    continue
            wavelength_delta = wavelength_delta_new
            arm_raster.append(wavelength_raster[i])
            arm_pixel_gaps.append(wavelength_delta_new)

        # Add final wavelength arm to self.wavelength_arms
        mean_interval = np.mean(np.asarray(arm_pixel_gaps))
        wavelength_arms.append([np.asarray(arm_raster), mean_interval])

        return {
            "wavelength_arms": wavelength_arms,
            "break_points": break_points
        }
