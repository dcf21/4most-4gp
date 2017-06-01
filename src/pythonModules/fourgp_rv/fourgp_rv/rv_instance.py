#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from multiprocessing import cpu_count
import logging
import itertools

import fourgp_speclib
import fourgp_degrade

logger = logging.getLogger(__name__)


class RvInstance(object):
    """
    A class which is adapted from Brani's code for forward-modelling the radial velocities of stars.
    """

    # Define the grid of template spectra supplied by Brani for use by RV code
    # [Axis variable, (min,max,step)]
    grid_axes = [["Teff", (4000, 8250, 250)],
                 ["Fe/H", (0.5, 3.0, 0.5)],
                 ["logg", (1.5, 5.5, 0.5)]
                 ]

    grid_axis_values = [np.arange(axis[1][0], axis[1][1], axis[1][2])
                        for axis in grid_axes]

    expected_number_spectra = np.prod([i.shape[0] for i in grid_axis_values])

    def __init__(self, spectrum_library, threads=None):
        """
        Instantiate the RV code, and read from disk the library of template spectra.
        
        :param spectrum_library:
            A SpectrumLibrary containing the template spectra we use for modelling.

        :type spectrum_library:
            SpectrumLibrary

        :param threads:
            The number of concurrent CPU threads to use.

        :type threads:
            int
        """

        assert isinstance(spectrum_library, fourgp_speclib.SpectrumLibrary), \
            "Argument to RvInstance should be a SpectrumLibrary."
        self._spectrum_library = spectrum_library

        # Work out how many CPUs we should allow the Cannon to use
        if threads is None:
            threads = cpu_count()
        self._threads = threads

        # Load template spectra
        spectrum_list = self._spectrum_library.search()

        # Check that we have the expected number of template spectra
        assert len(spectrum_list) == self.expected_number_spectra, \
            "Supplied SpectrumLibrary did not contain the expected number of spectra. " \
            "Was expecting {} spectra; got {} spectra.". \
                format(self.expected_number_spectra, len(spectrum_list))

        # Check that the template spectra are sampled on the right grid
        grid_axis_indices = [range(int((axis[1][1] - axis[1][0]) / axis[1][2])) for axis in self.grid_axes]
        grid_axis_index_combinations = itertools.product(*grid_axis_indices)

        # Check that each grid point in turn exists
        grid_spectrum_ids = []
        for axis_indices in grid_axis_index_combinations:
            search_criteria = {}
            for axis_counter, index in enumerate(axis_indices):
                metadata_key = self.grid_axes[axis_counter][0] + "_index"
                search_criteria[metadata_key] = index
            matches = self._spectrum_library.search(**search_criteria)
            assert len(matches) == 1, "Could not find spectrum matching {}".format(search_criteria)
            grid_spectrum_ids.append(matches[0]['specId'])

        # Load library of template spectra
        self._template_spectra = self._spectrum_library.open(ids=grid_spectrum_ids, shared_memory=True)

    @classmethod
    def from_spectrum_library_sqlite(cls, library_path, *args, **kwargs):
        """
        Open a Spectrum Library on disk at supplied file path, and instantiate an RvInstance using the template
        spectra in it.

        :param library_path:
            The file path of the library of template spectra to use.

        :type library_path:
            str

        :return:
            RvInstance
        """
        library = fourgp_speclib.SpectrumLibrarySqlite(path=library_path)
        return cls(spectrum_library=library, *args, **kwargs)

    @staticmethod
    def log_probability(theta, def_param):

        # Unpack stellar parameters from vector passed by optimiser
        velocity, t_eff, fe_h, log_g, sigma_gauss, c0, c1, c2 = theta

        # Unpack additional parameters
        template_library, observed_spectrum = def_param

        # Return a probability of minus infinity if we are outside bounds of valid parameter space
        if ((np.abs(velocity) > 500) or
                (t_eff < 3875) or (t_eff > 8125) or
                (fe_h < 0.25) or (fe_h > 2.75) or
                (log_g < 1.25) or (log_g > 5.25) or
                (sigma_gauss < 1) or (sigma_gauss > 1.3)
            ):
            return -np.inf

        # Pick a template
        template_spectrum = template_library[int(np.abs(np.rint((t_eff - 4000) / 250.))),
                            int(np.abs(np.rint((fe_h - 0.5) / 0.5))),
                            int(np.abs(np.rint((log_g - 1.5) / 0.5))), :]

        if any(np.isnan(template_spectrum)):
            return -np.inf

        # Convolve the template spectrum with a Gaussian
        convolver = fourgp_degrade.SpectrumConvolver(template_spectrum)
        template_convolved = convolver.gaussian_convolve(sigma_gauss)

        # Shift in wavelength
        template_observer_frame = template_convolved.apply_radial_velocity(velocity * 1000)

        # Interpolate the convolved, Doppler-shifted template onto the observed spectrum's wavelength
        interpolator = fourgp_degrade.SpectrumInterpolator(template_observer_frame)
        template_resampled = interpolator.match_to_other_spectrum(other=observed_spectrum,
                                                                  interpolate_errors=False, interpolate_mask=False)

        # Multiply the template spectrum by the observed spectrum's continuum
        flux_interp = model_spec([c0, c1, c2], wave, flux_interp)

        # log likelihood
        lnL = np.sum(norm.logpdf(flux, loc=flux_interp, scale=error))
        return lnL + norm.logpdf(sigma_gauss, loc=1.15, scale=0.02) + norm.logpdf(velocity, loc=0, scale=150) - 100000

    def fit_rv(self, spectrum):
        shared_spectrum = fourgp_speclib.SpectrumArray.from_spectrum(spectrum=spectrum, shared_memory=True)
