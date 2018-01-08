#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import norm
from multiprocessing import cpu_count
import logging
import itertools
import emcee
from emcee.interruptible_pool import InterruptiblePool
import random

import fourgp_speclib
import fourgp_degrade

logger = logging.getLogger(__name__)


class RvInstance(object):
    """
    A class which is adapted from Brani's code for forward-modelling the radial velocities of stars.
    """

    # Define the grid of template spectra supplied by Brani for use by RV code
    # [Axis variable, min, max, step]
    grid_axes = [["Teff", 4000, 8250, 250],
                 ["Fe/H", 0.5, 3.0, 0.5],
                 ["log_g", 1.5, 5.5, 0.5]
                 ]

    # Define initial guesses for each stellar label, together with minimum / maximum allowed values
    grid_axes_initial_guesses = {"sigma_gauss": 1.15, "Teff": 6000, "Fe/H": 1.5, "log_g": 4.5, "velocity": 0,
                                 "c0": 1, "c1": 0, "c2": 0}

    grid_axes_step_sizes = {"sigma_gauss": 0.02, "Teff": 250, "Fe/H": 0.05, "log_g": 0.1, "velocity": 100,
                            "c0": 100, "c1": 1, "c2": 0.00001}

    grid_axes_min = {"velocity": -500, "sigma_gauss": 1.001, "c0": -np.inf, "c1": -np.inf, "c2": -np.inf}
    grid_axes_max = {"velocity": 500, "sigma_gauss": 1.29, "c0": np.inf, "c1": np.inf, "c2": np.inf}

    for axis in grid_axes:
        grid_axes_min[axis[0]] = axis[1] - axis[3] / 2
        grid_axes_max[axis[0]] = axis[2] + axis[3] / 2

    # MCMC code requires a vector of parameters to explore. This is the ordering of parameters in the vector
    mcmc_parameter_order = ["velocity", "Teff", "Fe/H", "log_g", "sigma_gauss", "c0", "c1", "c2"]

    # Define the mesh of parameter values sampled in grid of template spectra
    grid_axis_values = [np.arange(axis[1], axis[2], axis[3]) for axis in grid_axes]

    # Total number of template spectra
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
        grid_axis_indices = [range(int((axis[2] - axis[1]) / axis[3])) for axis in self.grid_axes]
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

        # Default parameters for MCMC
        self.n_dim = 8  # Number of parameters in the model
        self.n_walkers = 160  # Number of MCMC walkers
        self.n_burn = 1000  # Length of the "burn-in" period to let chains stabilize
        self.n_steps = 1300

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

    # This method has to be static because class instances cannot be passed between threads if using a multiprocessing
    # pool of workers
    @staticmethod
    def pick_template_spectrum(template_library, grid_axes, axis_values):
        """
        Pick a template spectrum from a SpectrumArray which most closely matches the stellar parameters in the
        dictionary axis_values.

        :param template_library:
            A SpectrumArray object containing a library of template spectra. The grid of stellar parameters sampled
            should be specified by grid_axes.

        :param grid_axes:
            A list of the stellar parameters which are varied in template_library. Each axis of the grid is specified
            as [name, min, max, step size]. Stellar parameters are varied linearly. The final axis on the list varies
            the mostly rapidly with index in the SpectrumArray.

        :param axis_values:
            A dictionary of stellar parameters to which we should match a template spectrum.

        :return:
            A Spectrum object.
        """
        template_number = 0
        for axis in grid_axes:
            # How many samples do we have along this axis?
            axis_length = int(round((axis[2] - axis[1]) / axis[3]))
            template_number *= axis_length
            # Find the spectrum which most closely matches the requested value of this stellar parameter.
            axis_position = int(round(axis_values[axis[0]] - axis[1]) / axis[3])
            # Clip to the limits of ths axis
            axis_position = min(max(axis_position, 0), axis_length - 1)
            # Multiply-and-accumulate the index of the spectrum we want
            template_number += axis_position
        return template_library.extract_item(template_number)

    # This method has to be static because class instances cannot be passed between threads if using a multiprocessing
    # pool of workers
    @staticmethod
    def log_probability(theta, template_library, observed_spectrum, grid_axes):
        """
        This is the log-probability function which we use an MCMC chain to sample.

        :param theta:
            A vector containing floating point values for: [velocity, t_eff, fe_h, log_g, sigma_gauss, c0, c1, c2]

        :param template_library:
            A SpectrumArray object containing a grid of template continuum-normalised spectra at various stellar
            parameter values.

        :param observed_spectrum:
            A Spectrum object containing the observing spectrum whose radial velocity we're trying to estimate.

        :param grid_axes:
            A list of the stellar-parameter axes sampled by template_library.

        :return:
            A floating-point log-probability value.
        """

        # Unpack stellar parameters from vector passed by optimiser
        # This must match self.mcmc_parameter_order above
        # velocity has units of km/s
        velocity, t_eff, fe_h, log_g, sigma_gauss, c0, c1, c2 = theta

        # Return a probability of minus infinity if we are outside bounds of valid parameter space
        if (np.abs(velocity) > 500) or (sigma_gauss < 1) or (sigma_gauss > 1.3):
            return -np.inf

        # Check that stellar parameters are within the range of values spanned by template spectra
        for axis_no, axis_value in enumerate([t_eff, fe_h, log_g]):
            axis = grid_axes[axis_no]
            if ((axis_value <= axis[1] - axis[3] / 2) or
                    (axis_value >= axis[2] + axis[3] / 2)):
                return -np.inf

        # Pick a template
        template_spectrum = RvInstance.pick_template_spectrum(template_library=template_library,
                                                              grid_axes=grid_axes,
                                                              axis_values={"Teff": t_eff, "Fe/H": fe_h, "log_g": log_g})

        if any(np.isnan(template_spectrum.values)):
            return -np.inf

        # Convolve the template spectrum with a Gaussian
        convolver = fourgp_degrade.SpectrumConvolver(template_spectrum)
        template_convolved = convolver.gaussian_convolve(sigma_gauss)

        # Shift in wavelength
        template_observer_frame = template_convolved.apply_radial_velocity(velocity * 1000)  # Unit m/s

        # Interpolate the convolved, Doppler-shifted template onto the observed spectrum's wavelength
        interpolator = fourgp_degrade.SpectrumInterpolator(template_observer_frame)
        template_resampled = interpolator.match_to_other_spectrum(other=observed_spectrum,
                                                                  interpolate_errors=False, interpolate_mask=False)

        # Multiply the template spectrum by the observed spectrum's continuum
        template_continuum = fourgp_speclib.SpectrumPolynomial(wavelengths=template_resampled.wavelengths,
                                                               terms=3,
                                                               coefficients=(c0, c1, c2))
        template_with_continuum = template_continuum * template_resampled

        # Mask out bad data
        mask = (observed_spectrum.mask * np.isfinite(observed_spectrum.values) * (observed_spectrum.value_errors > 0) *
                np.isfinite(template_with_continuum.values))

        # log likelihood
        log_likelihood = np.sum(norm.logpdf(x=observed_spectrum.values[mask],
                                            loc=template_with_continuum.values[mask],
                                            scale=observed_spectrum.value_errors[mask]))
        log_likelihood += norm.logpdf(x=sigma_gauss, loc=1.15, scale=0.02)
        log_likelihood += norm.logpdf(x=velocity, loc=0, scale=150)

        return log_likelihood - 100000

    def fit_rv(self, observed_spectrum, initial_guesses=None):
        """
        Estimate the radial velocity of the observed specrtrum <observed_spectrum>.

        :param observed_spectrum:
            A Spectrum object containing the observed spectrum we are to fit.

        :param initial_guesses:
            An optional dictionary of initial guesses for the stellar parameters of this spectrum. If it is not
            supplied, default values are assumed.

        :return:
            A dictionary of stellar parameters, including a radial velocity with key "velocity" (units km/s).
        """

        # Initial guesses for stellar parameters
        stellar_labels = self.grid_axes_initial_guesses.copy()

        if initial_guesses is not None:
            stellar_labels.update(initial_guesses)

        # Look up the range of wavelengths spanned by the template spectra
        template_min_wavelength = self._template_spectra.wavelengths[0]
        template_max_wavelength = self._template_spectra.wavelengths[-1]

        # Truncate observed spectra so that even when the maximum allowed radial velocity is applied, it won't reach
        # outside the wavelength range of the template spectra
        v_max = self.grid_axes_max["velocity"] * 1.2  # Add a bit of a safety margin
        c = 299792458.0
        observed_spectrum.mask_exclude(wavelength_max=template_min_wavelength * (1 + v_max / c))
        observed_spectrum.mask_exclude(wavelength_min=template_max_wavelength / (1 + v_max / c))
        observed_truncated = observed_spectrum.truncate_to_mask()

        # Make spectrum shared so that it can be passed among threads in a worker pool
        observed_shared = fourgp_speclib.SpectrumArray.from_spectra(spectra=[observed_truncated], shared_memory=True). \
            extract_item(0)

        # Make initial continuum fit to observed spectrum, using small wavelength window at 5650 to 5820 A
        observed_continuum_fitter = fourgp_speclib.SpectrumSmoothFactory(
            wavelengths=observed_shared.wavelengths,
            terms=3,
            function_family=fourgp_speclib.SpectrumPolynomial)

        # Pick a template continuum normalised spectrum
        template_continuum_normalised = RvInstance.pick_template_spectrum(template_library=self._template_spectra,
                                                                          grid_axes=self.grid_axes,
                                                                          axis_values=stellar_labels)
        interpolator = fourgp_degrade.SpectrumInterpolator(template_continuum_normalised)
        template_interpolated = interpolator.match_to_other_spectrum(observed_shared)

        # Fit continuum to observed spectrum, assuming its absorption features matched by template at zero rv
        observed_continuum_fit = observed_continuum_fitter.fit_to_continuum_via_template(
            other=observed_shared,
            template=template_interpolated,
            lambda_min_norm=5650, lambda_max_norm=5820)

        # Add coefficients of initial continuum fit to dictionary of stellar parameters
        stellar_labels["c0"] = observed_continuum_fit.coefficients[0]
        stellar_labels["c1"] = observed_continuum_fit.coefficients[1]
        stellar_labels["c2"] = observed_continuum_fit.coefficients[2]

        # Initialise starting points for MCMC walkers
        theta0 = [stellar_labels[x] for x in self.mcmc_parameter_order]
        walker_positions = np.random.normal(loc=theta0,
                                            scale=[self.grid_axes_step_sizes[x] for x in self.mcmc_parameter_order],
                                            size=(self.n_walkers, self.n_dim))

        # Clip the positions of the walkers to appropriate ranges
        walker_positions = np.clip(a=walker_positions,
                                   a_min=[self.grid_axes_min[x] for x in self.mcmc_parameter_order],
                                   a_max=[self.grid_axes_max[x] for x in self.mcmc_parameter_order])
        # print [RvInstance.log_probability(walker_positions[i], self._template_spectra, observed_shared, self.grid_axes) for i in range(self.n_walkers)]

        # Start workers
        pool = InterruptiblePool(processes=self._threads)

        # initialize the sampler
        sampler = emcee.EnsembleSampler(nwalkers=self.n_walkers, dim=self.n_dim, lnpostfn=RvInstance.log_probability,
                                        # pool=pool,
                                        kwargs={"template_library": self._template_spectra,
                                                "observed_spectrum": observed_shared,
                                                "grid_axes": self.grid_axes})

        # burn-in the chains
        sampler.run_mcmc(pos0=walker_positions, N=self.n_burn)

        med = np.median(a=sampler.lnprobability[:, -1])
        rms = 0.741 * (np.percentile(a=sampler.lnprobability[:, -1], q=75) -
                       np.percentile(a=sampler.lnprobability[:, -1], q=25))

        # Determine the starting point after burn-in. Eliminate bad chains
        good_chains = sampler.lnprobability[:, -1] > (med - 3 * rms)

        median_params = np.median(a=sampler.chain[good_chains, -1, :], axis=0)
        rms_params = 0.741 * (np.percentile(a=sampler.chain[good_chains, -1, :], q=75, axis=0) -
                              np.percentile(a=sampler.chain[good_chains, -1, :], q=25, axis=0))
        best = np.random.normal(loc=median_params,
                                scale=rms_params,
                                size=(self.n_walkers, self.n_dim))

        # clip the guesses to appropriate ranges
        best = np.clip(a=best,
                       a_min=[self.grid_axes_min[x] for x in self.mcmc_parameter_order],
                       a_max=[self.grid_axes_max[x] for x in self.mcmc_parameter_order])

        # Reset the chains to remove the burn-in samples.
        sampler.reset()

        # Run the chains for real
        sampler.run_mcmc(pos0=best, N=self.n_steps)
        pool.terminate()

        max_prob = sampler.flatchain[np.argmax(sampler.flatlnprobability)]
        output = {}
        for i, par in enumerate(self.mcmc_parameter_order):
            output[par] = max_prob[i]

        return output


def random_radial_velocity():
    """
    Pick a random radial velocity in km/s, following the approximate probability distribution expected for 4MOST
    galactic stars.

    :return:
        Radial velocity in km/s
    """

    distribution_selector = random.uniform(a=0, b=100)

    if distribution_selector < 10:
        # Pick 10% of stars from a uniform distribution from -200 to 200
        return random.uniform(a=-200, b=200)  # Unit km/s
    else:
        return random.gauss(mu=0, sigma=25)
