# -*- coding: utf-8 -*-

"""
This provides a wrapper to the 2016 version of Andy Casey's Cannon.

The code for this Cannon is in the master branch of this repository: https://github.com/dcf21/AnniesLasso

This is the version of the Cannon used in all technical reports issued by Dominic Ford.
"""

import numpy as np
from multiprocessing import cpu_count
import logging
from astropy.table import Table
import AnniesLasso as tc

import fourgp_speclib
from fourgp_degrade import SpectrumResampler

logger = logging.getLogger(__name__)


class CannonInstanceCaseyOld(object):
    """
    A class which holds an instance of the Cannon, and provides convenience methods for training it on arrays of spectra
    loaded from 4GP SpectrumLibrary objects.
    """

    def __init__(self, training_set, label_names, wavelength_arms=None,
                 censors=None, progress_bar=False, threads=None, tolerance=None, polynomial_order=2,
                 load_from_file=None, debugging=False):
        """
        Instantiate the Cannon and train it on the spectra contained within a SpectrumArray.

        :param training_set:
            A SpectrumArray containing the spectra to train the Cannon on.

        :param label_names:
            A list of the names of the labels the Cannon is to estimate. We require that all of the training spectra
            have metadata fields defining all of these labels.

        :param wavelength_arms:
            A list of the wavelength break-points between arms which should have continuum fitted separately. For
            compatibility we accept this argument, but it is not used for continuum-normalised spectra.

        :param threads:
            The number of CPU cores we should use. If None, we look up how many cores this computer has.

        :param tolerance:
            The tolerance xtol which the method <scipy.optimize.fmin_powell> uses to determine convergence.

        :param polynomial_order:
            The order of polynomials to use as fitting functions within the Cannon.

        :param load_from_file:
            The filename of the internal state of a pre-trained Cannon, which we should load rather than doing
            training from scratch.

        :param debugging:
            Boolean flag determining whether we produce debugging output

        :type debugging:
            bool
        """

        self._debugging_output_counter = 0
        self._debugging = debugging
        self.cannon_version = tc.__version__
        self._wavelength_arms = wavelength_arms
        logger.info("Wavelength arm breakpoints: {}".format(self._wavelength_arms))

        assert isinstance(training_set, fourgp_speclib.SpectrumArray), \
            "Training set for the Cannon should be a SpectrumArray."

        # Hook for normalising input spectra
        training_set = self.normalise(training_set)

        self._training_set = training_set
        self._progress_bar = progress_bar

        # Work out how many CPUs we should allow the Cannon to use
        if threads is None:
            threads = cpu_count()

        # Turn error bars on fluxes into inverse variances
        inverse_variances = training_set.value_errors ** (-2)

        # Flag bad data points
        ignore = (training_set.values < 0) + ~np.isfinite(inverse_variances)
        inverse_variances[ignore] = 0
        training_set.values[ignore] = 1

        # Check that labels are correctly set in metadata
        for index in range(len(training_set)):
            metadata = training_set.get_metadata(index)
            for label in label_names:
                assert label in metadata, "Label <{}> not set on training spectrum number {}. " \
                                          "Labels on this spectrum are: {}.".format(
                    label, index, ", ".join(list(metadata.keys())))
                assert np.isfinite(metadata[label]), "Label <{}> is not finite on training spectrum number {}. " \
                                                     "Labels on this spectrum are: {}.".format(
                    label, index, metadata)

        # Compile table of training values of labels from metadata contained in SpectrumArray
        training_label_values = Table(names=label_names,
                                      rows=[[training_set.get_metadata(index)[label] for label in label_names]
                                            for index in range(len(training_set))])

        self._model = tc.L1RegularizedCannonModel(labelled_set=training_label_values,
                                                  normalized_flux=training_set.values,
                                                  normalized_ivar=inverse_variances,
                                                  dispersion=training_set.wavelengths,
                                                  threads=threads)

        self._model.vectorizer = tc.vectorizer.NormalizedPolynomialVectorizer(
            labelled_set=training_label_values,
            terms=tc.vectorizer.polynomial.terminator(label_names, polynomial_order)
        )

        if censors is not None:
            self._model.censors = censors

        self._model.s2 = 0
        self._model.regularization = 0

        if load_from_file is None:
            logger.info("Starting to train the Cannon")
            self._model.train(
                progressbar=self._progress_bar,
                op_kwargs={'xtol': tolerance, 'ftol': tolerance},
                op_bfgs_kwargs={}
            )
            logger.info("Cannon training completed")
        else:
            logger.info("Loading Cannon from disk")
            self._model.load(filename=load_from_file)
            logger.info("Cannon loaded successfully")
        self._model._set_s2_by_hogg_heuristic()

    def fit_spectrum(self, spectrum):
        """
        Fit stellar labels to a continuum-normalised spectrum.

        :param spectrum:
            A Spectrum object containing the spectrum for the Cannon to fit.

        :type spectrum:
            Spectrum

        :return:
        """

        assert isinstance(spectrum, fourgp_speclib.Spectrum), \
            "Supplied spectrum for the Cannon to fit is not a Spectrum object."

        assert spectrum.raster_hash == self._training_set.raster_hash, \
            "Supplied spectrum for the Cannon to fit is not sampled on the same raster as the training set."

        # Hook for normalising input spectra
        spectrum = self.normalise(spectrum)

        inverse_variances = spectrum.value_errors ** (-2)

        # Ignore bad pixels.
        bad = (spectrum.value_errors < 0) + (~np.isfinite(inverse_variances * spectrum.values))
        inverse_variances[bad] = 0
        spectrum.values[bad] = np.nan

        labels, cov, meta = self._model.fit(
            normalized_flux=spectrum.values,
            normalized_ivar=inverse_variances,
            progressbar=False,
            full_output=True)

        return labels, cov, meta

    def normalise(self, spectrum):
        """
        This is a hook for doing some kind of normalisation on spectra. Not implemented in this base class.

        :param spectrum:
            The spectrum to be normalised.
        :return:
            Normalised version of this spectrum.
        """

        return spectrum

    def save_model(self, filename, overwrite=True):
        """
        Save the parameters of a trained Cannon instance.

        :param filename:
            Output file in which to save trained Cannon parameters.

        :type filename:
            str

        :param overwrite:
            Do we overwrite any existing files?

        :type overwrite:
            bool

        :return:
            None
        """
        self._model.save(filename=filename,
                         include_training_data=False,
                         overwrite=overwrite)

    def __str__(self):
        return "<{module}.{name} instance".format(module=self.__module__,
                                                  name=type(self).__name__)

    def __repr__(self):
        return "<{0}.{1} object at {2}>".format(self.__module__,
                                                type(self).__name__, hex(id(self)))


class CannonInstanceCaseyOldWithRunningMeanNormalisation(CannonInstanceCaseyOld):
    """
    A class which holds an instance of the Cannon, and automatically normalises the test spectra using a running mean
    with a width of the certain number of pixels. This allows the Cannon to operate on spectra which have not been
    continuum normalised.

    The training and test spectra are both normalised in the same way.
    """

    def __init__(self, training_set, label_names, wavelength_arms, normalisation_window=300,
                 censors=None, progress_bar=False, threads=None, tolerance=None, polynomial_order=2,
                 load_from_file=None, debugging=False):
        """
        Instantiate the Cannon and train it on the spectra contained within a SpectrumArray.

        :param training_set:
            A SpectrumArray containing the spectra to train the Cannon on.

        :param label_names:
            A list of the names of the labels the Cannon is to estimate. We require that all of the training spectra
            have metadata fields defining all of these labels.

        :param wavelength_arms:
            A list of the wavelength break-points between arms which should have continuum fitted separately.

        :type wavelength_arms:
            List of wavelengths in A.

        :param normalisation_window:
            The width of the running mean normalisation window.

        :type normalisation_window:
            int

        :param threads:
            The number of CPU cores we should use. If None, we look up how many cores this computer has.

        :param tolerance:
            The tolerance xtol which the method <scipy.optimize.fmin_powell> uses to determine convergence.

        :param polynomial_order:
            The order of polynomials to use as fitting functions within the Cannon.

        :param load_from_file:
            The filename of the internal state of a pre-trained Cannon, which we should load rather than doing
            training from scratch.

        :param debugging:
            Boolean flag determining whether we produce debugging output

        :type debugging:
            bool
        """

        self._window_width = normalisation_window

        # Initialise
        super(CannonInstanceCaseyOldWithRunningMeanNormalisation, self).__init__(training_set=training_set,
                                                                                 label_names=label_names,
                                                                                 wavelength_arms=wavelength_arms,
                                                                                 censors=censors,
                                                                                 progress_bar=progress_bar,
                                                                                 threads=threads,
                                                                                 tolerance=tolerance,
                                                                                 polynomial_order=polynomial_order,
                                                                                 load_from_file=load_from_file,
                                                                                 debugging=debugging
                                                                                 )

    def normalise(self, spectrum):
        """
        This is a hook for doing some kind of normalisation on spectra.

        :param spectrum:
            The spectrum to be normalised.
        :return:
            Normalised version of this spectrum.
        """

        # If we're passed a spectrum array, normalise each spectrum in turn
        if isinstance(spectrum, fourgp_speclib.SpectrumArray):
            l = len(spectrum)
            for i in range(l):
                spectrum_item = spectrum.extract_item(i)
                spectrum_normalised = self.normalise(spectrum_item)
                spectrum_item.values[:] = spectrum_normalised.values
                spectrum_item.value_errors[:] = spectrum_normalised.value_errors
            return spectrum

        assert isinstance(spectrum, fourgp_speclib.Spectrum), \
            "The CannonInstance.normalise method requires a Spectrum object as input."

        if self._debugging:
            self._debugging_output_counter += 1

        # Returns an array of length len(x)-(N-1)
        def running_mean(x, n):
            cumulative_sum = np.cumsum(np.insert(x, 0, 0))
            return (cumulative_sum[n:] - cumulative_sum[:-n]) / float(n)

        # Work out the raster of pixels inside each wavelength arm
        raster = spectrum.wavelengths
        lower_cut = 0
        arm_rasters = []
        for break_point in self._wavelength_arms:
            arm_rasters.append((raster >= lower_cut) * (raster < break_point))
            lower_cut = break_point
        arm_rasters.append(raster >= lower_cut)

        output_wavelengths = []
        output_values = []
        output_value_errors = []

        for arm in arm_rasters:
            output_wavelengths.append(raster[arm])
            input_values = spectrum.values[arm]
            input_errors = spectrum.value_errors[arm]

            normalisation = running_mean(input_values, self._window_width)
            padding_needed = len(input_values) - len(normalisation)
            padding_left = int(padding_needed / 2)
            padding_right = padding_needed - padding_left
            normalisation_full = np.concatenate([np.repeat(normalisation[0], padding_left),
                                                 normalisation,
                                                 np.repeat(normalisation[-1], padding_right)
                                                 ])

            output_values.append(input_values / normalisation_full)
            output_value_errors.append(input_errors / normalisation_full)

        output = fourgp_speclib.Spectrum(wavelengths=np.concatenate(output_wavelengths),
                                         values=np.concatenate(output_values),
                                         value_errors=np.concatenate(output_value_errors),
                                         metadata=spectrum.metadata)

        # Produce debugging output if requested
        if self._debugging:
            np.savetxt("/tmp/debug_{:06d}.txt".format(self._debugging_output_counter),
                       np.transpose([raster, spectrum.values, spectrum.value_errors]))

        return output


class CannonInstanceCaseyOldWithContinuumNormalisation(CannonInstanceCaseyOld):
    """
    A class which holds an instance of the Cannon, and automatically continuum normalises the test spectra in
    an iterative fashion, using pixels which are known to have continuum normalised values close to one.

    The training spectra must be continuum normalised.

    This class is currently very bad, and should not be expected to produce meaningful results. For the time being,
    you should continuum normalise spectra yourself and use the CannonInstance class instead.
    """

    def __init__(self, training_set, label_names, wavelength_arms,
                 continuum_model_family=fourgp_speclib.SpectrumPolynomial,
                 censors=None, progress_bar=False, threads=None, tolerance=None, polynomial_order=2,
                 load_from_file=None, debugging=False):
        """
        Instantiate the Cannon and train it on the spectra contained within a SpectrumArray.

        :param training_set:
            A SpectrumArray containing the spectra to train the Cannon on.

        :param label_names:
            A list of the names of the labels the Cannon is to estimate. We require that all of the training spectra
            have metadata fields defining all of these labels.

        :param wavelength_arms:
            A list of the wavelength break-points between arms which should have continuum fitted separately.

        :type wavelength_arms:
            List of wavelengths in A.

        :param continuum_model_family:
            A class defining smooth functions of some form, which we use to fit the continuum in this spectrum.

        :type continuum_model_family:
            SpectrumSmooth subclass

        :param threads:
            The number of CPU cores we should use. If None, we look up how many cores this computer has.

        :param tolerance:
            The tolerance xtol which the method <scipy.optimize.fmin_powell> uses to determine convergence.

        :param polynomial_order:
            The order of polynomials to use as fitting functions within the Cannon.

        :param load_from_file:
            The filename of the internal state of a pre-trained Cannon, which we should load rather than doing
            training from scratch.

        :param debugging:
            Boolean flag determining whether we produce debugging output

        :type debugging:
            bool
        """

        assert issubclass(continuum_model_family, fourgp_speclib.SpectrumSmooth), \
            "Input continuum model family must be a subclass of <SpectrumSmooth>."

        self._continuum_model_family = continuum_model_family

        # Initialise
        super(CannonInstanceCaseyOldWithContinuumNormalisation, self).__init__(training_set=training_set,
                                                                               label_names=label_names,
                                                                               wavelength_arms=wavelength_arms,
                                                                               censors=censors,
                                                                               progress_bar=progress_bar,
                                                                               threads=threads,
                                                                               tolerance=tolerance,
                                                                               polynomial_order=polynomial_order,
                                                                               load_from_file=load_from_file,
                                                                               debugging=debugging
                                                                               )

    def fit_spectrum(self, spectrum):
        """
        Fit stellar labels to a spectrum which has not been continuum normalised.

        :param spectrum:
            A Spectrum object containing the spectrum for the Cannon to fit.

        :type spectrum:
            Spectrum

        :return:
        """

        assert isinstance(spectrum, fourgp_speclib.Spectrum), \
            "Supplied spectrum for the Cannon to fit is not a Spectrum object."

        assert spectrum.raster_hash == self._training_set.raster_hash, \
            "Supplied spectrum for the Cannon to fit is not sampled on the same raster as the training set."

        if self._debugging:
            self._debugging_output_counter += 1

        # Fitting tolerances
        max_iterations = 20  # Iterate a maximum number of times

        # Work out the raster of pixels inside each wavelength arm
        raster = spectrum.wavelengths
        lower_cut = 0
        arm_rasters = []
        for break_point in self._wavelength_arms:
            arm_rasters.append((raster >= lower_cut) * (raster < break_point))
            lower_cut = break_point
        arm_rasters.append(raster >= lower_cut)

        # Make initial continuum mask, which covers entire spectrum
        continuum_mask = np.ones_like(raster, dtype=bool)

        # Begin iterative fitting of continuum
        iteration = 0
        while True:
            iteration += 1

            # Treat each wavelength arm separately.
            continuum_models = []
            for i, arm_raster in enumerate(arm_rasters):
                # Make a mask of pixels which are both continuum and inside this wavelength arm
                pixel_mask = (arm_raster * continuum_mask *
                              np.isfinite(spectrum.value_errors) * (spectrum.value_errors > 0)
                              )
                # logger.info("Continuum pixels in arm {}: {} / {}".format(i, sum(pixel_mask), len(pixel_mask)))
                continuum_raster = raster[pixel_mask]
                continuum_values = spectrum.values[pixel_mask]
                continuum_value_errors = spectrum.value_errors[pixel_mask]

                # Make a new spectrum object containing only continuum pixels inside this wavelength arm
                continuum_spectrum = fourgp_speclib.Spectrum(wavelengths=continuum_raster,
                                                             values=continuum_values,
                                                             value_errors=continuum_value_errors,
                                                             )
                # logger.info("Continuum spectrum length: {}".format(len(continuum_spectrum)))

                # Fit a smooth function through these pixels
                continuum_model_factory = fourgp_speclib.SpectrumSmoothFactory(
                    function_family=self._continuum_model_family,
                    wavelengths=continuum_raster)

                continuum_smooth = continuum_model_factory.fit_to_continuum_via_mask(
                    other=continuum_spectrum,
                    mask=np.ones_like(continuum_raster, dtype=bool)
                )

                if isinstance(continuum_smooth, str):
                    logger.info(continuum_smooth)
                    return None, None, None, None, None, None

                # logger.info("Best-fit polynomial coefficients: {}".format(continuum_smooth.coefficients))

                # Resample smooth function onto the full raster of pixels within this wavelength arm
                resampler = SpectrumResampler(input_spectrum=continuum_smooth)
                continuum_models.append(resampler.onto_raster(raster[arm_raster]))

            # Splice together the continuum in all the wavelength arms
            continuum_model = fourgp_speclib.spectrum_splice(*continuum_models)

            # Create continuum-normalised spectrum using the continuum model we've just made
            cn_spectrum = spectrum / continuum_model

            # Run the Cannon
            labels, cov, meta = super(CannonInstanceCaseyOldWithContinuumNormalisation, self). \
                fit_spectrum(spectrum=cn_spectrum)

            # Fetch the Cannon's model spectrum
            model = fourgp_speclib.Spectrum(wavelengths=raster,
                                            values=self._model.predict(labels=labels),
                                            value_errors=np.zeros_like(raster))

            # Make new model of which pixels are continuum (based on Cannon's template being close to one)
            continuum_mask = (model.values > 0.99) * (model.values < 1.01)
            logger.info("Continuum pixels: {} / {}".format(sum(continuum_mask), len(continuum_mask)))
            logger.info("Best-fit labels: {}".format(list(labels[0])))

            # Produce debugging output if requested
            if self._debugging:
                np.savetxt("/tmp/debug_{:06d}_{:03d}.txt".format(self._debugging_output_counter, iteration),
                           np.transpose([raster, spectrum.values, spectrum.value_errors,
                                         continuum_model.values, model.values, continuum_mask])
                           )

            # Decide whether output is good enough for us to stop iterating
            if iteration >= max_iterations:
                break

        return labels, cov, meta
