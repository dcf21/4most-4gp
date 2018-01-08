#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from multiprocessing import cpu_count
import logging
import AnniesLasso as tc

import fourgp_speclib

logger = logging.getLogger(__name__)


# This is a drop-in replacement for astropy Tables. Astropy is a bit of a pain to install on lunarc because it needs
# cfitsio, so this saves a lot of trouble...
def dcf_table(names, rows):
    dtypes = [(str(i), np.float64) for i in names]

    # Numpy breaks with some incomprehensible error message "expected a readable buffer object" if it doesn't
    # receive a list of tuples...
    data = [tuple(i) for i in rows]

    output = np.array(data, dtype=dtypes)
    return output


class CannonInstance(object):
    """
    A class which holds an instance of the Cannon, and provides convenience methods for training it on arrays of spectra
    loaded from 4GP SpectrumLibrary objects.
    """

    def __init__(self, training_set, label_names, censors=None, progress_bar=False, threads=None, tolerance=1e-4,
                 load_from_file=None):
        """
        Instantiate the Cannon and train it on the spectra contained within a SpectrumArray.
        
        :param training_set:
            A SpectrumArray containing the spectra to train the Cannon on.

        :param label_names:
            A list of the names of the labels the Cannon is to estimate. We require that all of the training spectra
            have metadata fields defining all of these labels.

        :param threads:
            The number of CPU cores we should use. If None, we look up how many cores this computer has.

        :param tolerance:
            The tolerance xtol which the method <scipy.optimize.fmin_powell> uses to determine convergence.

        :param load_from_file:
            The filename of the internal state of a pre-trained Cannon, which we should load rather than doing
            training from scratch.
        """

        assert isinstance(training_set, fourgp_speclib.SpectrumArray), \
            "Training set for the Cannon should be a SpectrumArray."
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
                    label, index, ", ".join(metadata.keys()))
                assert np.isfinite(metadata[label]), "Label <{}> is not finite on training spectrum number {}. " \
                                          "Labels on this spectrum are: {}.".format(
                    label, index, metadata)

        # Compile table of training values of labels from metadata contained in SpectrumArray
        training_label_values = dcf_table(names=label_names,
                                          rows=[[training_set.get_metadata(index)[label] for label in label_names]
                                                for index in range(len(training_set))])

        self._model = tc.L1RegularizedCannonModel(labelled_set=training_label_values,
                                                  normalized_flux=training_set.values,
                                                  normalized_ivar=inverse_variances,
                                                  dispersion=training_set.wavelengths,
                                                  threads=threads)

        self._model.vectorizer = tc.vectorizer.NormalizedPolynomialVectorizer(
            labelled_set=training_label_values,
            terms=tc.vectorizer.polynomial.terminator(label_names, 2)
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
        Fit stellar labels to a spectrum.
        
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

    def fit_spectrum_with_continuum(self, spectrum, wavelength_arms, continuum_model_factory):
        """
        Fit stellar labels to a spectrum which has not been continuum normalised.

        :param spectrum:
            A Spectrum object containing the spectrum for the Cannon to fit.

        :type spectrum:
            Spectrum

        :param wavelength_arms:
            A list of the wavelength break-points between arms which should have continuum fitted separately.

        :type wavelength_arms:
            List of wavelengths in A.

        :param continuum_model_factory:
            A factory for smooth functions of some form, which we use to fit the continuum in this spectrum.

        :type continuum_model_factory:
            SpectrumSmoothFactory instance

        :return:
        """

        from fourgp_speclib import Spectrum, SpectrumSmoothFactory, spectrum_splice
        from fourgp_degrade import SpectrumInterpolator

        assert isinstance(spectrum, fourgp_speclib.Spectrum), \
            "Supplied spectrum for the Cannon to fit is not a Spectrum object."

        assert spectrum.raster_hash == self._training_set.raster_hash, \
            "Supplied spectrum for the Cannon to fit is not sampled on the same raster as the training set."

        assert isinstance(continuum_model_factory, SpectrumSmoothFactory), \
            "Input continuum model factory must be an instance of <SpectrumSmoothFactory>."

        # Fitting tolerances
        max_iterations = 3  # Iterate a maximum of 3 times

        # Work out the raster of pixels inside each wavelength arm
        raster = spectrum.wavelengths
        lower_cut = 0
        arm_rasters = []
        for break_point in wavelength_arms:
            arm_rasters.append((raster>=lower_cut) * (raster < break_point))
            lower_cut = break_point
        arm_rasters.append(raster>=lower_cut)

        # Make initial continuum mask, which covers entire spectrum
        continuum_mask = np.ones_like(raster)

        # Begin iterative fitting
        iteration = 0
        while True:
            iteration += 1

            continuum_models = []
            for i, arm_raster in enumerate(arm_rasters):
                pixel_mask = arm_raster[i] * continuum_mask
                continuum_raster = raster[pixel_mask]
                continuum_values = spectrum.values[pixel_mask]

                continuum_spectrum = Spectrum(wavelengths=continuum_raster, values=continuum_values)
                continuum_smooth = continuum_model_factory.fit_to_continuum_via_mask(other=continuum_spectrum,
                                                                                     mask=np.ones_like(continuum_raster)
                                                                                     )
                continuum_models.append(continuum_smooth)
            continuum_model = spectrum_splice(*continuum_models)

            # Create continuum-normalised spectrum
            cn_spectrum = spectrum / continuum_model

            # Run the Cannon
            labels, cov, meta = self.fit_spectrum(spectrum=cn_spectrum)

            # Fetch the Cannon's model spectrum
            model = Spectrum(wavelengths=raster,
                             values=self._model.predict(labels=labels),
                             value_errors=np.zeros_like(raster))

            # Make new model of which pixels are continuum
            continuum_mask = (model.values>0.99) * (model.values<1.01)

            # Decide whether output is good enough for us to stop iterating
            if iteration<max_iterations:
                break

            # Refine continuum model
            continuum = continuum_model * difference_spectrum
            lower_cut = 0
            arm_continuum = []
            for break_point in wavelength_arms:



        return labels, cov, meta

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
