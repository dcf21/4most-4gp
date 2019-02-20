# -*- coding: utf-8 -*-

"""
This provides a wrapper to Anna Ho's Cannon.

The code for this Cannon is in the master branch of this repository: https://github.com/annayqho/TheCannon

This is a version of the Cannon that Dominic Ford began investigating in November 2018.
"""

import logging
import pickle

import numpy as np

try:
    import TheCannon.dataset as ho_dataset
    import TheCannon.model as ho_model
except ImportError:
    print("!!! Warning! Could not find Anna Ho's Cannon in python environment.")

import fourgp_speclib

logger = logging.getLogger(__name__)


class CannonInstanceAnnaHo(object):
    """
    THIS CLASS IS CURRENTLY DEPRECATED (MARCH 2019). USE <CannonInstanceCaseyOld> FOR BEST RESULTS.

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

        assert polynomial_order == 2, "Anna Ho's Cannon only supports quadratic polynomials. " \
                                      "You requested <polynomial_order={}>.".format(polynomial_order)

        assert censors == None, "Anna Ho's Cannon does not support censoring. " \
                                "But you requested that it should be enabled."

        self._debugging_output_counter = 0
        self._debugging = debugging
        self.cannon_version = "AnnaHo"
        self._label_names = label_names
        self._wavelength_arms = wavelength_arms
        logger.info("Wavelength arm breakpoints: {}".format(self._wavelength_arms))

        assert isinstance(training_set, fourgp_speclib.SpectrumArray), \
            "Training set for the Cannon should be a SpectrumArray."

        # Hook for normalising input spectra
        training_set = self.normalise(training_set)

        self._training_set = training_set

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
        dataset = ho_dataset.Dataset(wl=training_set.wavelengths,
                                     tr_ID=range(len(training_set)),
                                     tr_flux=training_set.values,
                                     tr_ivar=inverse_variances,
                                     tr_label=np.array([np.array([training_set.get_metadata(index)[label]
                                                                  for label in label_names])
                                                        for index in range(len(training_set))]),
                                     test_ID=[],
                                     test_flux=[],
                                     test_ivar=[]
                                     )

        dataset.set_label_names(names=label_names)
        self._model = ho_model.CannonModel(order=2, useErrors=False)

        if load_from_file is None:
            logger.info("Starting to train the Cannon")
            self._model.train(ds=dataset)
            logger.info("Cannon training completed")
        else:
            logger.info("Loading Cannon from disk")
            self._model = pickle.load(file=open(load_from_file, "rb"))
            logger.info("Cannon loaded successfully")

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

        # Compile table of training values of labels from metadata contained in SpectrumArray
        dataset = ho_dataset.Dataset(wl=spectrum.wavelengths,
                                     tr_ID=[],
                                     tr_flux=[],
                                     tr_ivar=[],
                                     tr_label=[],
                                     test_ID=np.array((0,)),
                                     test_flux=np.array((spectrum.values,)),
                                     test_ivar=np.array((inverse_variances,))
                                     )

        dataset.set_label_names(names=self._label_names)
        errs_all, chisq_all = self._model.infer_labels(ds=dataset)

        labels = dataset.test_label_vals
        cov = errs_all
        meta = None

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

    def save_model(self, filename, overwrite=None):
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
        return
        pickle.dump(
            obj=self._model,
            file=open(filename, "wb")
        )

    def __str__(self):
        return "<{module}.{name} instance".format(module=self.__module__,
                                                  name=type(self).__name__)

    def __repr__(self):
        return "<{0}.{1} object at {2}>".format(self.__module__,
                                                type(self).__name__, hex(id(self)))
