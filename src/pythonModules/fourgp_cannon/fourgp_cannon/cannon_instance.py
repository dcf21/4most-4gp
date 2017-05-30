#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from multiprocessing import cpu_count
import logging
from astropy.table import Table
import AnniesLasso as tc

import fourgp_speclib

logger = logging.getLogger(__name__)


class CannonInstance(object):
    """
    A class which holds an instance of the Cannon, and provides convenience methods for training it on arrays of spectra
    loaded from 4GP SpectrumLibrary objects.
    """

    def __init__(self, training_set, label_names, censors=None, progress_bar=False, threads=None):
        """
        Instantiate the Cannon and train it on the spectra contained within a SpectrumArray.
        
        :param training_set:
            A SpectrumArray containing the spectra to train the Cannon on.
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
            terms=tc.vectorizer.polynomial.terminator(label_names, 2)
        )

        if censors is not None:
            self._model.censors = censors

        self._model.s2 = 0
        self._model.regularization = 0

        logger.info("Starting to train the Cannon")
        self._model.train(progressbar=self._progress_bar)
        logger.info("Cannon training completed")
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
        self._model.save(filename=filename, overwrite=overwrite)

    def __str__(self):
        return "<{module}.{name} instance".format(module=self.__module__,
                                                  name=type(self).__name__)

    def __repr__(self):
        return "<{0}.{1} object at {2}>".format(self.__module__,
                                                type(self).__name__, hex(id(self)))
