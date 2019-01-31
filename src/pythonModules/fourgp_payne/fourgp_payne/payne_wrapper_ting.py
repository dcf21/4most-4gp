# -*- coding: utf-8 -*-

"""
This provides a wrapper to Yuan-Sen Ting's version of the Payne.

The code for this Payne is as sent by email from Yuan-Sen Ting, 22/01/2019

"""

import logging
import os
import pickle
from multiprocessing import cpu_count

import fourgp_speclib
import numpy as np

from .ting_train_nn import train_nn
from .ting_test_nn import test_nn

logger = logging.getLogger(__name__)


class PayneInstanceTing(object):
    """
    A class which holds an instance of the Payne, and provides convenience methods for training it on arrays of spectra
    loaded from 4GP SpectrumLibrary objects.
    """

    def __init__(self, training_set, label_names, training_data_archive,
                 censors=None, threads=None, batch_number=0, batch_count=1,
                 debugging=False):
        """
        Instantiate the Payne and train it on the spectra contained within a SpectrumArray.

        :param training_set:
            A SpectrumArray containing the spectra to train the Payne on.

        :param label_names:
            A list of the names of the labels the Payne is to estimate. We require that all of the training spectra
            have metadata fields defining all of these labels.

        :param threads:
            The number of CPU cores we should use. If None, we look up how many cores this computer has.

        :param batch_number:
            If training pixels in multiple batches on different machines, then this is the number of the batch of pixels
            we are to train. It should be in the range 0 .. batch_count-1 inclusive.

        :param batch_count:
            If training pixels in multiple batches on different machines, then this is the number of batches.

        :param training_data_archive:
            The filename of a directory where we store the internal state of a pre-trained neural network, which we
            reload when we want to test the Payne. Each batch of pixels stores its weights in a separate file, so that
            we can bring together the trained pixels from many different nodes of a cluster.

        :param debugging:
            Boolean flag determining whether we produce debugging output

        :type debugging:
            bool
        """

        self._debugging_output_counter = 0
        self._debugging = debugging
        self.payne_version = "Ting"
        self._label_names = label_names

        assert isinstance(training_set, fourgp_speclib.SpectrumArray), \
            "Training set for the Payne should be a SpectrumArray."

        # Hook for normalising input spectra
        training_set = self.normalise(training_set)

        self._training_set = training_set

        # Work out how many CPUs we should allow the Cannon to use
        if threads is None:
            threads = cpu_count()
        self.threads = threads

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
                                          "Labels on this spectrum are: {}.".\
                    format(label, index, ", ".join(list(metadata.keys())))
                assert np.isfinite(metadata[label]), "Label <{}> is not finite on training spectrum number {}. " \
                                                     "Labels on this spectrum are: {}.".\
                    format(label, index, metadata)

        # Compile table of training values of labels from metadata contained in SpectrumArray
        training_label_values = np.array([[training_set.get_metadata(index)[label] for label in label_names]
                                           for index in range(len(training_set))])

        # If we need to train a batch of pixels, do that now
        if batch_number >= 0:
            logger.info("Starting to train the Payne")
            training_data_batch = train_nn(
                batch_number=batch_number,
                batch_count=batch_count,
                labelled_set=training_label_values,
                normalized_flux=training_set.values,
                normalized_ivar=inverse_variances,
                dispersion=training_set.wavelengths
            )

            # Save weights of neural network
            logger.info("Payne training completed")
            os.system("mkdir -p {}".format(training_data_archive))

            logger.info("Saving Payne to disk")
            batch_pickle_filename = os.path.join(training_data_archive,
                                                 "batch_{:04d}_of_{:04d}.pkl".format(batch_number, batch_count))
            with open(batch_pickle_filename, "wb") as f:
                pickle.dump(training_data_batch, f)
            logger.info("Saving completed")

        # Reload training data from all batches in preparation for testing
        logger.info("Loading Payne from disk")
        payne_batches = []
        for i in range(batch_count):
            batch_pickle_filename = os.path.join(training_data_archive,
                                                 "batch_{:04d}_of_{:04d}.pkl".format(i, batch_count))

            assert os.path.exists(batch_pickle_filename), "Could not proceed from training to testing, because " \
                                                          "training data for batch {:d} of pixels is not present " \
                                                          "on this server.".format(i)

            payne_batches.append(pickle.load(open(batch_pickle_filename, 'rb')))
        logger.info("Payne loaded successfully")

        self._payne_status = {}
        for keyword in payne_batches[0]:
            self._payne_status[keyword] = np.concatenate(
                [payne_batches[i][keyword] for i in range(len(payne_batches))]
            )
        del payne_batches  # Free up memory

        logger.info("Payne batches merged successfully")

    def fit_spectrum(self, spectrum):
        """
        Fit stellar labels to a continuum-normalised spectrum.

        :param spectrum:
            A Spectrum object containing the spectrum for the Payne to fit.

        :type spectrum:
            Spectrum

        :return:
        """

        assert isinstance(spectrum, fourgp_speclib.Spectrum), \
            "Supplied spectrum for the Payne to fit is not a Spectrum object."

        assert spectrum.raster_hash == self._training_set.raster_hash, \
            "Supplied spectrum for the Payne to fit is not sampled on the same raster as the training set."

        # Hook for normalising input spectra
        spectrum = self.normalise(spectrum)

        inverse_variances = spectrum.value_errors ** (-2)

        # Ignore bad pixels.
        bad = (spectrum.value_errors < 0) + (~np.isfinite(inverse_variances * spectrum.values))
        inverse_variances[bad] = 0
        spectrum.values[bad] = np.nan

        fit_data = test_nn(
            payne_status=self._payne_status,
            threads=self.threads,
            num_labels=len(self._label_names),
            test_spectra=np.array([spectrum.values]),  # We are testing spectra one at a time, but wrapper allowing parallel testing
            test_spectra_errors=np.array([inverse_variances]),
        )

        return fit_data

    def normalise(self, spectrum):
        """
        This is a hook for doing some kind of normalisation on spectra. Not implemented in this base class.

        :param spectrum:
            The spectrum to be normalised.
        :return:
            Normalised version of this spectrum.
        """

        return spectrum

    def __str__(self):
        return "<{module}.{name} instance".format(module=self.__module__,
                                                  name=type(self).__name__)

    def __repr__(self):
        return "<{0}.{1} object at {2}>".format(self.__module__,
                                                type(self).__name__, hex(id(self)))
