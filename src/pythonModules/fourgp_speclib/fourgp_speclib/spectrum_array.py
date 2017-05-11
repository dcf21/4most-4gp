#!/usr/bin/env python
# -*- coding: utf-8 -*-

from os import path as os_path
import numpy as np
import logging
from ctypes import c_double
from multiprocessing.sharedctypes import RawArray

from spectrum import hash_numpy_array, Spectrum

logger = logging.getLogger(__name__)


class SpectrumArray(object):
    """
    An object representing an array of spectra.
    
    :ivar np.ndarray wavelengths:
        A 1D array listing the wavelengths at which this array of spectra are sampled.
        
    :ivar np.ndarray values:
        A 2D array listing the value measurements for each spectrum in this SpectrumArray.
        
    :ivar np.ndarray value_errors:
        A 2D array listing the standard errors in the value measurements for each spectrum in this SpectrumArray.
        
    :ivar str raster_hash:
        A string hash of the wavelength raster, used to quickly check whether spectra are sampled on a common raster.
        
    :ivar list[dict] metadata_list:
        A list of dictionaries of metadata about each of the spectra in this SpectrumArray
        
    :ivar bool shared_memory:
        Boolean flag indicating whether this SpectrumArray uses multiprocessing shared memory.
    """

    def __init__(self, wavelengths, values, value_errors, metadata_list, shared_memory=False):
        """
        Instantiate new SpectrumArray object.
        
        :param wavelengths: 
            A 1D array listing the wavelengths at which this array of spectra are sampled.
        
        :param values: 
            A 2D array listing the value measurements for each spectrum in this SpectrumArray.
            
        :param value_errors: 
            A 2D array listing the standard errors in the value measurements for each spectrum in this SpectrumArray.
            
        :param metadata_list:
            A list of dictionaries of metadata about each of the spectra in this SpectrumArray
            
        :param shared_memory:
            Boolean flag indicating whether the data in this SpectrumArray is supplied in multiprocessing shared memory.
            
        :type shared_memory:
            bool
        """

        # Sanity check inputs
        assert wavelengths.shape[0] == values.shape[1], "Inconsistent number of wavelength samples."
        assert wavelengths.shape[0] == value_errors.shape[1], "Inconsistent number of wavelength samples."
        assert values.shape[0] == value_errors.shape[0], "Inconsistent number of spectra in SpectrumArray."
        assert len(metadata_list) == values.shape[0], "Inconsistent number of spectra in SpectrumArray."

        # Store inputs as instance variables
        self.wavelengths = wavelengths
        self.values = values
        self.value_errors = value_errors
        self.metadata_list = metadata_list
        self.shared_memory = shared_memory
        self._update_raster_hash()

    def __len__(self):
        """
        Return the number of spectra in this SpectrumArray.
        
        :return:
            Integer number of spectra in this SpectrumArray.
        """
        return self.values.shape[0]

    @classmethod
    def from_files(cls, filenames, metadata_list, path="", shared_memory=False):
        """
        Instantiate new SpectrumArray object, using data in a list of text files.
        
        :param filenames: 
            List of the filenames of the text files from which to import spectra. Each file should have three columns:
            wavelength, value, and error in value.
            
        :type filenames:
            List[str]
            
        :param path:
            The file path from which to load the list of spectrum files.
            
        :type path:
            str
            
        :param metadata_list:
            A list of dictionaries of metadata about each of the spectra in this SpectrumArray
            
        :param shared_memory:
            Boolean flag indicating whether this SpectrumArray should use multiprocessing shared memory.
            
        :type shared_memory:
            bool
         
        :return:
            SpectrumArray object
        """

        # Load first spectrum to work out what wavelength raster we're using
        assert len(filenames) > 0, "Cannot open a SpectrumArray with no members: there is no wavelength raster"
        wavelengths, item_values, item_value_errors = np.loadtxt(os_path.join(path, filenames[0])).T
        raster_hash = hash_numpy_array(wavelengths)

        # Allocate numpy array to store this SpectrumArray into
        if not shared_memory:

            # If we're not using shared memory (which the multiprocessing module can share between threads),
            # we allocate a simple numpy array
            values = np.empty([len(filenames), len(wavelengths)])
            value_errors = np.empty([len(filenames), len(wavelengths)])
        else:

            # If we need to shared this array between threads (read only!), then we allocate the memory as a
            # multiprocessing RawArray
            wavelengths_shared_base = RawArray(c_double, wavelengths.size)
            wavelengths_shared = np.frombuffer(wavelengths_shared_base)
            wavelengths_shared[:] = wavelengths[:]
            wavelengths = wavelengths_shared

            values_shared_base = RawArray(c_double, wavelengths.size * len(filenames))
            values = np.frombuffer(values_shared_base)
            values = values.reshape([len(filenames), len(wavelengths)])

            value_errors_shared_base = RawArray(c_double, wavelengths.size * len(filenames))
            value_errors = np.frombuffer(value_errors_shared_base)
            value_errors = value_errors.reshape([len(filenames), len(wavelengths)])

        # Load spectra one by one
        for i, filename in enumerate(filenames):
            filename = os_path.join(path, filename)
            assert os_path.exists(filename), "File <{}> does not exist.".format(filename)
            item_wavelengths, item_values, item_value_errors = np.loadtxt(filename).T
            assert hash_numpy_array(item_wavelengths) == raster_hash, \
                "Item <{}> has a different wavelength raster from preceding spectra in SpectrumArray.".format(
                    filename)
            values[i, :] = item_values
            value_errors[i, :] = item_values

        # Instantiate a SpectrumArray object
        return cls(wavelengths=wavelengths,
                   values=values,
                   value_errors=value_errors,
                   metadata_list=metadata_list,
                   shared_memory=shared_memory)

    def __str__(self):
        return "<{module}.{name} instance".format(module=self.__module__,
                                                  name=type(self).__name__)

    def __repr__(self):
        return "<{0}.{1} object at {2}>".format(self.__module__,
                                                type(self).__name__, hex(id(self)))

    def _update_raster_hash(self):
        """
        Update the internal string hash of the wavelength raster that this spectrum array is sampled on.
        
        This hash is used to quickly check whether two spectra are sampled on the same raster before doing arithmetic
        operations on them.
        
        :return:
            None
        """
        self.raster_hash = hash_numpy_array(self.wavelengths)

    def get_metadata(self, index):
        """
        Extract the metadata which we have on a single spectrum from a SpectrumArray.
        
        :param index:
            Index of the spectrum to extract
            
        :type index:
            int
            
        :return:
            dict
        """
        return self.metadata_list[index]

    def extract_item(self, index):
        """
        Extract a single spectrum from a SpectrumArray. This creates a numpy view of the spectrum, without copying the
        data.
        
        :param index:
            Index of the spectrum to extract
            
        :type index:
            int
            
        :return:
            Spectrum object
        """

        # Check that requested index is within range
        assert 0 <= index < len(self), "Index of SpectrumArray out of range."

        return Spectrum(wavelengths=self.wavelengths,
                        values=self.values[index, :],
                        value_errors=self.value_errors[index, :],
                        metadata=self.metadata_list[index])
