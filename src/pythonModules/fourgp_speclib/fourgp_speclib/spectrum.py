# -*- coding: utf-8 -*-

from os import path as os_path
from math import sqrt
import numpy as np
import scipy.integrate
import hashlib
import logging

logger = logging.getLogger(__name__)


def hash_numpy_array(item):
    """
    Efficiently produce a string hash of a numpy array, for quick checking of whether arrays match.
    
    :param item:
        Any numpy array
        
    :type item:
        np.ndarray
        
    :return:
        String hash
    """
    raw = item.copy().view(np.uint8)
    return hashlib.sha1(raw).hexdigest()


def requires_common_raster(method):
    """
    A decorator for spectrum methods that require that another spectrum as an input and require it to be sampled on the
    same wavelength raster as us.

    :param method:
        A method belonging to a sub-class of Spectrum.
    """

    def wrapper(spectrum, other, *args, **kwargs):
        assert isinstance(other, Spectrum), "Can only copy mask from another Spectrum object."
        assert spectrum.raster_hash == other.raster_hash, \
            "Cannot do arithmetic on spectra sampled on a different wavelength rasters"
        return method(spectrum, other=other, *args, **kwargs)

    return wrapper


def spectrum_splice(*spectra):
    """
    Splice together a number of input spectra.

    :param spectra:
        A list of spectra to be spliced together.

    :return:
        Spliced Spectrum object
    """

    wavelengths = np.concatenate([x.wavelengths for x in spectra])
    values = np.concatenate([x.values for x in spectra])
    value_errors = np.concatenate([x.value_errors for x in spectra])

    return Spectrum(wavelengths=wavelengths,
                    values=values,
                    value_errors=value_errors)


class Spectrum(object):
    """
    An object representing a single spectrum.
    
    :ivar np.ndarray wavelengths:
        A 1D array listing the wavelengths at which this array of spectra are sampled.
        
    :ivar np.ndarray values:
        A 1D array listing the value measurements in this spectrum.
        
    :ivar np.ndarray value_errors:
        A 1D array listing the standard errors in the value measurements.
        
    :ivar np.ndarray mask:
        A 1D array listing which wavelength samples we've currently selected to use.
        
    :ivar bool mask_set:
        Boolean selecting whether any wavelengths are currently masked out.
        
    :ivar dict metadata:
        Dictionary of metadata about this spectrum.
        
    :ivar str raster_hash:
        A string hash of the wavelength raster, used to quickly check whether spectra are sampled on a common raster.
    """

    def __init__(self, wavelengths, values, value_errors, metadata=None):
        """
        Instantiate a new Spectrum object.
        
        :param wavelengths: 
            A 1D array listing the wavelengths at which this array of spectra are sampled.
            
        :type wavelengths:
            np.ndarray
        
        :param values: 
            A 2D array listing the value measurements for each spectrum in this SpectrumArray.
            
        :type values:
            np.ndarray
            
        :param value_errors: 
            A 2D array listing the standard errors in the value measurements for each spectrum in this SpectrumArray.
            
        :type value_errors:
            np.ndarray
            
        :param metadata:
            Dictionary of metadata about this spectrum.
            
        :type metadata:
            dict
        
        """

        if metadata is None:
            metadata = {}

        self.wavelengths = wavelengths
        self.values = values
        self.value_errors = value_errors
        self.mask = np.ones_like(self.wavelengths, dtype=np.bool)
        self.mask_set = False
        self.metadata = metadata

        self._validate_data_dimensions()
        self._update_raster_hash()

    @classmethod
    def from_file(cls, filename, binary=True, columns=None, *args, **kwargs):
        """
        Factory method to read a spectrum in from a text file.
        
        :param filename: 
            Filename of text file, containing wavelengths, values and errors in three columns
            
        :type filename:
            str

        :param binary:
            Boolean specifying whether we store spectra on disk in binary format or plain text.

        :type binary:
            bool

        :param columns:
            When reading text files, this specifies which columns the wavelength, value and value_error are in.
            The default is (0, 1, 2).

        :type columns:
            list or tuple
            
        :return:
            Spectrum object
        """

        assert os_path.exists(filename), "File <{}> does not exist.".format(filename)

        # Extract spectrum from a text file
        if not binary:
            # Work out which columns we are reading data from
            if columns is None:
                columns = (0, 1, 2)
            assert isinstance(columns, (list, tuple)), "Columns must be specified as a list or a tuple."
            assert len(columns) > 0, "At least one column must be specified to give us a wavelength raster."

            # Read the contents of the text file
            data = np.loadtxt(filename).T

            # Read wavelengths from text file
            wavelengths = data[columns[0]]

            # Read values from text file, or fill with zeros if they are not specified
            if len(columns) > 1:
                values = data[columns[1]]
            else:
                values = np.zeros_like(wavelengths)

            # Read value errors from text file, or fill with zeros if they are not specified
            if len(columns) > 2:
                value_errors = data[columns[2]]
            else:
                value_errors = np.zeros_like(wavelengths)

        # Extract spectrum from a binary file
        else:
            wavelengths, values, value_errors = np.load(filename)

        return cls(wavelengths=wavelengths, values=values, value_errors=value_errors, *args, **kwargs)

    def to_file(self, filename, binary=True, overwrite=False):
        """
        Dump a spectrum object to a text file, with three columns containing wavelengths, data values, and errors.
        
        :param filename:
            Filename of text file to write.
            
        :type filename:
            str
            
        :param binary:
            Boolean specifying whether we store spectra on disk in binary format or plain text.
            
        :type binary:
            bool

        :param overwrite:
            Boolean specifying whether we're allowed to overwrite existing files.

        :type overwrite:
            bool
            
        :return:
            bool: Success
        """

        if os_path.exists(filename) and not overwrite:
            logger.error("File <{}> already exists. Set overwrite API option to force overwriting of it.".format(
                filename))
            return False

        if not binary:
            np.savetxt(filename, np.transpose([self.wavelengths, self.values, self.value_errors]))
        else:
            np.save(filename, np.asarray([self.wavelengths, self.values, self.value_errors]))
        return True

    def __str__(self):
        return "<{module}.{name} instance".format(module=self.__module__,
                                                  name=type(self).__name__)

    def __repr__(self):
        return "<{0}.{1} object at {2}>".format(self.__module__,
                                                type(self).__name__, hex(id(self)))

    def __len__(self):
        return self.wavelengths.shape[0]

    def __eq__(self, other):
        return np.array_equal(self.wavelengths, other.wavelengths) and np.array_equal(self.values, other.values)

    def __ne__(self, other):
        return not self.__eq__(other)

    def _validate_data_dimensions(self):
        """
        Validate that the wavelength raster, value array, and error array are of matching sizes.
        
        :return: 
            None
        """
        for i in [self.wavelengths, self.values, self.value_errors]:
            assert isinstance(i, np.ndarray), "Input argument to Spectrum class was not a numpy array"

        for i in [self.wavelengths, self.values, self.value_errors]:
            assert len(i.shape) == 1, "Input argument to Spectrum class was not a 1D numpy array"
            assert i.shape[0] == self.wavelengths.shape[0], "Input argument to Spectrum class were of differing lengths"

    def _update_raster_hash(self):
        """
        Update the internal string hash of the wavelength raster that this spectrum array is sampled on.
        
        This hash is used to quickly check whether two spectra are sampled on the same raster before doing arithmetic
        operations on them.
        
        :return:
            None
        """
        self.raster_hash = hash_numpy_array(self.wavelengths)

    def copy(self):
        """
        Duplicate a spectrum, allocating new memory to hold a new copy of the data within in.

        :return:
            A new Spectrum object
        """

        new_wavelengths = self.wavelengths.copy()
        new_values = self.values.copy()
        new_value_errors = self.value_errors.copy()
        output = Spectrum(wavelengths=new_wavelengths, values=new_values, value_errors=new_value_errors,
                          metadata=self.metadata.copy())

        if self.mask_set:
            output.copy_mask_from(self)

        return output

    def integral(self):
        """
        Evaluate the integrated flux in a Spectrum

        :return:
            float
        """

        return scipy.integrate.trapz(x=self.wavelengths, y=self.values)

    def photometry(self, band):
        """
        Evaluate photometry in a named photometric band, using pyphot.

        :param band:
            Name of photometric band, as recognised by pyphot. See
            <http://mfouesneau.github.io/docs/pyphot/libcontent.html> for a list of recognised bands.

        :type band:
            str

        :return:
            float AB magnitude
        """

        import pyphot, tables
        lib = pyphot.get_library()

        try:
            photometer = lib[band]
        except tables.exceptions.NoSuchNodeError:
            logger.error("Could not find photometric band <{}>".format(band))
            raise

        wavelengths = self.wavelengths * pyphot.unit['AA']
        fluxes = self.values * pyphot.unit['erg/s/cm**2/AA']

        flux = photometer.get_flux(slamb=wavelengths, sflux=fluxes)
        mag = -2.5 * np.log10(flux) - photometer.AB_zero_mag
        return mag

    def apply_redshift(self, z):
        """
        Apply a redshift of z to this spectrum, and return a new Spectrum object.

        :param z:
            The redshift to apply.

        :type z:
            float

        :return:
            Spectrum object containing the redshifted spectrum.
        """

        new_wavelengths = self.wavelengths * (1 + z)
        new_values = self.values.copy()
        new_value_errors = self.value_errors.copy()
        output = Spectrum(wavelengths=new_wavelengths, values=new_values, value_errors=new_value_errors,
                          metadata=self.metadata.copy())

        if self.mask_set:
            output.copy_mask_from(self)

        return output

    def apply_radial_velocity(self, v):
        """
        Apply a radial velocity of v to this spectrum, and return a new Spectrum object. A positive radial velocity
        means that the object is receding from the observer.

        :param z:
            The radial velocity to apply (units m/s)

        :type z:
            float

        :return:
            Spectrum object containing the redshifted (receding) spectrum
        """
        # https://ned.ipac.caltech.edu/level5/Hogg/Hogg3.html
        c = 299792458.0
        return self.apply_redshift(sqrt((1 + v / c) / (1 - v / c)) - 1)

    def remove_redshift(self, z):
        """
        Undo a redshift of z from this spectrum, turning an observed (redshifted) spectrum into an object-rest-frame
        spectrum. The rest frame spectrum is returned as a new Spectrum object.

        :param z:
            The redshift to undo.

        :type z:
            float

        :return:
            Spectrum object containing the object-rest-frame spectrum
        """
        return self.apply_redshift(-z)

    def correct_radial_velocity(self, v):
        """
        Undo a radial velocity of v from this spectrum, turning an observed spectrum of an object receding at velocity v
        into an object-rest-frame spectrum. The rest frame spectrum is returned as a new Spectrum object. A positive
        radial velocity means that the object is receding from the observer.

        :param z:
            The radial velocity to undo (units m/s)

        :type z:
            float

        :return:
            Spectrum object containing the object-rest-frame spectrum
        """

        return self.apply_radial_velocity(-v)

    @requires_common_raster
    def copy_mask_from(self, other):
        """
        Duplicate the wavelength mask set up on another spectrum, and apply it to this one.
        
        :param other:
            The spectrum to copy the mask from
            
        :type other:
            Spectrum
        
        :return:
            None
        """
        self.mask = other.mask.copy()
        self.mask_set = other.mask_set

    def mask_include(self, wavelength_min=0, wavelength_max=np.inf):
        """
        Update the wavelength mask to include wavelengths in the specified range.
        
        :param wavelength_min: 
            The shortest wavelength in the range to be added into the mask.
            
        :type wavelength_min:
            float
        
        :param wavelength_max:
            The longest wavelength in the range to be added into the mask.
            
        :type wavelength_max:
            float
        
        :return: 
            None
        """

        window = (self.wavelengths >= wavelength_min) * (self.wavelengths <= wavelength_max)
        self.mask[window] = True
        self.mask_set = not np.all(self.mask)

    def mask_exclude(self, wavelength_min=0, wavelength_max=np.inf):
        """
        Update the wavelength mask to exclude wavelengths in the specified range.
        
        :param wavelength_min: 
            The shortest wavelength in the range to be removed from the mask.
            
        :type wavelength_min:
            float
        
        :param wavelength_max:
            The longest wavelength in the range to be removed from the mask.
            
        :type wavelength_max:
            float
        
        :return: 
            None
        """

        window = (self.wavelengths >= wavelength_min) * (self.wavelengths <= wavelength_max)
        self.mask[window] = False
        self.mask_set = not np.all(self.mask)

    def truncate_to_mask(self):
        """
        Return a new Spectrum object which contains only the data values within the masked region of this spectrum.

        :return:
            A new Spectrum object
        """

        new_wavelengths = self.wavelengths[self.mask]
        new_values = self.values[self.mask]
        new_value_errors = self.value_errors[self.mask]

        output = Spectrum(wavelengths=new_wavelengths, values=new_values, value_errors=new_value_errors)
        return output

    @requires_common_raster
    def __add__(self, other):
        """
        Add the values in another spectrum to the values in this one, and return a new Spectrum object.
        
        :param other:
            The Spectrum object to add to this one.
            
        :type other:
            Spectrum
        
        :return:
            Spectrum object containing the sum of the two spectra.
        """

        new_value_errors = np.hypot(self.value_errors, other.value_errors)
        new_values = self.values + other.values

        output = Spectrum(wavelengths=self.wavelengths, values=new_values, value_errors=new_value_errors)
        if self.mask_set or other.mask_set:
            output.mask = self.mask * other.mask  # Logical AND
            output.mask_set = True
        return output

    @requires_common_raster
    def __sub__(self, other):
        """
        Subtract the values in another spectrum from the values in this one, and return a new Spectrum object.
        
        :param other:
            The Spectrum object to subtract from this one.
            
        :type other:
            Spectrum
        
        :return:
            Spectrum object containing the difference of the two spectra.
        """

        new_value_errors = np.hypot(self.value_errors, other.value_errors)
        new_values = self.values - other.values

        output = Spectrum(wavelengths=self.wavelengths, values=new_values, value_errors=new_value_errors)
        if self.mask_set or other.mask_set:
            output.mask = self.mask * other.mask  # Logical AND
            output.mask_set = True
        return output

    @requires_common_raster
    def __iadd__(self, other):
        """
        Add the values in another spectrum to the values in this one.
        
        :param other:
            The Spectrum object to add to this one.
            
        :type other:
            Spectrum
        
        :return:
            self
        """

        self.value_errors = np.hypot(self.value_errors, other.value_errors)
        self.values = self.values + other.values

        if other.mask_set:
            self.mask *= other.mask  # Logical AND
            self.mask_set = True

        return self

    @requires_common_raster
    def __isub__(self, other):
        """
        Subtract the values in another spectrum from the values in this one.
        
        :param other:
            The Spectrum object to subtract from this one.
            
        :type other:
            Spectrum
        
        :return:
            self
        """

        self.value_errors = np.hypot(self.value_errors, other.value_errors)
        self.values = self.values - other.values

        if other.mask_set:
            self.mask *= other.mask  # Logical AND
            self.mask_set = True

        return self

    @requires_common_raster
    def __mul__(self, other):
        """
        Multiply the values in another spectrum by the values in this one, and return a new Spectrum object.
        
        :param other:
            The Spectrum object to multiply by this one.
            
        :type other:
            Spectrum
        
        :return:
            Spectrum object containing the sum of the two spectra.
        """

        new_values = self.values * other.values
        new_value_errors = np.hypot(self.value_errors / self.values, other.value_errors / other.values) * \
                           np.abs(new_values)

        output = Spectrum(wavelengths=self.wavelengths, values=new_values, value_errors=new_value_errors)
        if self.mask_set or other.mask_set:
            output.mask = self.mask * other.mask  # Logical AND
            output.mask_set = True
        return output

    @requires_common_raster
    def __div__(self, other):
        """
        Divide the values in this spectrum by the values in another one, and return a new Spectrum object.
        
        :param other:
            The Spectrum object to divide this one by.
            
        :type other:
            Spectrum
        
        :return:
            Spectrum object containing the quotient of the two spectra.
        """

        new_values = self.values / other.values
        new_value_errors = np.hypot(self.value_errors / self.values, other.value_errors / other.values) * \
                           np.abs(new_values)

        output = Spectrum(wavelengths=self.wavelengths, values=new_values, value_errors=new_value_errors)
        if self.mask_set or other.mask_set:
            output.mask = self.mask * other.mask  # Logical AND
            output.mask_set = True
        return output

    @requires_common_raster
    def __imul__(self, other):
        """
        Multiply the values in this spectrum by the values in another.
        
        :param other:
            The Spectrum object to multiply this one by.
            
        :type other:
            Spectrum
        
        :return:
            self
        """

        new_values = self.values * other.values
        self.value_errors = np.hypot(self.value_errors / self.values, other.value_errors / other.values) * \
                            np.abs(new_values)
        self.values = new_values

        if other.mask_set:
            self.mask *= other.mask  # Logical AND
            self.mask_set = True

        return self

    @requires_common_raster
    def __idiv__(self, other):
        """
        Divide the values in this spectrum by the values in another.
        
        :param other:
            The Spectrum object to divide this one by.
            
        :type other:
            Spectrum
        
        :return:
            self
        """

        new_values = self.values / other.values
        self.value_errors = np.hypot(self.value_errors / self.values, other.value_errors / other.values) * \
                            np.abs(new_values)
        self.values = new_values

        if other.mask_set:
            self.mask *= other.mask  # Logical AND
            self.mask_set = True

        return self
