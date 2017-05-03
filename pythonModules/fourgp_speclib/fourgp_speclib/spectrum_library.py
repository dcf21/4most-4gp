#!/usr/bin/env python
# -*- coding: utf-8 -*-


def requires_ids_or_filenames(method):
    """
    A decorator for spectrum library methods that require either a list of Ids or a list of filenames.

    :param method:
        A method belonging to a sub-class of SpectrumLibrary.
    """

    def wrapper(model, *args, **kwargs):
        have_ids = ("ids" in kwargs) and (kwargs["ids"] is not None)
        have_filenames = ("filenames" in kwargs) and (kwargs["filenames"] is not None)
        assert have_ids or have_filenames, "Must supply a list of Ids or a list of filenames"
        assert not (have_ids and have_filenames), "Must supply either a list of Ids or a list of filenames, not both."

        # If a single Id is supplied, rather than a list of Ids, turn it into a one-entry tuple
        if have_ids and not isinstance(kwargs["ids"], (list, tuple)):
            kwargs["ids"] = (kwargs["ids"],)

        # If a single filename is supplied, turn it into a one-entry tuple
        if have_filenames and not isinstance(kwargs["filenames"], (list, tuple)):
            kwargs["filenames"] = (kwargs["filenames"],)

        return method(model, *args, **kwargs)

    return wrapper


class SpectrumLibrary(object):
    """
    An abstract spectrum library object
    """

    def __init__(self):
        self._metadata_fields = None

    def __str__(self):
        return "<{module}.{name} instance".format(module=self.__module__,
                                                  name=type(self).__name__)

    def __repr__(self):
        return "<{0}.{1} object at {2}>".format(self.__module__,
                                                type(self).__name__, hex(id(self)))

    def purge(self):
        """
        This irrevocably deletes the spectrum library from the database and from your disk. You have been warned.
         
        :return:
            None
        """

        raise NotImplementedError("The purge method must be implemented by each SpectrumLibrary implementation.")

    def list_metadata_fields(self):
        """
        List all of the metadata fields set on spectra in this spectrum library.
        
        :return:
            List of strings
        """

        return self._metadata_fields

    def search(self, **kwargs):
        """
        Search for spectra within this SpectrumLibrary which fall within some metadata constraints.
        
        :param kwargs:
            A dictionary of metadata constraints. Constraints can be specified either as <key: value> pairs, in
            which case the value must match exactly, or as <key: [min,max]> in which case the value must fall within
            the specified range.
         
        :return:
            A tuple of objects, each representing a spectrum which matches the search criteria. Within each object,
            the properties <specId> and <filename> are defined as integers and strings respectively.
        """

        raise NotImplementedError("The search method must be implemented by each SpectrumLibrary implementation.")

    def open(self, ids=None, filenames=None):
        """
        Open some spectra from this spectrum library, and return them as a SpectrumArray object.
        
        :param ids: 
            List of the integer ids of the spectra to receive this metadata, or None to select them by filename.
            
        :type ids:
            List of int, or None
            
        :param filenames:
            List of the filenames of the spectra to receive this metadata, or None to select them by integer id.
            
        :type filenames:
            List of str, or None
            
        :return:
            A SpectrumArray object.
        """

        raise NotImplementedError("The open method must be implemented by each SpectrumLibrary implementation.")

    def insert(self, spectra, filenames, origin="Undefined", metadata=None, overwrite=False):
        """
        Insert the spectra from a SpectrumArray object into this spectrum library.
        
        :param spectra: 
            A SpectrumArray object contain the spectra to be inserted into this spectrum library.
            
        :type spectra:
            SpectrumArray
            
        :param filenames:
            A list of the filenames with which to save the spectra contained within this SpectrumArray
            
        :type filenames:
            List of str
            
        :param origin: 
            A string describing where these spectra are being imported from. Normally the name of the module which is
            importing them.
            
        :type origin:
            str
            
        :param metadata: 
            A list of dictionaries of metadata to set on each of the spectra in this SpectrumArray.
            
        :type metadata:
            List of dict
            
        :param overwrite:
            Boolean flag indicating whether we're allowed to overwrite pre-existing spectra with the same filenames
            
        :type overwrite:
            bool
            
        :return:
            None
        """

        raise NotImplementedError("The insert method must be implemented by each SpectrumLibrary implementation.")

    def import_from(self, other, overwrite=False, **kwargs):
        """
        Search for spectra within another SpectrumLibrary, and import all matching spectra into this library.
        
        :param other:
            The SpectrumLibrary from which we should import spectra.
            
        :type other:
            SpectrumLibrary
        
        :param overwrite:
            Boolean flag indicating whether we're allowed to overwrite pre-existing spectra with the same filenames
            
        :type overwrite:
            bool
            
        :param kwargs:
            A dictionary of metadata constraints. Constraints can be specified either as <key: value> pairs, in
            which case the value must match exactly, or as <key: [min,max]> in which case the value must fall within
            the specified range.
         
        :return:
            A tuple of objects, each representing a spectrum which matches the search criteria. Within each object,
            the properties <specId> and <filename> are defined as integers and strings respectively.
        """

        spectra = other.search(**kwargs)

        for spectrum in spectra:
            uid = spectrum["specId"]
            filename = spectrum["filename"]
            origin = spectrum["origin"]
            obj = other.open(ids=[uid])
            self.insert(spectra=obj, filenames=[filename], origin=origin, overwrite=overwrite)
