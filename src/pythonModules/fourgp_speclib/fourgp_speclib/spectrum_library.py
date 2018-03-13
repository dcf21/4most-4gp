#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

from os import path as os_path
import re


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
    An abstract spectrum library object.
    
    Spectrum libraries are a bit like having a directory full of data files on disk, each containing a spectrum.
    However, they also include a database which can store arbitrary metadata about each spectrum -- for example,
    stellar parameters and abundances. It is possible to search a spectrum library based on metadata constraints.

    Various implementations of the SpectrumLibrary class are provided, storing the metadata in different flavours of
    SQL database. SQLite is probably the simplest and creates portable libraries that you can transfer to a different
    machine with all metadata intact. MySQL is a faster database engine, and probably a better option for data which
    doesn't need to move around.

    
    :ivar list[string] _metadata_fields:
        A list of the metadata fields set on spectra in this SpectrumLibrary
    """

    def __init__(self, path=None):
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

    @classmethod
    def open_and_search(cls, library_spec, workspace, extra_constraints):
        """
        Helper function which allows you to open a spectrum library and search it in a single function call.

        The argument <library_spec> can take the form of the name of a library <my_library>, or can contain
        metadata constraints as a comma-separated list in square brackets, e.g. <my_library[continuum_normalised=1]>.
        This is useful if you want to write a script which allows the user to specify parameter cuts on to command line
        when they specify which spectrum library to operate on.

        Multiple constraints should be specified as a comma-separated list, e.g.:

        my_library[continuum_normalised=1,5000<Teff<6000]

        Constraints can be specified in two formats as shown above. You can either require equality, or require that
        the parameter falls within a range using the < operator.

        If you use the < operator, you must specify both lower and upper limit.

        :param library_spec:
            The name of the spectrum library to open, suffixed with any metadata constraints in [] brackets.

        :param workspace:
            The path of disk to where spectrum libraries are stored.

        :param extra_constraints:
            A dictionary containing any additional metadata constraints to be added to the ones supplied in
            library_spec.

        :return:
            Dictionary, containing:
                library -- a SpectrumLibrary instance
                items -- the result of the metadata search within the SpectrumLibrary
                constraints -- a list of the metadata constraints applied
        """

        test = re.match("([^\[]*)\[(.*)\]$", library_spec)
        constraints = {}
        if test is None:
            library_name = library_spec
        else:
            library_name = test.group(1)
            for constraint in test.group(2).split(","):
                words_1 = constraint.split("=")
                words_2 = constraint.split("<")
                if len(words_1) == 2:
                    constraint_name = words_1[0]
                    try:
                        constraint_value = float(words_1[1])
                    except ValueError:
                        constraint_value = words_1[1]
                    constraints[constraint_name] = constraint_value
                elif len(words_2) == 3:
                    constraint_name = words_2[1]
                    try:
                        constraint_value_a = float(words_2[0])
                        constraint_value_b = float(words_2[2])
                    except ValueError:
                        constraint_value_a = words_2[0]
                        constraint_value_b = words_2[2]
                    constraints[constraint_name] = (constraint_value_a, constraint_value_b)
                else:
                    assert False, "Could not parse constraint <{}>".format(constraint)
        constraints.update(extra_constraints)
        library_path = os_path.join(workspace, library_name)
        input_library = cls(path=library_path)
        library_items = input_library.search(**constraints)
        return {
            "library": input_library,
            "items": library_items
        }
