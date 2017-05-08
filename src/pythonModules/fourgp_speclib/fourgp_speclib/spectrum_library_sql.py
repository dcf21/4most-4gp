#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from os import path as os_path
import re
import time
import hashlib
import logging

from spectrum_library import SpectrumLibrary, requires_ids_or_filenames
from spectrum_array import SpectrumArray
from spectrum import Spectrum

logger = logging.getLogger(__name__)


class SpectrumLibrarySql(SpectrumLibrary):
    """
    A spectrum library implementation that uses SQL database to store metadata about each spectrum.
    
    :cvar string _schema:
        The SQL schema used for storing SpectrumLibraries
        
    :ivar _db:
        Database handle object
        
    :ivar _db_cursor:
        Database cursor object
        
    :ivar _path:
        Path to the directory where this SpectrumLibrary is stored on disk
        
    :ivar _unique_id:
        The unique string identifier used to identify this SpectrumLibrary in the <libraries> table in the database
        
    :ivar _library_id:
        The numerical identifier for this SpectrumLibrary in the <libraries> table in the database
    """

    _schema = """

-- Table of spectrum libraries using this database
CREATE TABLE libraries (
    libraryId INTEGER PRIMARY KEY AUTO_INCREMENT,
    name VARCHAR(1024) UNIQUE NOT NULL
);

-- Table of string descriptions of which tools imported particular spectra into the library 
CREATE TABLE origins (
    originId INTEGER PRIMARY KEY AUTO_INCREMENT,
    name VARCHAR(256) UNIQUE NOT NULL
);

-- Table of spectra within this library
CREATE TABLE spectra (
    specId INTEGER PRIMARY KEY AUTO_INCREMENT,
    libraryId INTEGER NOT NULL,
    filename VARCHAR(256) UNIQUE NOT NULL,
    originId INTEGER NOT NULL,
    importTime REAL,
    FOREIGN KEY (libraryId) REFERENCES libraries (libraryId) ON DELETE CASCADE,
    FOREIGN KEY (originId) REFERENCES origins (originId) ON DELETE CASCADE
);

CREATE INDEX search_by_filename ON spectra (libraryId, filename);
CREATE INDEX search_by_id ON spectra (libraryId, specId);

-- Table of metadata fields which have been set on at least one spectrum
CREATE TABLE metadata_fields (
    fieldId INTEGER PRIMARY KEY AUTO_INCREMENT,
    name VARCHAR(32) UNIQUE NOT NULL
);

-- Table of all metadata items set of spectra
CREATE TABLE spectrum_metadata (
    specId INTEGER NOT NULL,
    fieldId INTEGER NOT NULL,
    libraryId INTEGER NOT NULL,  -- this field is also in the record in <spectra>, but copied here for fast searching
    valueFloat REAL,
    valueString VARCHAR(256),
    PRIMARY KEY (specId, fieldId),
    FOREIGN KEY (specId) REFERENCES spectra(specId) ON DELETE CASCADE,
    FOREIGN KEY (fieldId) REFERENCES metadata_fields(fieldId) ON DELETE CASCADE,
    FOREIGN KEY (libraryId) REFERENCES libraries (libraryId) ON DELETE CASCADE
);

CREATE INDEX search_metadata_floats ON spectrum_metadata (libraryId, fieldId, valueFloat);
CREATE INDEX search_metadata_strings ON spectrum_metadata (libraryId, fieldId, valueString);
    
    """

    def __init__(self, path, create=False):
        """
        Create a new SpectrumLibrary object, storing metadata about the spectra in an SQL database.
        
        :param path:
            The file path to use for storing spectra.
         
        :type path:
            str
         
        :param create:
            If true, create a new empty spectrum library. If false, it is an error to attempt to open a library which
            doesn't exist.
        
        :type create:
            bool
        """

        # Create new spectrum library if requested
        self._path = path
        if create:
            self._create()

        # Check that we're not overwriting an existing library
        assert os_path.exists(path), \
            "Could not open spectrum library <{}>: directory not found".format(path)

        if not create:
            self._db, self._db_cursor = self._open_database()

        try:
            # Check that this library uses the right flavour of SQL
            with open(os_path.join(self._path, "type_id")) as f:
                expected_library_type = type(self).__name__
                library_type = f.read()
                assert library_type == expected_library_type, \
                    "This library was created with class <%{}>. Cannot open with class <{}>.".format(
                        library_type, expected_library_type)

            # Look up numeric Id for this particular spectrum library in the database
            # Some SQL implementations (e.g. MySQL) share multiple spectrum libraries in a single database. Others
            # (e.g. SQLite) have a separate database for each spectrum library
            with open(os_path.join(self._path, "unique_id")) as f:
                self._unique_id = f.read().strip()
                self._library_id = self._fetch_library_id(self._unique_id, False)

        except IOError:
            logger.error("Spectrum library did not have required header files.")
            raise

        # Initialise
        self._metadata_init()

        super(SpectrumLibrarySql, self).__init__()

    def _create(self):
        """
        Create a new, empty spectrum library.
            
        :return:
            None
        """

        # Check that we're not overwriting an existing library
        assert not os_path.exists(self._path), \
            "Could not create spectrum library <{}>: file already exists".format(self._path)

        # Check that parent directory exists
        parent_path, file_name = os_path.split(self._path)
        assert os_path.exists(parent_path), \
            "Could not create spectrum library <{}>: parent directory does not exist".format(self._path)

        # Create directory to hold spectra in this library
        try:
            os.mkdir(self._path)
        except IOError:
            raise

        # Create SQL database to hold metadata about these spectra
        self._create_database()
        self._db, self._db_cursor = self._open_database()

        # Record the database format used by this library
        with open(os_path.join(self._path, "type_id"), "w") as f:
            library_type = type(self).__name__
            f.write(library_type)

        # Create a random unique id for this library
        with open(os_path.join(self._path, "unique_id"), "w") as f:
            unique_id = hashlib.md5(os.urandom(32).encode("hex")).hexdigest()
            self._library_id = self._fetch_library_id(unique_id, True)
            f.write(unique_id)

        self._db.commit()

    def _open_database(self):
        """
        Open a connection to the SQL database.
        
        :return:
            Return a tuple containing a database object and a cursor.
        """

        raise NotImplementedError("The _open_database method must be implemented separately for each SQL "
                                  "implementation")

    def _create_database(self):
        """
        Create a new spectrum library in the SQL database. In some implementations, this may create a whole SQL
        database. In other implementations, multiple spectrum libraries may share a single SQL database.
        
        :return:
            None
        """

        raise NotImplementedError("The _create_database method must be implemented separately for each SQL "
                                  "implementation")

    def close(self):
        self._db.commit()
        self._db.close()

    def refresh_database(self):
        self._db.commit()
        self._db.close()
        self._open_database()

    def __str__(self):
        return "<{module}.{name} instance with path <{path}>".format(module=self.__module__,
                                                                     name=type(self).__name__,
                                                                     path=self._path)

    def _metadata_init(self):
        """
        Create an internal look-up list of the metadata fields defined on spectra in this library.
        
        :return:
            None
        """

        # Create list of available metadata fields
        self._metadata_fields = []
        self._parameterised_query("SELECT fieldId, name FROM metadata_fields;")
        for item in self._db_cursor:
            self._metadata_fields.append(item["name"])

    def _fetch_library_id(self, name, add_record=False):
        """
        Look up the database internal id used to represent a particular spectrum library.
        
        :param name:
            A named spectrum library.
            
        :type name:
            str
            
        :param add_record:
            If true, this is a new library that we are expecting to add to the database. If false, we expect this
            library to already exist.
            
        :type add_record:
            bool
            
        :return:
            Integer origin id.
        """

        while True:
            # Look up whether this origin name already exists in the database
            self._parameterised_query("SELECT libraryId FROM libraries WHERE name=?;", (name,))
            results = self._db_cursor.fetchall()

            assert not (len(results) > 0 and add_record), \
                "Attempting to create a library which already exists in the database."

            assert not (len(results) == 0 and (not add_record)), \
                "Attempting to access a library which doesn't exist in the database."

            if results:
                return results[0]["libraryId"]

            # If not, add it into the database
            self._parameterised_query("INSERT INTO libraries (name) VALUES (?);", (name,))
            self._db.commit()
            add_record = False

    def _fetch_origin_id(self, name):
        """
        Look up the database internal id used to represent a particular named data origin. All spectra are required
        to have a named origin, usually the name of the module used to generate it.
        
        :param name:
            A named data origin.
            
        :type name:
            str
            
        :return:
            Integer origin id.
        """

        while True:
            # Look up whether this origin name already exists in the database
            self._parameterised_query("SELECT originId FROM origins WHERE name=?;", (name,))
            results = self._db_cursor.fetchall()
            if results:
                return results[0]["originId"]

            # If not, add it into the database
            self._parameterised_query("INSERT INTO origins (name) VALUES (?);", (name,))
            self._db.commit()

    def _fetch_metadata_field_id(self, name):
        """
        Look up the database internal id used to represent a particular named metadata field. Spectra can have arbitrary
        metadata associated with them.
        
        :param name:
            A named metadata field.
            
        :type name:
            str
            
        :return:
            Integer field id.
        """

        while True:
            # Look up whether this metadata field already exists in the database
            self._parameterised_query("SELECT fieldId FROM metadata_fields WHERE name=?;", (name,))
            results = self._db_cursor.fetchall()
            if results:
                return results[0]["fieldId"]

            # If not, add it into the database
            self._parameterised_query("INSERT INTO metadata_fields (name) VALUES (?);", (name,))
            self._db.commit()
            self._metadata_init()

    def _filenames_to_ids(self, filenames):
        """
        Convert a list of spectra filenames into database Ids. This helper function is used by various methods which
        can act on spectra referenced either by database id or by filename.
         
        :param filenames:
            A list of filenames whose database ids should be looked up.
            
        :type filenames:
            list of str
            
        :return:
            list of int
        """

        if len(filenames) == 0:
            return []

        output = []
        filename_strings = ["\"%s\"" % re.sub("[^a-zA-Z0-9_-]", "", i) for i in filenames]
        self._parameterised_query("""
SELECT specId FROM spectra WHERE libraryId=%s AND filename IN (%s);
""" % (self._library_id, ",".join(filename_strings)))

        for item in self._db_cursor:
            output.append(item["specId"])
        return output

    def _ids_to_filenames(self, ids):
        """
        Convert a list of database Ids into spectra filenames. This helper function is used by various methods which
        can act on spectra referenced either by database id or by filename.
         
        :param ids:
            A list of database ids whose filenames should be looked up.
            
        :type ids:
            list of int
            
        :return:
            list of str
        """

        if len(ids) == 0:
            return []

        output = []
        id_strings = [str(int(i)) for i in ids]
        self._parameterised_query("""
SELECT filename FROM spectra WHERE libraryId=%s AND specId IN (%s);
""" % (self._library_id, ",".join(id_strings)))

        for item in self._db_cursor:
            output.append(item["filename"])
        return output

    def purge(self):
        """
        This irrevocably deletes the spectrum library from the database and from your disk. You have been warned.
         
        :return:
            None
        """

        # Delete spectra
        self._parameterised_query("SELECT filename FROM spectra WHERE libraryId=?;", (self._library_id,))
        for item in self._db_cursor:
            os.unlink(os_path.join(self._path, item["filename"]))

        # Delete id files
        os.unlink(os_path.join(self._path, "unique_id"))
        os.unlink(os_path.join(self._path, "type_id"))

        # Delete database entries
        self._parameterised_query("DELETE FROM libraries WHERE libraryId=?;", (self._library_id,))
        self._db.commit()

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

        # Start building a list of SQL search criteria as string fragments
        criteria = ["libraryId = ?"]

        # List of parameters to substitute into the SQL query we are building
        criteria_params = [self._library_id]

        # Loop over metadata constraints
        for key, search_range in kwargs.iteritems():

            # Check that requested metadata field exists
            assert key in self._metadata_fields, "Unknown metadata field <{}>.".format(key)

            constraint = """s.specId IN (SELECT specId FROM spectrum_metadata 
                                         WHERE libraryId=? AND
                                               fieldId=(SELECT fieldId FROM metadata_fields WHERE name=?) AND
                                               ({}) )"""
            criteria_params.append(self._library_id)
            criteria_params.append(key)

            # If constraint is specified as a list, it should be of the form [min, max]
            if isinstance(search_range, (list, tuple)):
                assert len(search_range) == 2, \
                    "Search ranges must have two items, a minimum and a maximum. Supplied range has %d items." % \
                    (format(len(search_range)))

                if isinstance(search_range[0], (int, float)):
                    criteria.append(constraint.format("valueFloat BETWEEN ? AND ?"))
                else:
                    criteria.append(constraint.format("valueString BETWEEN ? AND ?"))
                criteria_params.append(min(search_range))
                criteria_params.append(max(search_range))

            # If constraint is not a list or tuple, we must match its exact value
            else:
                if isinstance(search_range, (int, float)):
                    criteria.append(constraint.format("valueFloat = ?"))
                else:
                    criteria.append(constraint.format("valueString = ?"))
                criteria_params.append(search_range)

        # Assemble our list of search criteria into an SQL query
        query = """
SELECT s.specId, s.filename, o.name AS origin
FROM spectra s
INNER JOIN origins o ON s.originId = o.originId
WHERE {};""".format(" AND ".join(criteria))
        self._parameterised_query(query, criteria_params)
        return self._db_cursor.fetchall()

    @requires_ids_or_filenames
    def get_metadata(self, ids=None, filenames=None):
        """
        Fetch dictionaries of the metadata set on a list of spectra in this library. The list of spectra can be
        specified either as a list of ids or a list of filenames.
        
        :param ids:
            A list of integer ids of the spectra to be queried. Set to None to search by filename instead.
            
        :type ids:
            List of int, or None
            
        :param filenames:
            A list of the filenames of the spectra to be queried. Set to None to search by integer id instead.
            
        :type filenames:
            List of str, or None
            
        :return:
            List of dictionaries containing metadata on the requested spectra
        """

        # If we are searching by filename, turn the list of filenames into a list of ids
        if filenames is not None:
            ids = self._filenames_to_ids(filenames=filenames)

        # Start building a list of output
        output = []

        # Loop over the spectra we are querying
        for id in ids:

            # Start building a dictionary of metadata
            item = {}

            # Search the database for metadata
            self._parameterised_query("""
SELECT f.name, CASE WHEN i.valueFloat IS NOT NULL THEN i.valueFloat ELSE i.valueString END AS value
FROM spectrum_metadata i
INNER JOIN metadata_fields f ON f.fieldId=i.fieldId
WHERE i.libraryId=? AND i.specId=?;
            """, (self._library_id, id))

            # Enter metadata into dictionary
            for entry in self._db_cursor:
                item[entry["name"]] = entry["value"]

            output.append(item)

        # Return list of output
        return output

    @requires_ids_or_filenames
    def set_metadata(self, metadata, ids=None, filenames=None):
        """
        Set metadata fields on a list of spectra within this spectrum library.
        
        :param metadata:
            Dictionary of the metadata fields to be set.
            
        :type metadata:
            dict
            
        :param ids:
            List of the integer ids of the spectra to receive this metadata, or None to select them by filename.
            
        :type ids:
            List of int, or None
            
        :param filenames:
            List of the filenames of the spectra to receive this metadata, or None to select them by integer id.
            
        :type filenames:
            List of str, or None
            
        :return:
            None
        """

        # If we are searching by filename, turn the list of filenames into a list of ids
        if filenames is not None:
            ids = self._filenames_to_ids(filenames=filenames)

        # Loop over the metadata fields we are to set
        for key, value in metadata.iteritems():

            # Look up the numeric id for this metadata field
            keyId = self._fetch_metadata_field_id(name=key)

            # Create a big table of parameters to substitute into an SQL query, acting on every spectrum at once
            query_data = zip([self._library_id] * len(ids), ids, [keyId] * len(ids), [value] * len(ids))

            # If this metadata item has a numeric value, we store it in the SQL field <valueFloat>
            if isinstance(value, (int, float)):
                self._parameterised_query_many("""
REPLACE INTO spectrum_metadata (libraryId, specId, fieldId, valueFloat) VALUES 
(?, ?, ?, ?)""", query_data)

            # ... otherwise we store it in the SQL field <valueString>
            else:
                self._parameterised_query_many("""
REPLACE INTO spectrum_metadata (libraryId, specId, fieldId, valueString) VALUES 
(?, ?, ?, ?)""", query_data)

        # Commit changes into database
        self._db.commit()

    @requires_ids_or_filenames
    def open(self, ids=None, filenames=None, shared_memory=False):
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
            
        :param shared_memory:
            Boolean flag indicating whether this SpectrumArray should use multiprocessing shared memory.
            
        :type shared_memory:
            bool
            
        :return:
            A SpectrumArray object.
        """

        if ids is not None:
            filenames = self._ids_to_filenames(ids=ids)

        metadata_list = []
        for filename in filenames:
            metadata_list.append(self.get_metadata(filenames=(filename,)))
        return SpectrumArray.from_files(path=self._path,
                                        filenames=filenames,
                                        metadata_list=metadata_list,
                                        shared_memory=shared_memory)

    def insert(self, spectra, filenames, origin="Undefined", metadata_list=None, overwrite=False):
        """
        Insert the spectra from a SpectrumArray object into this spectrum library.
        
        :param spectra: 
            A SpectrumArray or single Spectrum object containing the spectra to be inserted into this spectrum library.
            
        :type spectra:
            SpectrumArray or Spectrum
            
        :param filenames:
            A list of the filenames with which to save the spectra contained within this SpectrumArray, or a single
            string if only one spectrum is being inserted.
            
        :type filenames:
            List[str] or str
            
        :param origin: 
            A string describing where these spectra are being imported from. Normally the name of the module which is
            importing them.
            
        :type origin:
            str
            
        :param metadata_list: 
            A list of dictionaries of metadata to set on each of the spectra in this SpectrumArray, or a single
            dictionary to set the same metadata on all spectra, or None to set no metadata.
            
        :type metadata_list:
            List[dict] or Dict or None
            
        :param overwrite:
            Boolean flag indicating whether we're allowed to overwrite pre-existing spectra with the same filenames
            
        :type overwrite:
            bool
            
        :return:
            None
        """

        # Sanity check input
        if not isinstance(filenames, (list, tuple)):
            filenames = [filenames]
        if not isinstance(metadata_list, (list, tuple)):
            metadata_list = [metadata_list] * len(filenames)

        assert len(filenames) == len(metadata_list), "Inconsistent number of items being inserted."

        if isinstance(spectra, Spectrum):
            assert len(filenames) == 1
        elif isinstance(spectra, SpectrumArray):
            assert len(spectra) == len(filenames), "Inconsistent number of items being inserted."
        else:
            raise TypeError("Argument 'spectra' must be either a Spectrum or a SpectrumArray.")

        # Fetch the numerical id of the origin of these spectra
        origin_id = self._fetch_origin_id(origin)

        # Insert each spectrum in turn
        for index, (filename, metadata) in enumerate(zip(filenames, metadata_list)):

            # Write spectrum to text file
            spectrum = spectra if isinstance(spectra, Spectrum) else spectra.extract_item(index)
            spectrum.to_file(filename=os_path.join(self._path, filename), overwrite=overwrite)

            # Create database entry a spectrum
            self._parameterised_query("""
REPLACE INTO spectra (filename, originId, libraryId, importTime)
 VALUES (?, ?, ?, ?);
            """, (filename, origin_id, self._library_id, time.time()))

            # Set metadata on this spectrum
            self.set_metadata(filenames=(filename,), metadata=spectrum.metadata)
            if metadata is not None:
                self.set_metadata(filenames=(filename,), metadata=metadata)

        # Commit changes into database
        self._db.commit()
