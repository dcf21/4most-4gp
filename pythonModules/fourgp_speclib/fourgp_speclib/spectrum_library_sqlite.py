#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from os import path as os_path
import sqlite3

from spectrum_library import SpectrumLibrary, requires_ids_or_filenames


class SpectrumLibrarySqlite(SpectrumLibrary):
    """
    A spectrum library implementation that uses SQLite3 to store metadata about each spectrum.
    """

    _index_file_name = "index.db"

    _schema = """

-- Table of spectrum libraries using this database
CREATE TABLE libraries (
    libraryId INTEGER PRIMARY KEY,
    name VARCHAR(1024) UNIQUE NOT NULL
);

-- Table of string descriptions of which tools imported particular spectra into the library 
CREATE TABLE origins (
    originId INTEGER PRIMARY KEY,
    name VARCHAR(256) UNIQUE NOT NULL
);

-- Table of spectra within this library
CREATE TABLE spectra (
    specId INTEGER PRIMARY KEY,
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
    fieldId INTEGER PRIMARY KEY,
    name VARCHAR(32) UNIQUE NOT NULL
);

-- Table of all metadata items set of spectra
CREATE TABLE spectrum_metadata (
    specId INTEGER NOT NULL,
    fieldId INTEGER NOT NULL,
    libraryId INTEGER NOT NULL,  -- this field is also in the record in <spectra>, but copied here for fast searching
    valueFloat REAL,
    valueString VARCHAR(256),
    FOREIGN KEY (specId) REFERENCES spectra(specId) ON DELETE CASCADE,
    FOREIGN KEY (fieldId) REFERENCES metadata_fields(fieldId) ON DELETE CASCADE,
    FOREIGN KEY (libraryId) REFERENCES libraries (libraryId) ON DELETE CASCADE
);

CREATE INDEX search_metadata_floats ON spectrum_metadata (libraryId, fieldId, valueFloat);
CREATE INDEX search_metadata_strings ON spectrum_metadata (libraryId, fieldId, valueString);
    
    """

    def __init__(self, path, create=False):
        """
        Create a new SpectrumLibrary object, storing metadata about the spectra in an SQLite database.
        
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

        self._library_id = "0"

        # Create new spectrum library if requested
        if create:
            self._create(path)

        # Check that we're not overwriting an existing library
        assert os_path.exists(path), \
            "Could not open spectrum library <{}>: directory not found".format(path)

        self._path = path
        self._path_db = os_path.join(path, self._index_file_name)
        self._db = sqlite3.connect(self._path_db)
        self._db_cursor = self._db.cursor()

        # Initialise
        self._metadata_init()

        super(SpectrumLibrarySqlite, self).__init__()

    def _create(self, path):
        """
        Create a new, empty spectrum library.
        
        :param path:
            The file path to use for storing spectra in this library.
            
        :type path:
            str
            
        :return:
            None
        """

        # Check that we're not overwriting an existing library
        assert not os_path.exists(path), \
            "Could not create spectrum library <{}>: file already exists".format(path)

        # Check that parent directory exists
        parent_path, file_name = os_path.split(path)
        assert os_path.exists(parent_path), \
            "Could not create spectrum library <{}>: parent directory does not exist".format(path)

        # Create directory to hold spectra in this library
        try:
            os.mkdir(path)
        except IOError:
            raise

        # Create SQLite database to hold metadata about these spectra
        db_path = os_path.join(path, self._index_file_name)
        db = sqlite3.connect(db_path)
        c = db.cursor()
        c.executescript(self._schema)
        c.execute("INSERT INTO libraries (name) VALUES (%s)", (self._library_id,))
        db.commit()
        db.close()

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
        self._db_cursor.execute("SELECT fieldId, name FROM metadata_fields;")
        for item in self._db_cursor:
            self._metadata_fields.append(item['name'])

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
            self._db_cursor.execute("SELECT originId FROM origins WHERE name=%s;", (name,))
            results = self._db_cursor.fetchall()
            if results:
                return results[0]['originId']

            # If not, add it into the database
            self._db_cursor.execute("INSERT INTO origins (name) VALUES (%s);", (name,))

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
            self._db_cursor.execute("SELECT fieldId FROM metadata_fields WHERE name=%s;", (name,))
            results = self._db_cursor.fetchall()
            if results:
                return results[0]['fieldId']

            # If not, add it into the database
            self._db_cursor.execute("INSERT INTO metadata_fields (name) VALUES (%s);", (name,))
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

        output = []
        self._db_cursor.execute("""
SELECT specId FROM spectra WHERE libraryId=%s AND filename IN %s;
""", (self._library_id, filenames))

        for item in self._db_cursor:
            output.append(item["specId"])
        return output

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
        criteria = ["libraryId = %s"]

        # List of parameters to substitute into the SQL query we are building
        criteria_params = [self._library_id]

        # Loop over metadata constraints
        for key, search_range in kwargs.iteritems():

            # Check that requested metadata field exists
            assert key in self._metadata_fields, "Unknown metadata field <{}>.".format(key)

            # If constraint is specified as a list, it should be of the form [min, max]
            if isinstance(search_range, (list, tuple)):
                assert len(search_range) == 2, \
                    "Search ranges must have two items, a minimum and a maximum. Supplied range has %d items." % \
                    (format(len(search_range)))
                criteria.append("""
EXISTS (SELECT 1
    FROM spectrum_metadata i
    INNER JOIN metadata_fields f ON f.fieldId = i.fieldId
    WHERE f.libraryId={} AND f.name={} AND ((i.valueFloat BETWEEN %s AND %s) OR (i.valueString BETWEEN %s AND %s)) )
                """.format(self._library_id, key))
                criteria_params.append(min(search_range))
                criteria_params.append(max(search_range))
                criteria_params.append(min(search_range))
                criteria_params.append(max(search_range))

            # If constraint is not a list or tuple, we must match its exact value
            else:
                criteria.append("""
EXISTS (SELECT 1
    FROM spectrum_metadata i
    INNER JOIN metadata_fields f ON f.fieldId = i.fieldId
    WHERE f.libraryId={} AND f.name={} AND ((i.valueFloat = %s) OR (i.valueString = %s)) )
                """.format(self._library_id, key))
                criteria_params.append(search_range)
                criteria_params.append(search_range)

        # Assemble our list of search criteria into an SQL query
        query = "SELECT specId, filename FROM spectra WHERE {};".format(" AND ".join(criteria))
        self._db_cursor.execute(query, criteria_params)
        return self._db_cursor.fetchall()

    def list_metadata_fields(self):
        return self._metadata_fields

    @requires_ids_or_filenames
    def get_metadata(self, ids=None, filenames=None):
        if filenames is not None:
            ids = self._filenames_to_ids(filenames=filenames)
        output = []
        for id in ids:
            item = {}
            self._db_cursor.execute("""
SELECT f.name, (CASE WHEN i.valueFloat IS NOT NULL THEN i.valueFloat ELSE i.valueString) AS value
FROM spectrum_metadata i
INNER JOIN metadata_fields f ON f.fieldId=i.fieldId
WHERE i.libraryId=%s AND i.specId=%s;
            """, (self._library_id, id))
            for entry in self._db_cursor:
                item[entry["name"]] = entry["value"]
            output.append(item)
        return output

    @requires_ids_or_filenames
    def set_metadata(self, metadata, ids=None, filenames=None):
        if filenames is not None:
            ids = self._filenames_to_ids(filenames=filenames)
        for key, value in metadata.iteritems():
            keyId = self._fetch_metadata_field_id(name=key)
            query_data = zip([self._library_id] * len(ids), [keyId] * len(ids), [value] * len(ids), ids)
            if isinstance(value, (int, float)):
                self._db_cursor.executemany("""
REPLACE INTO spectrum_metadata (libraryId, specId, fieldId, valueFloat) VALUES 
(%s, %s, (SELECT fieldId FROM metadata_fields WHERE name=%s), %s)""", query_data)

            else:
                self._db_cursor.executemany("""
REPLACE INTO spectrum_metadata (libraryId, specId, fieldId, valueString) VALUES 
(%s, %s, (SELECT fieldId FROM metadata_fields WHERE name=%s), %s)""", query_data)

    @requires_ids_or_filenames
    def open(self, ids=None, filenames=None):
        if filenames is not None:
            ids = self._filenames_to_ids(filenames=filenames)
        metadata = self.get_metadata(filenames=(filename,))

    def insert(self, spectrum, filename, origin="Undefined", metadata=None, overwrite=False):
        origin_id = self._fetch_origin_id(origin)

        self._db_cursor.execute("""
REPLACE INTO spectra (filename, originId, importTime)
 VALUES (%s, %s, (JULIANDAY('now') - 2440587.5) * 86400.0);
        """, (filename, origin_id))

        if metadata is not None:
            self.set_metadata(filenames=(filename,), metadata=metadata)

    def refresh_database(self):
        self._db.commit()
        self._db.close()
        self._db = sqlite3.connect(self._path_db)
