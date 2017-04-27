#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path
import sqlite3

from spectrum_library import SpectrumLibrary


class SpectrumLibrarySqlite(SpectrumLibrary):
    """
    A spectrum library implementation that uses SQLite3 to store metadata about each spectrum.
    """

    _index_file_name = "index.db"

    _schema = """

-- Table of string descriptions of which tools imported particular spectra into the library 
CREATE TABLE origins (
    originId INTEGER PRIMARY KEY,
    name VARCHAR(256)
);

-- Table of spectra within this library
CREATE TABLE spectra (
    specId INTEGER PRIMARY KEY,
    filename VARCHAR(256),
    originId INTEGER,
    importTime REAL,
    FOREIGN KEY (originId) REFERENCES origins (originId) ON DELETE CASCADE
);

-- Table of metadata fields which have been set on at least one spectrum
CREATE TABLE metadata_fields (
    fieldId INTEGER PRIMARY KEY,
    name VARCHAR(32)
);

-- Table of all metadata items set of spectra
CREATE TABLE spectrum_metadata (
    specId INTEGER,
    fieldId INTEGER,
    valueFloat REAL,
    valueString VARCHAR(256),
    FOREIGN KEY (specId) REFERENCES spectra(specId) ON DELETE CASCADE,
    FOREIGN KEY (fieldId) REFERENCES spectrum_metadata(fieldId) ON DELETE CASCADE
);
    
    """

    def __init__(self, path, create=False, *args, **kwargs):

        # Create new spectrum library if requested
        if create:
            self._create(path)

        # Check that we're not overwriting an existing library
        assert os.path.exists(path), \
            "Could not open spectrum library <{}>: directory not found".format(path)

        self._path = path
        self._path_db = os.path.join(path, self._index_file_name)
        self._db = sqlite3.open(self._path_db)
        self._db_cursor = self._db.cursor()

        super(SpectrumLibrarySqlite, self).__init__(*args, **kwargs)

    def _create(self, path):
        # Check that we're not overwriting an existing library
        assert not os.path.exists(path), \
            "Could not create spectrum library <{}>: file already exists".format(path)

        # Check that parent directory exists
        parent_path, file_name = os.path.split(path)
        assert os.path.exists(parent_path), \
            "Could not create spectrum library <{}>: parent directory does not exist".format(path)

        # Create directory to hold spectra in this library
        try:
            os.mkdir(path)
        except IOError:
            raise

        # Create SQLite database to hold metadata about these spectra
        db_path = os.path.join(path, self._index_file_name)
        db = sqlite3.connect(db_path)
        c = db.cursor()
        c.executescript(self._schema)
        db.commit()
        db.close()

    def refresh_database(self):
        self._db.commit()
        self._db.close()
        self._db = sqlite3.open(self._path_db)
