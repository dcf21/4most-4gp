#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from os import path as os_path
import MySQLdb

from spectrum_library_sql import SpectrumLibrarySQL


class SpectrumLibraryMySql(SpectrumLibrarySQL):
    """
    A spectrum library implementation that uses MySQL to store metadata about each spectrum.
    """

    def __init__(self, path, create=False, db_user="fourgp", db_passwd="fourgp", db_name="fourgp"):
        """
        Create a new SpectrumLibrary object, storing metadata about the spectra in a MySQL database.
        
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

        super(SpectrumLibraryMySql, self).__init__()

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

    def refresh_database(self):
        self._db.commit()
        self._db.close()
        self._db = sqlite3.connect(self._path_db)
