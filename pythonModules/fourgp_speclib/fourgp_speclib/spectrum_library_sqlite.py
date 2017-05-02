#!/usr/bin/env python
# -*- coding: utf-8 -*-

from os import path as os_path
import sqlite3

from spectrum_library_sql import SpectrumLibrarySql


class SpectrumLibrarySqlite(SpectrumLibrarySql):
    """
    A spectrum library implementation that uses SQLite3 to store metadata about each spectrum.
    """

    _index_file_name = "index.db"

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
        super(SpectrumLibrarySqlite, self).__init__(path=path, create=create)

    def _create_database(self):
        """
        Create a new, empty spectrum library.
        
        :param path:
            The file path to use for storing spectra in this library.
            
        :type path:
            str
            
        :return:
            None
        """

        # Create SQLite database to hold metadata about the spectra in this library
        db_path = os_path.join(self._path, self._index_file_name)
        db = sqlite3.connect(db_path)
        c = db.cursor()
        c.executescript(self._schema)
        c.execute("INSERT INTO libraries (name) VALUES (%s)", (self._library_id,))
        db.commit()
        db.close()

    def _open_database(self):
        self._path_db = os_path.join(self._path, self._index_file_name)
        self._db = sqlite3.connect(self._path_db)
        self._db_cursor = self._db.cursor()
        return (self._db, self._db_cursor)

    def refresh_database(self):
        self._db.commit()
        self._db.close()
        self._db = sqlite3.connect(self._path_db)