#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from os import path as os_path
import sqlite3

from spectrum_library_sql import SpectrumLibrarySql


class SpectrumLibrarySqlite(SpectrumLibrarySql):
    """
    A spectrum library implementation that uses SQLite3 to store metadata about each spectrum.
    
    :cvar string _index_file_name:
        The filename used to store the SQLite database within the directory holding this SpectrumLibrary.
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

        self._db = None
        self._db_cursor = None

        super(SpectrumLibrarySqlite, self).__init__(path=path, create=create)

    def _create_database(self):
        """
        Create a database record for a new, empty spectrum library.
            
        :return:
            None
        """

        # Create SQLite database to hold metadata about the spectra in this library
        db_path = os_path.join(self._path, self._index_file_name)

        assert not os_path.exists(db_path),\
            "Attempting to overwrite SQLite database <{}> that already exists.".format(db_path)

        db = sqlite3.connect(db_path)
        c = db.cursor()
        c.executescript(self._schema)
        c.execute("INSERT INTO libraries (name) VALUES (%s)", (self._library_id,))
        db.commit()
        db.close()

    def _open_database(self):
        self._path_db = os_path.join(self._path, self._index_file_name)

        assert os_path.exists(self._path_db),\
            "Attempting to open an SQLite database <{}> that doesn't exist.".format(db_path)

        self._db = sqlite3.connect(self._path_db)
        self._db_cursor = self._db.cursor()
        return self._db, self._db_cursor

    def purge(self):
        """
        This irrevocably deletes the spectrum library from the database and from your disk. You have been warned.
         
        :return:
            None
        """

        super(SpectrumLibrarySqlite, self).purge()

        # Delete SQLite file
        self._db.close()
        os.unlink(os_path.join(self._path_db))