#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from os import path as os_path
import re
import sqlite3

from spectrum_library_sql import SpectrumLibrarySql


class SpectrumLibrarySqlite(SpectrumLibrarySql):
    """
    A spectrum library implementation that uses SQLite3 to store metadata about each spectrum.
    
    :cvar string _index_file_name:
        The filename used to store the SQLite database within the directory holding this SpectrumLibrary.
    """

    _index_file_name = "index.db"

    def __init__(self, path, create=False, gzip_spectra=True, binary_spectra=True):
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

        :param gzip_spectra:
            If true, we store spectra on disk in gzipped text files. This reduces file size by 90%.
            This setting is a property of stored when new libraries are created, and the argument is ignored if
            we are not creating a new library.

        :type gzip_spectra:
            bool

        :param binary_spectra:
            If true, we store spectra on disk in binary format.
            This setting is a property of stored when new libraries are created, and the argument is ignored if
            we are not creating a new library.

        :type binary_spectra:
            bool
        """

        self._db = None
        self._db_cursor = None

        super(SpectrumLibrarySqlite, self).__init__(path=path, create=create,
                                                    gzip_spectra=gzip_spectra, binary_spectra=binary_spectra)

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

        # SQLite databases work faster if primary keys don't auto increment, so remove keyword from schema
        db = sqlite3.connect(db_path)
        c = db.cursor()
        schema = re.sub("AUTO_INCREMENT", "", self._schema)
        c.executescript(schema)
        db.commit()
        db.close()

    def _open_database(self):
        self._path_db = os_path.join(self._path, self._index_file_name)

        assert os_path.exists(self._path_db),\
            "Attempting to open an SQLite database <{}> that doesn't exist.".format(self._path_db)

        self._db = sqlite3.connect(self._path_db)
        self._db.row_factory = sqlite3.Row
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

    def _parameterised_query(self, sql, parameters=None):
        if parameters is None:
            parameters = ()
        self._db_cursor.execute(sql, parameters)

    def _parameterised_query_many(self, sql, parameters=None):
        self._db_cursor.executemany(sql, parameters)