#!/usr/bin/env python
# -*- coding: utf-8 -*-

import MySQLdb

from spectrum_library_sql import SpectrumLibrarySql


class SpectrumLibraryMySql(SpectrumLibrarySql):
    """
    A spectrum library implementation that uses MySQL to store metadata about each spectrum.
    
    :ivar string _db_user:
            Username for connecting to MySQL server
            
    :ivar string _db_passwd:
        Password for connecting into MySQL server
        
    :ivar string _db_name:
        Name of the database to use on the MySQL server
        
    :ivar string _db_host:
        Hostname of the MySQL server
    """

    def __init__(self, path, create=False, db_user="fourgp", db_passwd="fourgp", db_name="fourgp", db_host="localhost"):
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
            
        :param db_user:
            Username for connecting to MySQL server
            
        :type db_user:
            str
            
        :param db_passwd:
            Password for connecting into MySQL server
            
        :type db_passwd:
            str
            
        :param db_name:
            Name of the database to use on the MySQL server
            
        :type db_name:
            str
            
        :param db_host:
            Hostname of the MySQL server
            
        :type db_host:
            str
        """

        self._db_user = db_user
        self._db_passwd = db_passwd
        self._db_name = db_name
        self._db_host = db_host

        self._db = None
        self._db_cursor = None

        super(SpectrumLibraryMySql, self).__init__(path=path, create=create)

    def _create_database(self):
        """
        Create a database record for a new, empty spectrum library.
            
        :return:
            None
        """

        # Create MySQL database to hold metadata about the spectra in this library
        db = MySQLdb.connect(host=self._db_host, user=self._db_user, passwd=self._db_passwd, db=self._db_name)
        c = db.cursor(cursorclass=MySQLdb.cursors.DictCursor)
        c.executescript(self._schema)
        c.execute("INSERT INTO libraries (name) VALUES (%s)", (self._library_id,))
        db.commit()
        db.close()

    def _open_database(self):
        self._db = MySQLdb.connect(host=self._db_host, user=self._db_user, passwd=self._db_passwd, db=self._db_name)
        self._db_cursor = self._db.cursor(cursorclass=MySQLdb.cursors.DictCursor)
        return self._db, self._db_cursor