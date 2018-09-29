# -*- coding: utf-8 -*-

import re
import MySQLdb

from .spectrum_library_sql import SpectrumLibrarySql


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

    def __init__(self, path, create=False, gzip_spectra=True, binary_spectra=True, purge_db=False,
                 db_user="fourgp", db_passwd="fourgp", db_name="fourgp", db_host="localhost"):
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
         
        :param purge_db:
            If true, wipe the database clean and start a new schema. Warning: This will trash everything in the
            database, including any other SpectrumLibraries which share the same MySQL database.
        
        :type purge_db:
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
        self._purge_db = purge_db

        self._db = None
        self._db_cursor = None

        super(SpectrumLibraryMySql, self).__init__(path=path, create=create,
                                                   gzip_spectra=gzip_spectra, binary_spectra=binary_spectra)

    def _create_database(self):
        """
        Create a database record for a new, empty spectrum library.
            
        :return:
            None
        """

        # If we've been asked the purge the database, do that now
        if self._purge_db:
            db = MySQLdb.connect(host=self._db_host, user=self._db_user, passwd=self._db_passwd, db=self._db_name)
            c = db.cursor(cursorclass=MySQLdb.cursors.DictCursor)
            c.execute("DROP DATABASE IF EXISTS {};".format(self._db_name))
            c.execute("CREATE DATABASE {};".format(self._db_name))
            db.commit()
            db.close()

        # Create MySQL database to hold metadata about the spectra in this library
        db = MySQLdb.connect(host=self._db_host, user=self._db_user, passwd=self._db_passwd, db=self._db_name)
        c = db.cursor(cursorclass=MySQLdb.cursors.DictCursor)

        # Test if database tables have already been set up. If not, build database from scratch
        try:
            c.execute("SELECT 1 FROM libraries;")
        except MySQLdb.ProgrammingError:
            for line in self._schema.split(";"):
                line = line.strip()
                if line:
                    c.execute(line)
        db.commit()
        db.close()

    def _open_database(self):
        self._db = MySQLdb.connect(host=self._db_host, user=self._db_user, passwd=self._db_passwd, db=self._db_name)
        self._db_cursor = self._db.cursor(cursorclass=MySQLdb.cursors.Cursor)
        return self._db, self._db_cursor

    def _parameterised_query(self, sql, parameters=None):
        if parameters is None:
            parameters = ()
        sql = re.sub(r"\?", r"%s", sql)
        self._db_cursor.execute(sql, parameters)

    def _parameterised_query_many(self, sql, parameters=None):
        sql = re.sub(r"\?", r"%s", sql)
        self._db_cursor.executemany(sql, parameters)
