#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
Unit tests for the MySQL implementation of spectrum libraries
"""

from os import path as os_path
import uuid
import unittest
import fourgp_speclib

from test_spectrum_library_sql import TestSpectrumLibrarySQL

# These tests require a test MySQL database to be present
db_host = "localhost"
db_user = "fourgp_unittest"
db_passwd = "fourgp_unittest"
db_name = "fourgp_unittest"


class TestSpectrumLibraryMySqlCreation(unittest.TestCase):
    def test_database_creation(self):
        """
        Test that we can create a new SpectrumLibrary based on an MySQL database.
        """
        unique_filename = uuid.uuid4()
        db_path = os_path.join("/tmp", "speclib_test_{}".format(unique_filename))
        lib = fourgp_speclib.SpectrumLibraryMySql(path=db_path, create=True, purge_db=True,
                                                  db_user=db_user, db_passwd=db_passwd,
                                                  db_name=db_name, db_host=db_host)
        lib.purge()

    def test_multiple_libraries(self):
        """
        Test that we can create multiple SpectrumLibraries at once.
        """
        unique_filename_1 = uuid.uuid4()
        db_path_1 = os_path.join("/tmp", "speclib_test_{}".format(unique_filename_1))
        unique_filename_2 = uuid.uuid4()
        db_path_2 = os_path.join("/tmp", "speclib_test_{}".format(unique_filename_2))
        lib_1 = fourgp_speclib.SpectrumLibraryMySql(path=db_path_1, create=True, purge_db=True,
                                                    db_user=db_user, db_passwd=db_passwd,
                                                    db_name=db_name, db_host=db_host)
        lib_2 = fourgp_speclib.SpectrumLibraryMySql(path=db_path_2, create=True,
                                                    db_user=db_user, db_passwd=db_passwd,
                                                    db_name=db_name, db_host=db_host)
        lib_3 = fourgp_speclib.SpectrumLibraryMySql(path=db_path_1,
                                                    db_user=db_user, db_passwd=db_passwd,
                                                    db_name=db_name, db_host=db_host)
        lib_4 = fourgp_speclib.SpectrumLibraryMySql(path=db_path_2,
                                                    db_user=db_user, db_passwd=db_passwd,
                                                    db_name=db_name, db_host=db_host)
        lib_3.close()
        lib_4.close()
        lib_1.purge()
        lib_2.purge()

    def test_non_existent_database(self):
        """
        Test that we get an exception if we try to open a SpectrumLibrary that doesn't exist.
        """
        unique_filename = uuid.uuid4()
        db_path = os_path.join("/tmp", "speclib_test_{}".format(unique_filename))
        with self.assertRaises(AssertionError):
            fourgp_speclib.SpectrumLibraryMySql(path=db_path, create=False, purge_db=True,
                                                db_user=db_user, db_passwd=db_passwd,
                                                db_name=db_name, db_host=db_host)


class TestSpectrumLibraryMySQLBinary(unittest.TestCase, TestSpectrumLibrarySQL):
    def setUp(self):
        """
        Open connection to a clean SpectrumLibrary based on MySQL.
        """
        unique_filename = uuid.uuid4()
        self._db_path = os_path.join("/tmp", "speclib_test_{}".format(unique_filename))
        self._lib = fourgp_speclib.SpectrumLibraryMySql(path=self._db_path, create=True, purge_db=True,
                                                        binary_spectra=True,
                                                        db_user=db_user, db_passwd=db_passwd,
                                                        db_name=db_name, db_host=db_host)

    def tearDown(self):
        """
        Tear down SpectrumLibrary based on MySQL.
        """
        self._lib.purge()


class TestSpectrumLibraryMySQLGzip(unittest.TestCase, TestSpectrumLibrarySQL):
    def setUp(self):
        """
        Open connection to a clean SpectrumLibrary based on MySQL.
        """
        unique_filename = uuid.uuid4()
        self._db_path = os_path.join("/tmp", "speclib_test_{}".format(unique_filename))
        self._lib = fourgp_speclib.SpectrumLibraryMySql(path=self._db_path, create=True, purge_db=True,
                                                        gzip_spectra=True, binary_spectra=False,
                                                        db_user=db_user, db_passwd=db_passwd,
                                                        db_name=db_name, db_host=db_host)

    def tearDown(self):
        """
        Tear down SpectrumLibrary based on MySQL.
        """
        self._lib.purge()


class TestSpectrumLibraryMySQLText(unittest.TestCase, TestSpectrumLibrarySQL):
    def setUp(self):
        """
        Open connection to a clean SpectrumLibrary based on MySQL.
        """
        unique_filename = uuid.uuid4()
        self._db_path = os_path.join("/tmp", "speclib_test_{}".format(unique_filename))
        self._lib = fourgp_speclib.SpectrumLibraryMySql(path=self._db_path, create=True, purge_db=True,
                                                        binary_spectra=False, gzip_spectra=False,
                                                        db_user=db_user, db_passwd=db_passwd,
                                                        db_name=db_name, db_host=db_host)

    def tearDown(self):
        """
        Tear down SpectrumLibrary based on MySQL.
        """
        self._lib.purge()


# Run tests if we are run from command line
if __name__ == '__main__':
    unittest.main()
