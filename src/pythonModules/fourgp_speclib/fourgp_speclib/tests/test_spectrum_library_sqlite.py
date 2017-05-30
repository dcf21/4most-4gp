#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for the SQLite implementation of spectrum libraries
"""

from os import path as os_path
import uuid
import unittest
import fourgp_speclib

from test_spectrum_library_sql import TestSpectrumLibrarySQL


class TestSpectrumLibrarySQLiteCreation(unittest.TestCase):
    def test_database_creation(self):
        """
        Test that we can create a new SpectrumLibrary based on an SQLite database.
        """
        unique_filename = uuid.uuid4()
        db_path = os_path.join("/tmp", "speclib_test_{}".format(unique_filename))
        lib = fourgp_speclib.SpectrumLibrarySqlite(path=db_path, create=True)
        lib.purge()

    def test_multiple_libraries(self):
        """
        Test that we can create multiple SpectrumLibraries at once.
        """
        unique_filename_1 = str(uuid.uuid4())
        db_path_1 = os_path.join("/tmp", "speclib_test_{}".format(unique_filename_1))
        unique_filename_2 = str(uuid.uuid4())
        db_path_2 = os_path.join("/tmp", "speclib_test_{}".format(unique_filename_2))
        lib_1 = fourgp_speclib.SpectrumLibrarySqlite(path=db_path_1, create=True)
        lib_2 = fourgp_speclib.SpectrumLibrarySqlite(path=db_path_2, create=True)
        lib_3 = fourgp_speclib.SpectrumLibrarySqlite(path=db_path_1)
        lib_4 = fourgp_speclib.SpectrumLibrarySqlite(path=db_path_2)
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
            fourgp_speclib.SpectrumLibrarySqlite(path=db_path, create=False)


class TestSpectrumLibrarySQLiteBinary(unittest.TestCase, TestSpectrumLibrarySQL):
    def setUp(self):
        """
        Open connection to a clean SpectrumLibrary based on SQLite.
        """
        unique_filename = uuid.uuid4()
        self._db_path = os_path.join("/tmp", "speclib_test_{}".format(unique_filename))
        self._lib = fourgp_speclib.SpectrumLibrarySqlite(path=self._db_path, create=True,
                                                         binary_spectra=True)

    def tearDown(self):
        """
        Tear down SpectrumLibrary based on SQLite.
        """
        self._lib.purge()


class TestSpectrumLibrarySQLiteGzip(unittest.TestCase, TestSpectrumLibrarySQL):
    def setUp(self):
        """
        Open connection to a clean SpectrumLibrary based on SQLite.
        """
        unique_filename = uuid.uuid4()
        self._db_path = os_path.join("/tmp", "speclib_test_{}".format(unique_filename))
        self._lib = fourgp_speclib.SpectrumLibrarySqlite(path=self._db_path, create=True,
                                                         binary_spectra=False, gzip_spectra=True)

    def tearDown(self):
        """
        Tear down SpectrumLibrary based on SQLite.
        """
        self._lib.purge()


class TestSpectrumLibrarySQLiteText(unittest.TestCase, TestSpectrumLibrarySQL):
    def setUp(self):
        """
        Open connection to a clean SpectrumLibrary based on SQLite.
        """
        unique_filename = uuid.uuid4()
        self._db_path = os_path.join("/tmp", "speclib_test_{}".format(unique_filename))
        self._lib = fourgp_speclib.SpectrumLibrarySqlite(path=self._db_path, create=True,
                                                         binary_spectra=False, gzip_spectra=False)

    def tearDown(self):
        """
        Tear down SpectrumLibrary based on SQLite.
        """
        self._lib.purge()


# Run tests if we are run from command line
if __name__ == '__main__':
    unittest.main()