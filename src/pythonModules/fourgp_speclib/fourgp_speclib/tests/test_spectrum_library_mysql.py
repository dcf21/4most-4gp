#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for the MySQL implementation of spectrum libraries
"""

from os import path as os_path
import uuid
import unittest
import numpy as np
import fourgp_speclib

from .test_spectrum_library_sql import TestSpectrumLibrarySQL

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


class TestSpectrumLibraryMySQL(unittest.TestCase, TestSpectrumLibrarySQL):
    def setUp(self):
        """
        Open connection to a clean SpectrumLibrary based on MySQL.
        """
        unique_filename = uuid.uuid4()
        self._db_path = os_path.join("/tmp", "speclib_test_{}".format(unique_filename))
        self._lib = fourgp_speclib.SpectrumLibraryMySql(path=self._db_path, create=True, purge_db=True,
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