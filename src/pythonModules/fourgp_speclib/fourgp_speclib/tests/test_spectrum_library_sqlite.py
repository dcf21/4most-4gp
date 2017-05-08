#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for the SQLite implementation of spectrum libraries
"""

from os import path as os_path
import uuid
import unittest
import numpy as np
import fourgp_speclib


class TestSpectrumLibrarySQLiteCreation(unittest.TestCase):
    def test_database_creation(self):
        unique_filename = uuid.uuid4()
        db_path = os_path.join("/tmp", "speclib_test_{}".format(unique_filename))
        lib = fourgp_speclib.SpectrumLibrarySqlite(path=db_path, create=True)
        lib.purge()

    def test_non_existent_database(self):
        unique_filename = uuid.uuid4()
        db_path = os_path.join("/tmp", "speclib_test_{}".format(unique_filename))
        with self.assertRaises(AssertionError):
            fourgp_speclib.SpectrumLibrarySqlite(path=db_path, create=False)


class TestSpectrumLibrarySQLite(unittest.TestCase):
    def setUp(self):
        unique_filename = uuid.uuid4()
        self._db_path = os_path.join("/tmp", "speclib_test_{}".format(unique_filename))
        self._lib = fourgp_speclib.SpectrumLibrarySqlite(path=self._db_path, create=True)

    def tearDown(self):
        self._lib.purge()


# Run tests if we are run from command line
if __name__ == '__main__':
    unittest.main()