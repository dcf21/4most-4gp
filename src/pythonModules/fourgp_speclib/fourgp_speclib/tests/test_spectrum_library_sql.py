#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for all SQL implementations of spectrum libraries
"""

from os import path as os_path
import uuid
import unittest
import numpy as np
import fourgp_speclib


class TestSpectrumLibrarySQL(object):
    """
    This class is a mixin which adds lots of standard tests to any SQL-based SpectrumLibrary unittest class.
    """

    def test_refresh(self):
        """
        Check that we can refresh database connection.
        """
        self._lib.refresh_database()
