#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `neteval` package."""


import unittest
from neteval.runner import NetevalRunner


class TestNetevalrunner(unittest.TestCase):
    """Tests for `neteval` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_constructor(self):
        """Tests constructor"""
        myobj = NetevalRunner(0)

        self.assertIsNotNone(myobj)

    def test_run(self):
        """ Tests run()"""
        myobj = NetevalRunner(4)
        self.assertEqual(4, myobj.run())
