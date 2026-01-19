#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Integration Tests for `neteval` package."""

import os

import unittest
#from neteval import netevalcmd

SKIP_REASON = 'NETEVAL_INTEGRATION_TEST ' \
              'environment variable not set, cannot run integration ' \
              'tests'

@unittest.skipUnless(os.getenv('NETEVAL_INTEGRATION_TEST') is not None, SKIP_REASON)
class TestIntegrationNeteval(unittest.TestCase):
    """Tests for `neteval` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_something(self):
        """Tests parse arguments"""
        self.assertEqual(1, 1)
