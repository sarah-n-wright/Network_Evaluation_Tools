#! /usr/bin/env python

import logging


logger = logging.getLogger(__name__)


class NetevalRunner(object):
    """
    Class to run algorithm
    """
    def __init__(self, exitcode):
        """
        Constructor

        :param exitcode: value to return via :py:meth:`.NetevalRunner.run` method
        :type int:
        """
        self._exitcode = exitcode
        logger.debug('In constructor')

    def run(self):
        """
        Runs Network Evaluations Tools v2


        :return:
        """
        logger.debug('In run method')
        return self._exitcode
