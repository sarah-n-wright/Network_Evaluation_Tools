Integration testing
=======================

To support `integration testing <https://en.wikipedia.org/wiki/Integration_testing>`__ the unit tests in **neteval**
include a parallel set of tests reside in the existing test framework and
can be activated if ``NETEVAL_INTEGRATION_TEST`` environment
variable is set to any value:

Example variable:

.. code-block::

    export NETEVAL_INTEGRATION_TEST="true"
    make test
