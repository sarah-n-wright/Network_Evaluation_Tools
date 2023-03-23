============================
Network Evaluations Tools v2
============================


.. image:: https://img.shields.io/pypi/v/neteval.svg
        :target: https://pypi.python.org/pypi/neteval

.. image:: https://img.shields.io/travis/sarah-n-wright/neteval.svg
        :target: https://travis-ci.com/sarah-n-wright/neteval

.. image:: https://readthedocs.org/projects/neteval/badge/?version=latest
        :target: https://neteval.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




Python Boilerplate contains all the boilerplate you need to create a Python package with command line


* Free software: MIT license
* Documentation: https://neteval.readthedocs.io.



Dependencies
------------

* TODO add

Compatibility
-------------

* Python 3.3+

Installation
------------

.. code-block::

   git clone https://github.com/sarah-n-wright/neteval
   cd neteval
   make dist
   pip install dist/netevalcmd*whl


Run **make** command with no arguments to see other build/deploy options including creation of Docker image 

.. code-block::

   make

Output:

.. code-block::

   clean                remove all build, test, coverage and Python artifacts
   clean-build          remove build artifacts
   clean-pyc            remove Python file artifacts
   clean-test           remove test and coverage artifacts
   lint                 check style with flake8
   test                 run tests quickly with the default Python
   test-all             run tests on every Python version with tox
   coverage             check code coverage quickly with the default Python
   docs                 generate Sphinx HTML documentation, including API docs
   servedocs            compile the docs watching for changes
   testrelease          package and upload a TEST release
   release              package and upload a release
   dist                 builds source and wheel package
   install              install the package to the active Python's site-packages
   dockerbuild          build docker image and store in local repository
   dockerpush           push image to dockerhub

For developers
-------------------------------------------

To deploy development versions of this package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below are steps to make changes to this code base, deploy, and then run
against those changes.

#. Make changes

   Modify code in this repo as desired

#. Build and deploy

.. code-block::

    # From base directory of this repo neteval
    pip uninstall neteval -y ; make clean dist; pip install dist/neteval*whl



Needed files
------------

**TODO:** Add description of needed files


Usage
-----

For information invoke :code:`netevalcmd.py -h`

**Example usage**

**TODO:** Add information about example usage

.. code-block::

   netevalcmd.py # TODO Add other needed arguments here


Via Docker
~~~~~~~~~~~~~~~~~~~~~~

**Example usage**

**TODO:** Add information about example usage


.. code-block::

   docker run -v `pwd`:`pwd` -w `pwd` sarah-n-wright/neteval:0.1.0 netevalcmd.py # TODO Add other needed arguments here


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _NDEx: http://www.ndexbio.org
