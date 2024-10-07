============================
Network Evaluations Tools v2
============================


.. image:: https://img.shields.io/pypi/v/neteval.svg
        :target: https://pypi.python.org/pypi/neteval

.. image:: https://readthedocs.org/projects/neteval/badge/?version=latest
        :target: https://neteval.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


Pyhton package for processing and evaluation biological networks. All included analysis tasks are described in 
Wright, SN., et al. *State of The Interactomes: an evaluation of molecular networks for generating biological insights.* 

Free software: MIT license

Installation
------------
To install the package, run the following command:

   .. code-block::

      pip install neteval

To install the package from source, clone the repository and run the following command:

   .. code-block::
   
      git clone https://github.com/sarah-n-wright/Network_Evaluation_Tools
      cd Network_Evaluation_Tools
      make dist
      pip install dist/netevalcmd*whl

For example usage of command line scripts see `Example Usage <https://github.com/sarah-n-wright/Network_Evaluation_Tools/ExampleUsage>`__.  
For example usage of all other funcitonality see `State of the Interactomes Notebooks <https://github.com/sarah-n-wright/Network_Evaluation_Tools/StateOfTheInteractomes_Notebooks>`__.  

Dependencies
------------

* `goatools >=1.3.1 <https://pypi.org/project/goatools>`__
* `hidef >=1.1.5 <https://pypi.org/project/hidef>`__
* `httplib2 >=0.20.2 <https://pypi.org/project/httplib2>`__
* `matplotlib >=3.5.0 <https://pypi.org/project/matplotlib>`__
* `mygene >=3.2.2 <https://pypi.org/project/mygene>`__
* `ndex2 >=3.5.0 <https://pypi.org/project/ndex2>`__
* `networkx >=2.6.3,<3.0 <https://pypi.org/project/networkx/2.6.3>`__
* `numpy >=1.21.4 <https://pypi.org/project/numpy>`__
* `obonet >=1.0.0 <https://pypi.org/project/obonet>`__
* `pandas >=1.3.4,<2.0 <https://pypi.org/project/pandas/1.3.4>`__
* `requests >=2.26.0 <https://pypi.org/project/requests>`__
* `scikit-learn >=1.0.1 <https://pypi.org/project/scikit-learn>`__
* `scipy >=1.7.2 <https://pypi.org/project/scipy>`__
* `seaborn >=0.13.0 <https://pypi.org/project/seaborn>`__
* `statsmodels >=0.13.5 <https://pypi.org/project/statsmodels>`__
* `tqdm >=4.62.3 <https://pypi.org/project/tqdm>`__

Compatibility
-------------

* Python 3.10+

Citing neteval
--------------

If you use neteval in your research, please cite the following publication:

Wright, SN., et al. *State of The Interactomes: an evaluation of molecular networks for generating biological insights.*


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _NDEx: http://www.ndexbio.org
