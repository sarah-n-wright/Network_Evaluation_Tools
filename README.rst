============================
Network Evaluations Tools v2
============================


.. image:: https://img.shields.io/pypi/v/neteval.svg
        :target: https://pypi.python.org/pypi/neteval

.. image:: https://readthedocs.org/projects/neteval/badge/?version=latest
        :target: https://neteval.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


Network Evaluation Tools 2 is a Python 3.10 package accompanying the manuscript:
*State of the Interactomes: an evaluation of molecular networks for generating biological insights.* Sarah N Wright et al. biorxiv.org (2024).

This repository contains all python code used in the study, as well as example usage on a small network example. 
This package updates and expands the work developed as part of 
`Huang and Carlin et al. 2018 <http://www.cell.com/cell-systems/fulltext/S2405-4712(18)30095-4>`_.


Installation
============
To install the package, run the following command:

   .. code-block:: bash

      pip install neteval

To install the package from source, clone the repository and run the following command:

   .. code-block:: bash
   
      git clone https://github.com/sarah-n-wright/Network_Evaluation_Tools
      cd Network_Evaluation_Tools
      make dist
      pip install dist/netevalcmd*whl

For example usage of command line scripts see `Example Usage <https://github.com/sarah-n-wright/Network_Evaluation_Tools/ExampleUsage>`__.  
For example usage of all other funcitonality see `State of the Interactomes Notebooks <https://github.com/sarah-n-wright/Network_Evaluation_Tools/StateOfTheInteractomes_Notebooks>`__.  

Optional External Packages
--------------------------

Interaction Prediction algorithms
"""""""""""""""""""""""""""""""""
* **L3:** The L3 algorithm can be installed from `https://github.com/kpisti/L3 <https://github.com/kpisti/L3>`_. The path to the executable can then be used with `edge_prediction.py`.

* **MPS(T):** The MPS(T) algorithm can be installed from `https://github.com/spxuw/PPI-Prediction-Project <https://github.com/spxuw/PPI-Prediction-Project>`_, and should be installed and run from a separate environment.

AlphaFold Multimer
""""""""""""""""""
* **AlphaFold-Multimer** can be installed from `https://github.com/YoshitakaMo/localcolabfold <https://github.com/YoshitakaMo/localcolabfold>`_, and should be installed and run from a separate environment.

Dependencies
============

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


Modules in this package
=======================

Interactome processing, standardization and annotation
------------------------------------------------------
Modules:

* ``processing_functions``: Module utilized by process_data, contains the NetworkData class.__
* ``gene_mapper``: Includes query_hgnc, query_ensembl, query_uniprot. Modules to map identifiers to NCBI Gene IDs from a variety of identifier types including: Uniprot, Ensembl, EnsemblProtein, RefSeq, Symbol. Incorporates methods to handle out of date identifiers.
* ``node_annotation``: Module for downloading and extracting gene annotation data from HGNC, NCBI, Uniprot and Ensembl. Includes the class ExpressionData for loading and analyzing mRNA and protein expression data.
* ``network_statistics``: Module to extract summary statistics for a set of network, such as node and edge counts.
* ``gsea_functions``: Module for performing downloading, processing and analyzing Gene Ontology data.

**Script(s):**

* ``process_data.py``: Python script to process a raw interactome based on a configuration file

For detailed usage see ``ExampleUsage/NP_NetworkProcessing.README.md``

Gene set recovery performance evaluation
----------------------------------------

Evaluation data collection and processing
"""""""""""""""""""""""""""""""""""""""""

Modules:

* ``get_disgen_associations``: Module to download DisGeNET associations and generate genesets
* ``get_gwas_associations``: Module to download GWAS Catalog associations and generate genesets

**Script(s):**

* ``prepare_evaluation_data.py``: Python script to convert gene identifiers within gene sets and filter based on network coverage.

For detailed usage, see ``ExampleUsage/GSR_GeneSetRecovery.README.md``

Gene set recovery
"""""""""""""""""

Modules:

* ``network_evaluation_functions``: Module for performing and evaluating gene set recovery.
* ``network_propagation``: Underlying network propagation methodology.
* ``shuffle_networks``: Module for creating degree-matched shuffled networks
* ``gene_set_recovery_results``: Module to load, evaluate, and visualize gene set recovery results. Includes the class EvaluationResults.

**Script(s):**

* ``run_network_evaluation.py``: Python script to perform gene set recovery performance evaluation

For detailed usage, see ``ExampleUsage/GSR_GeneSetRecovery.README.md``

Parsimonious Composite Networks (PCNets)
-----------------------------------------

**Script(s):**

* ``network_constructor``: Python script to create composite networks using the *global composite* and *ranked composite* approaches. See ExampleUsage/run_composite.sh.

For detailed usage see ``ExampleUsage/PC_PCNets.README.md``

Interaction & complex prediction
--------------------------------

Modules:

* ``community_annotation``: Module for assessing the quality of gene communities in a network.
* ``edge_prediction``: Module for performing and analyzing edge prediction results.

**Script(s):**

* ``edge_prediction.py``: Script for performing edge prediciton evaluation.
* ``alphafold_results.py``: Script for parsing and analyzing AlphaFold results.
* ``complex_evaluation.py``: Script for evaluating hierarchical complex prediction results.

For detailed usage see ``ExampleUsage/IP_InteractionPrediction.README.md`` and ``ExampleUsage/AF_AlphaFold.README.md``

General utilities
-----------------

* ``data_import_export_tools``: Module of functions for importing and exporting the various data formats used by this package.
* ``Timer``: Class that measures the elapsed time of various processing steps and outputs a summary.


Provided Data and Implementation Examples
=========================================

ExampleUsage
------------

This directory contains README and bash scripts for implemenation of each stage of the network evaluation pipeline. 
All examples utilize three small interactomes (DIP, PID2, and Wan). While most of the pipeline is designed to run in a 
high-performance computing environment, most of these examples can be run on a local machine.

* ``NP_NetworkProcessing.README.md``
* ``GSR_GeneSetRecovery.README.md``
* ``PC_PCNets.README.md``
* ``IP_InteractionPrediction.README.md``
* ``AF_AlphaFold.README.md``

Data
----

This directory contains key data sets used for the evaluation of interactomes for prosperity, including: 

* Annotation data from HGNC, Ensembl, NCBIm, and Uniprot
* Gene sets analyzed
* Gene conservation scores
* Example networks (Wan, DIP, PID2)
* CORUM and PANTHER edge lists


StateOfTheInteractomes_Notebooks
--------------------------------

This directory contains code and guidelines for reproducing data and figures contained
in the manuscript.

Notebooks
"""""""""

* 1_Statistics_and_Representation.ipynb
* 2_GO_analysis.ipynb
* 3_Gene_Set_Recovery.ipynb
* 4_Composite_networks.ipynb
* 5_Interaction_and_Complex_Prediction.ipynb
* 6_AlphaFold_Assessment.ipynb

Due to the computational requirements of the underlying analyses, these notebooks
leverage pre-computed data and example implementations with small networks. Much
of the State of the Interactomes pipeline is designed to run in a high-performance
computing environment. Please see ``ExampleUsage`` for guidelines on implementing
each stage of the pipeline.

To run the State Of the Interactomes Notebooks, install the required dependencies:

   .. code-block:: bash

        pip install -r requirements_stateoftheinteractomes.txt


**Inputs/Outputs**

* `StateOfTheInteractomes_Notebooks/Data/` contains pre-computed data for visualization
* Other data neccessary for analysis is contained in `Data/`
* Generated figures are saved to `StateOfTheInteractomes_Notebooks/Figures`
* Generated data is saved to `Data/example_outputs/`

StateOfTheInteractomes_Notebooks/Data
"""""""""""""""""""""""""""""""""""""

This directory contains data necessary for recreating the manuscript figures, including *Supplemental Tables 2-5,7-8*, and other precomputed results.

StateOfTheInteractomes_Notebooks/Supplemental_Code
""""""""""""""""""""""""""""""""""""""""""""""""""

This directory contains code used in the generation of the manuscript results that is not included in the primary ``neteval`` package. 
This includes implementation of EGAD (Extending Guilt by Association by Degree) and HiDeF (Hierarchical community Decoding Framework), as well 
as processing of PDB files and gene conservation scores. 

To run Supplemental Code, see additional dependencies in the associated README files:

* Gene Function Prediction by GBA (``EGAD_README.md``)
* Processing of Gene Conservation Scores (``phyloP_README.md``)

Compatibility
=============

* Python 3.10+


Citing neteval
==============

If you use neteval in your research, please cite the following publication:

Wright, SN., et al. *State of The Interactomes: an evaluation of molecular networks for generating biological insights.*


Credits
=======

This package is built from the original `Network Evaluation Tools <https://github.com/idekerlab/Network_Evaluation_Tools>`_ developed by `Huang and Carlin et al. 2018 <http://www.cell.com/cell-systems/fulltext/S2405-4712(18)30095-4>`_.


This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _NDEx: http://www.ndexbio.org
