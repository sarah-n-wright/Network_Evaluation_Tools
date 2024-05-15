# Network Evaluation Tools 2

Network Evaluation Tools 2 is a Python 3.10 package accompanying the manuscript:
*State of the Interactomes: an evaluation of molecular networks for generating biological insights.* Sarah N Wright et al. biorxiv.org (2024). 

This repository contains all python code used in the study, as well as example usage on a small network example. This package updates and expands the work developed as part of [Huang and Carlin et al. 2018](http://www.cell.com/cell-systems/fulltext/S2405-4712(18)30095-4). 

## Modules in this package

### Interactome processing, standardization and annotation
* `process_data`. Python script to process a raw interactome based on a configuration file
* `processing_functions`. Module utilized by `process_data`, contains the `NetworkData` class.
* `gene_mapper`. Includes `query_hgnc`, `query_ensembl`, `query_uniprot`. Modules to map identifiers to NCBI Gene IDs from a variety of identifier types including: Uniprot, Ensembl, EnsemblProtein, RefSeq, Symbol. Incorporates methods to handle out of date identifiers.
* `node_annotation`. Module for downloading and extracting gene annotation data from HGNC, NCBI, Uniprot and Ensembl. Includes the class `ExpressionData` for loading and analyzing mRNA and protein expression data.
* `network_statistics`. Module to extract summary statistics for a set of network, such as node and edge counts.
### Gene set recovery performance evaluation
#### Evaluation data collection and processing
* `get_disgen_associations` Module to download DisGeNET associations and generate genesets
* `get_gwas_associations` Module to download GWAS Catalog associations and generate genesets
* `prepare_evaluation_data` Python script to convert gene identifiers within gene sets and filter based on network coverage. 
#### Gene set recovery
* `network_evaluation_functions` Module for perfoming and evaluating gene set recovery.
* `network_propagation` Underlying network propagation methodology.
* `run_network_evaluation` Python script to perform gene set recovery performance evaluation
* `shuffle_networks` Module for creating degree-matched shuffled networks
* `gene_set_recovery_results` Module to load, evaluate, and visualize gene set recovery results. Includes the class `EvaluationResults`

### Parsimonious Composite Networks (PCNets)
* `network_constructor` Python script to create composite networks using the *global composite* and *ranked composite* approaches. See `ExampleUsage/run_composite.sh`.

### Interaction prediction performance evaluation
* `L3_prediction.sh` Bash script to perform L3 edge prediction, called by neteval.edge_prediction.l3_edge_prediction
* `edge_prediction`---------
* 


### General utilities
* `data_import_export_tools` Module of functions for importing and exporting the various data formats used by this package. 
* `Timer` Class that measures the elapsed time of various processing steps and outputs a summary.

   
### ??
* `network_filter`
* `node_permute_analysis`
* `performance_comparison`

## Example Usage
This repository contains an example implementation of the evaluation pipeline using the Database of Interacting Proteins (DIP). The raw interactome file is contained in `Data/PathwayCommons12.dip.hgnc.txt.clean.gz`.
### 1. Network Processing
Processing parameters are defined in `data_configs/dip.config`
```
### DATA ARGUMENTS ###
## REQUIRED ##

# name to be given to the network
name=dip
# full filepath of the input data
input_datafile=../Data/PathwayCommons12.dip.hgnc.txt.clean.gz
# full directory for saving processed data
outpath=../Data/example_outputs/
...
```
The network is standardized using `process_data_job.sh`. This converts all gene identifiers to NCBI Gene IDs (Entrezgenes), removes duplicate and self interactions, and filters based on score and species (where applicable). See `ExampleUsage/process_data_job.sh`.  
```
sbatch process_data_job.sh dip.config 0
```
Three files are created in `/Data/example_outputs/`:
* `dip_net.txt` - Processed edge list with gene identifiers, index in original datafile, and score (if applicable)
* `dip.nodelist` - List of unique nodes in the interactome.
* `dip.stats` - Statistics showing the number of nodes and edges remaining after each processing step.

### 2. Download and preparation of gene sets from DisGeNet and the GWAS catalog.

**DisGeNet**  
First create an account at DisGeNET, then download associations and create a `.genesets` file.
```
python neteval/get_disgen_associations.py -m 10 -M 500 -o <output_directory> -E <email> -P <password> -S BEFREE -u
```
Where:  
* `m`: minimum number of genes per set
* `M`: maximum number of genes per set
* `S`: source of gene associations ([https://www.disgenet.org/dbinfo])
* `u`: update flag to force download of latest associations, even if file already exists in output_directory

**GWAS Catalog**  
Download association data from GWAS Catalog, this will create a `.genesets` file
```
python neteval/get_gwas_associations.py -m <min_genes> -M <max_genes> -o <output_directory> -D <split_date>
```
**Format and Filter**  
Convert gene identifiers to NCBI Gene IDs, and filter based on network coverage of (a set of) interactome(s).
```
python neteval/prepare_evaluation_data.py -s <setfile> -i <id_type> -n <prefix_file> -F -C -m <min_genes> -M <max_genes> -o <output_file> -n <network_directory>
```
### 3. Gene Set Recovery
Perform gene set recovery evaluation using network propagation. See `ExampleUsage/run_single_evaluation.py`.
```
python neteval/run_network_evaluation.py -i <number_of_shuffles> -n <number_of_subsamples> -a <alpha> -p <subsampling_parameter> \
                                        --min_genes <min_geneset_genes> --null_network_outdir <path_for_null_networks>
                                        --null_AUPRCs_save_path <path_for_null_AUPRCs>
                                        --performance_save_path <path_for_performance_results>
                                        --performance_gain_save_path <path_for_performance_gain_results>
                                        <network_file> <geneset_file> <path_for_AUPRC results>
```
This will generate shuffled versions of the network and calculate gene set recovery statistics for the network against the shuffled distribution. 

Load gene set recovery results
```
from neteval.gene_set_recovery_results import *
from neteval.network_statistics import *
# Load the evaluation results
eval_results = EvaluationResults('Data/example_outputs/', 'Data/example_prefix_file.txt', metrics=['Performance', 'AUPRC'], genesets=['disgen', 'gwas'])
# Load the network statistics
stats = NetworkStats('Data/example_prefix_file.txt', 'Data/example_outputs/')
# Perform size adjustment and ranking
eval_results.size_adjsuted_performances(sizes=stats.edge_counts) # requires >1 network evaluated
eval_results.rank_all(na_option='center') # requires >1 network evaluated
```
### 4. Interaction prediction
**L3**
* Prediction and evaluation were performed via `neteval/edge_prediction.py`, utilizing the L3 algorithm from [https://github.com/kpisti/L3].
* See `ExampleUsage/run_L3_prediction.sh`

**MPS(T)** (installed and run separately)
* Networks files were first converted to csv using `neteval/net_to_csv.sh`.
* Topology calculations were performed using the 'external' topological feature extractor across 10 network folds. See `ExampleUsage/run_MPS_topology.sh`.
* ```java -Xms12g -Xmx12g algos.TopologicalFeaturesExtractor e ...```
* For consistency with evaluation metrics with L3, evaluation of predicted results was performed using `neteval/edge_prediction.py`, rather than the inbuilt evaluation functions. See `ExampleUsage/run_MPS_prediction.sh`

## Version and Dependencies
neteval requires python 3.10+, and:
  - Argparse >= 1.4.0
  - Httplib2 >= 0.20.2
  - Mygene >= 3.2.2
  - NetworkX >= 2.6.3, may not be compatible with 3+
  - Ndex2 >= 3.5.0
  - Numpy >= 1.21.4
  - Matplotlib >= 3.5.0
  - Pandas >= 1.3.4
  - Requests >= 2.26.0
  - Scipy >= 1.7.2
  - Seaborn >= 0.13.0
  - statsmodels >= 0.13.5
  - tqdm >= 4.62.3
  - Scikit-learn >= 1.0.1

## Installation
### `neteval` package
1. Clone the repository 
2. cd to new respository
3. Execute following command:  
```python setup.py install```
Note: pypi setup in progress

### Interaction prediction methods
**L3**
The L3 algortihm was installed from [https://github.com/kpisti/L3]. The path to the executable can then be used with `edge_prediction.py`

**MPS(T)**
The MPS(T) algortihm was installed from [https://github.com/spxuw/PPI-Prediction-Project], and should be installed and run from a separate environment.  

### AlphaFold Multimer
AlphaFold-Multimer was installed from [https://github.com/YoshitakaMo/localcolabfold], and should be installed and run from a separate environment.

## Data provided in this repository (see ```Data``` Folder)

  - _data_import_tools_ - This module contains functions for helping import network files and gene set files for analysis.
  - _gene_conversion_tools_ - This module contains functions for helping convert, filter, and save networks from their raw database form. Used in the Network Processing Jupyter Notebooks.
  - _miscellaneous_functions_ - This module contains various functions developed to help with analysis along the way. These functions are not well tested and may contain bugs. These functions were generally used to determine other network performance metrics on network recovery of gene sets.
  - _network_evaluation_functions_ - This module contains many of the core functions of the set-based network evaluation algorithm.
  - _network_propagation_ - This module contains functions to help with network propagation steps used in the set-based network evaluation algorithm.

## Version and Dendencies
Currently, the network_evaluation_tools package requires Python 2.7 - Python 2.7.13. Note that some functions in this package may not work with Python 3.0+.
network_evaluation_tools requires: 
  - Argparse >= 1.1
  - NetworkX >= 2.1
  - Numpy >= 1.11.0
  - Matplotlib >= 1.5.1
  - Pandas >= 0.19.0
  - Requests >= 2.13.0
  - Scipy >= 0.17.0
  - Scikit-learn >= 0.17.1

## Issues
Please feel free to post issues/bug reports. Questions can be sent to snwright@ucsd.edu

## License
See the [LICENSE](https://github.com/snwright/Network_Evaluation_Tools/blob/master/LICENSE.txt) file for license rights and limitations (MIT).


