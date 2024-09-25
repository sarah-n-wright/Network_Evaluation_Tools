# Interaction Prediction with MPS and L3

Note that implementation of these algorithms with large networks can require up to 100GB of memory. 

## Requirements & Installation

* [L3](https://github.com/kpisti/L3)
* [MPS](https://github.com/spxuw/PPI-Prediction-Project) - Should be installed and run from a separate environment to avoid dependency conflictes
* [GNU parallel](https://www.gnu.org/software/parallel/)

Identify the paths to L3, MPS and parallel executables for use in later steps.  
* For L3: `execdir` - Directory containing `L3.cpp`  
* For MPS: `mps_dir` - Directory containing `topological_feature_extractor`
* For parallel: full path to `parallel` executable  

## Step 1: Create cross-validation folds
Create 10 folds from the largest connected component of a network, where 10% of the interactions are removed for each fold.

Inputs:
* `pref` - Prefix of the network file (`<pref>_net.txt`)
* `datadir` - Directory containing the network file 
* `outdir` - Directory to write the cross-validation networks

Outputs:
* `.fold[n]` files - Tab separated edge list for each fold
* `.fold[n]_test` files - Tab separated edge list of interations held out for each fold

**Usage:**
```
# Run directly:
python ~/Git/Network_Evaluation_Tools/neteval/edge_prediction.py --networkprefix <pref> --runwhat Folds --datadir <datadir> --outdir <outdir>
# Iterate over a list of network prefixes
sh ~/Git/Network_Evaluation_Tools/ExampleUsage/IP_createFolds.sh ~/Git/Network_Evaluation_Tools/Data/example_prefix_file.txt
```

## Step 2: Generate predictions using L3
For cross-validation we run L3 for each of the 10 folds, for external evaluation we run L3 for the full network.

### Run a single prediction

Inputs:
* `<pref>` - Prefix of the network file (`<pref>_net.txt`)
* `<datadir>` - Directory containing the input files
* `<outdir>` - Directory to write the prediction results
* `<execdir>` - Path to the L3 executable

Outputs:
* `<.dat>` - Tab separated file listing all possible interactions, sorted by prediction score.

**Usage:**
```
# To run prediction for fold 1:
python ~/Git/Network_Evaluation_Tools/neteval/edge_prediction.py --networkprefix <pref> --runwhat Predict \
    --networksuffix .fold1 --datadir <datadir> --execdir <execdir> --outdir <outdir> --pred_method L3
```

### Run all predictions for all networks
Inputs:
* `<pref_file>` - File containing the list of network prefixes to be evaluated
* `<l3_execdir>` - Path to the L3 executable
* `<parallel_path>` - Path to the parallel executable

Outputs:
* `...<.dat>` - Tab separated file listing all possible interactions, sorted by prediction score.

**Usage:**
```
sh ~/Git/Network_Evaluation_Tools/ExampleUsage/IP_run_L3.sh <pref_file> <l3_execdir> <parallel_path>
```

## Step 3: Generate predictions using MPS

### Reformat network files to .csv format

Inputs:
* `<pref_file>` - File containing all network prefixes to be evaluated
* `<datadir>` - Directory containing the input files

Outputs:
* `<datadir>/csvs/<pref>.csv` - .csv formated versions of all input network files, included all folds

**Usage:**
```
sh ~/Git/Network_Evaluation_Tools/ExampleUsage/IP_make_csvs.sh <pref_file> <datadir>
```
### Generate Topological Scores using MPS
Inputs:
* `<pref>` - Prefix of the network file (`<pref>_net.txt`)
* `<datadir>` - Directory containing the input files
* `<outdir>` - Directory to write the prediction results
* `<mps_dir>` - Path to the MPS functions
* `<parallel_path>` - Path to the parallel executable

Outputs:
* `.tsv` - Tab separated file containing topological metrics for all possible interactions

**Usage:**
```
# Run a single analysis
cd ${mps_dir}/topological_feature_extractor/bin
java -Xms12g -Xmx12g algos.TopologicalFeaturesExtractor e $datadir/csvs <pref>.csv $outdir 
# Run for all networks and folds
sh ~/Git/Network_Evaluation_Tools/ExampleUsage/IP_run_MPS_topology.sh <pref_file> <datadir> <outdir> <mps_dir> <parallel_path>
```

## Step 4: Evaluate predictive performance for held-out self interactions

Inputs:
* `<pref>` - Prefix of the network file (`<pref>_net.txt`)
* `<datadir>` - Directory containing the input files
* `<outdir>` - Directory to write the prediction results
* `<pred_method>` - "L3" or "MPS"

Outputs:
* `..._results.tsv` -  tab separated file containing AURPC and P@k metrics

**Usage:**
```
# Single Evaluation (e.g. for fold1)
python ~/Git/Network_Evaluation_Tools/neteval/edge_prediction.py --networkprefix <pref> --runwhat EvaluateHeldOut \
        --networksuffix .fold1 --datadir <datadir> --outdir <outdir> --pred_method <pred_method>
# All Evalutations
sh ~/Git/Network_Evaluation_Tools/ExampleUsage/IP_run_heldout.sh <pref_file> <pred_method> <parallel_path>
```

## Step 5: Perform evaluation with external interaciton sets

This analysis should be performed on the prediction results from the full networks.

Inputs:
* `<pref>` - Prefix of the network file (`<pref>_net.txt`)
* `<datadir>` - Directory containing the input files
* `<outdir>` - Directory to write the prediction results
* `<pred_method>` - "L3" or "MPS"

Outputs:
* `...<benchmark>_results.tsv` -  tab separated file containing AURPC and P@k metrics. One file for each benchmark

**Usage:**
```
# Single Evaluation (e.g. for fold1)
python ~/Git/Network_Evaluation_Tools/neteval/edge_prediction.py --networkprefix <pref> --runwhat EvaluateExternal \
        --benchmarks corum panther --datadir <datadir> --outdir <outdir> --pred_method <pred_method>
# All Evalutations
sh ~/Git/Network_Evaluation_Tools/ExampleUsage/IP_run_external.sh <pref_file> <pred_method>
```

## Step 6: Perform evaluation of network coverage of predicted interaction

Inputs:
* `<pref>` - Prefix of the network file (`<pref>_net.txt`)
* `<datadir>` - Directory containing the input files
* `<outdir>` - Directory to write the prediction results
* `<pred_method>` - "L3" or "MPS"
* `<coverage_file>` - Path to the pickle object containing network coverage for all interactions 

Outputs:
* `...coverage_results.tsv` -  tab separated file with the top 100 predictions binned by network coverage
* `...unverified_edges.txt - tab separated edge list of top predicted interactions with network coverage = 0

**Usage:**
```
# Single Evaluation (e.g. for fold1)
python ~/Git/Network_Evaluation_Tools/neteval/edge_prediction.py --networkprefix <pref> --runwhat EvaluateCoverage \
        --coverage_file <coverage_file> --datadir <datadir> --outdir <outdir> --pred_method <pred_method>
# All Evalutations
sh ~/Git/Network_Evaluation_Tools/ExampleUsage/IP_run_external.sh <pref_file> <pred_method> <coverage_file>
```