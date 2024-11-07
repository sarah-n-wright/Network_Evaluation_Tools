# Perform gene set recovery evaluation

## Pre-requisites

* Establish an account at [DisGeNET.com](disgenet.com)
* Processed networks (See `NP_NetworkProcessing.README.md`)
* If evaluating large networks (>500K interactions), at least 64GB of available memory

All example code below is intended to be run from the base directory of `Network_Evaluation_Tools/`. All python scripts should be installed automatically in bin/ after installaion via PyPi.  

## Step 1: Download DisGeNET and GWAS Catalog data

### DisGeNET

Associations are downloaded from [DisGeNET](DisGeNET.com) using the API

Inputs:
* `min_genes` - Minimum number of gene associations to report a trait (default 5)
* `max_genes` - Maximum number of gene associations to report a trait (default 500)
* `outdir` - path to save the downloaded file and output genesets 
* `email` - Email of DisGeNET account for API
* `password` - Password of DisGeNET account for API
* `source` - Dataset to source gene-disease associations. See [DisGeNET Sources](https://disgenet.com/About#sources) (default 'BEFREE')
* `disgenfile` - Previously downloaded file to process. If not provided, one will be downloaded. (Optional)

Outputs:
* `<outdir>disease_list_<source><current date>.tsv` - The downloaded associations
* `<outdir>disease_list_<source><current date>.tsv.genesets` - The genesets derived from associations

**Usage:**
```
python get_disgen_associations.py -m <min_genes> -M <max_genes> -o <outdir> -E <email> -P <password> -S <source> -d <disgenfile>
```

### GWAS Catalog

All associations are downloaded from the GWAS catalog from [https://www.ebi.ac.uk/gwas/docs/file-downloads](https://www.ebi.ac.uk/gwas/docs/file-downloads) (All associations **v1.0.2 - with added ontology annotations, GWAS Catalog study accession numbers and genotyping technology**):

**Download, clean and process the GWAS data**

Inputs:
* `min_genes` - Minimum number of gene associations to report a trait (default 5)
* `max_genes` - Maximum number of gene associations to report a trait (default 500)
* `outdir` - path to save the downloaded file and output genesets 
* `p_th` - P-value threshold for accepting an association (default 5e-8)
* `split_date` - (optional) Cut off date to create before/after gene sets, 'YYYY-MM-DD' format
* `gwas_file` - (Optional). If not provided, data will be downloaded.

Outputs:
* `gwas_catalog_{current_date}.txt` - Downloaded associations
* `<gwas_file>.genesets` - Gene sets derived from associations

**Usage:**
```
python get_gwas_associations.py -m <min_genes> -M <max_genes> -o <outdir> -p <p_th> -D <split_date> -G <gwas_file> 
```

## Step 2: Preprocess downloaded genesets

Geneset processing performs:
 1. Conversion of gene identifiers to NCBI Gene IDs if needed
 2. Filters the genesets against all networks to ensure all networks have a minimum number of overlapping genes. Strict filtering requires that the same minimum set of genes is present in all networks. 

Inputs:
* `geneset_file` - A `.genesets` file from Step 1.
* `prefix_file` - File containing the list of network prefixes to base filtering on
* `output_file` - Full file path for writing the final genesets
* `datadir` - Directory containing the processed networks
* `id_type` - Input gene identifier type in `geneset_file`
* `min_genes` - (default 10)
* `max_genes` - (default 500)
* `-C` - Flag to perform conversion of gene set identifiers.
* `-F` - Flag to perform filtering of genesets based on min/max genes
* `--strict` - Flag to perform strict filtering. 

Outputs:
* `output_file` - Filtered and converted genesets

**Usage:**
```
python prepare_evaluation_data.py -s <geneset_file> -n <prefix_file> -o <out_file> -d <datadir> -i <id_type> \ 
                -m <min_genes> -M <max_genes> [-F] [-C] [--strict]
```


## Step 3: Run gene set recovery evaluation

Inputs:
* `net_path` - Full path to the network file to be analyzed 
* `node_sets_file` - Full path to the file of genesets
* `nshuffle` - Number of network shufflings to perform (Default 30)
* `nsamples` - Number of geneset samples to perform (Default 30)
* `alpha` - Network propagation constant. If None, will be fit based on parameter optimization. 
* `sampp` - Subsampling proportion for genesets. If None, will be fit based on paramter optimization
* `min_genes` - Minimum number of genes to evaluate each gene set
* `outdir` - Directory for writing output files

Outputs:
* `<pref>.<set>.aurpcs.csv` - Acutal observed auprcs for all sets
* `<pref>.<set>.null_aurpcs.csv` - Null auprcs from shuffled networks
* `<pref>.<set>.performance.csv` - Performance scores for all sets
* `<pref>.<set>.performance_gain.csv` - Performance gain results for all sets
* `shuffNet_[n]_<pref>_net_shuffled.txt` for `n` in `1:nshuffle` - Randomized network shufflings


**Usage:**
To run a single evaluation with a network a geneset file:
```
python run_network_evaluation.py \
    --cores 2 -i <nshuffle> -n <nsamples> -a <alpha> -p <sampp> --min_genes <min_genes> \
    -o <outdir> <net_path> <node_sets_file>
```
e.g. for DIP with DisGeNET genesets:
```
python run_network_evaluation.py \
    --cores 2 -i 50 -n 50 -a 0.64 -p 0.3 --min_genes 20 \
    -o Data/example_outputs/GeneSetRecovery/ Data/example_outputs/dip_net.txt Data/disgen.genesets
```
To run analysis with all included genesets and default parameters, run:

Inputs:
* `pref` - The prefix of the example network to run
```
sh ExampleUsage/GSR_run_evaluation.sh <pref>
```