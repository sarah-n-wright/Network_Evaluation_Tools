# In Silico Interaction Assessment using AlphaFold-Multimer

## Requirements

* [localcolabfold](https://github.com/YoshitakaMo/localcolabfold)

It is recommended to install `localcolabfold` in a separate environment from `neteval`.

All example code below is intended to be run from the base directory of `Network_Evaluation_Tools/`. All python scripts should be installed automatically in bin/ after installaion via PyPi.

## Step 1: Download protein sequences
Environment: `neteval`

Using neteval, download UniProt protein sequencesm, where `gitdir` is the path of the `Network_Evaluation_Tools` repo. Note that this can take up to 15 minutes.

```
from neteval.node_annotation import get_uniprot_annotation_data

uni_df = get_uniprot_annotation_data(fields=['id', 'sequence', 'xref_geneid'], outdir=os.path.join(gitdir, 'Data/'), 
                                     index_on=2, page_size=500, max_retries=5, taxid=9606, verbose=False)
```

## Step 2: Create list of protein pairs to assess

Implementation of AlphaFold-Multimer is computationally expensive, so careful prioritization of protein pairs is advised. Our AlphaFold-Multimer evaluations were conducted using NVIDIA Tesla V100 GPUs operating under Rocky Linux release 8.9. Each evaluation was allocated one CPU and 14 GB of memory and was subject to a 5-hour time limit. Most protein pairs completed evaluation in under 1 hour. 

Protein pairs to be assessed should be saved in a comma separated edgelist, e.g. 
```
> pairs.csv
10667,85015
142,5422
2017,2241
3091,5594
```
## Step 3: Create fasta files

Inputs:
* `pair_file` - full file path of the comma separated edgelist of protein pairs to be evaluated
* `ref_file` - full file path of the downloaded UniProt sequence data

Outputs:
* `<ProteinA>_<ProteinB>.fasta` for each pair of proteins. 

To extract the protein sequences from the downloaded UniProt file and create a `.fasta` file, run:
```
sh ExampleUsage/AF_fasta_pairs.sh <pair_file> <ref_file>
```

## Step 4: Perform AlphaFold-Multimer Analysis
Environment: `colabfold`
Inputs:
* `input_fasta` - `.fasta` file for a protein pair as generated in Step 3
* `outpath` - Directory for saving all outputs of the AF modeling

Outputs:
AlphaFold-Multimer produces a set of files detailing the MSA, modeling, and QC results. The ipTM and pTM scores are contained in the `.json` files
e.g. `<ProteinA>_<ProteinB>_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json`

**Usage:**
```
colabfold_batch --num-recycle 3 --model-type alphafold2_multimer_v3 <input_fasta> <outpath>
```

## Step 5: Assess results
Environment: `neteval`

Extract the pTM and ipTM scores, and calculate the model confidence for all successfully evaluated pairs

Inputs:
* `pair_file` - full file path of the comma separated edgelist of protein pairs to be evaluated
* `datadir` - directory containing the model outputs
* `outdir` - directory for saving the extracted results

Outputs:
* `<outdir>/<pair_file>_alphafold_results.tsv`

**Usage:**
```
python alphafold_results.py <pair_file> <datadir> <outdir>
```