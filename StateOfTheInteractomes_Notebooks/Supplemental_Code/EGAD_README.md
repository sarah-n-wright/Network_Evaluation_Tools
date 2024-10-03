# Gene Function Prediction by GBA

## Requirements

* R v4.4.1
* [EGAD v1.32.0](https://bioconductor.org/packages/release/bioc/html/EGAD.html)
* data.table >= v1.15.4

## Step 1: Create input GO sets

GO - gene associations should be formatted as:

```
GO      Branch  GeneID
GO:0008038      Process 2048
GO:0008038      Process 2049
GO:0008038      Process 152330
GO:0008038      Process 54538
GO:0008038      Process 8851
GO:0008038      Process 1949
GO:0008038      Process 4897
GO:0008038      Process 1954
GO:0008038      Process 1826
```
Where each row corresponds to a single association between a GO term and an NCBI Gene ID

## Step 2: Run EGAD for all networks
Perform the guilt-by-association analysis for the example networks.

* `datadir` - directory containing network files, go data, and outputs
* `prefix_file` - filepath for list of network prefixes
* `stat` - 'PR' for AUPRC, 'AUROC' for AUROC
* `min_genes` - Minimum intersection between network genes and gene set
* `max_genes` - Maximum intersection between network genes and gene set
* `go_file` - Filepath to go association data

```
#sh ega_run.sh <prefix_file> <stat> <min_genes> <max_genes> <datadir> <go_file>
sh egad_run.sh ../../Data/example_prefix_file.txt \ 
	PR 20 250 ../../Data/example_outputs/ \
	../../Data/go_dataframe_for_EGAD.tsv
```
Prediction performance metrics will be written to the file:
`<datadir>/<pref>_GO_min<min_genes>_max<max_genes>_cv5.<stat>.roc.txt`

