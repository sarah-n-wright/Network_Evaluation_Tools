# Processing and Standardization of Networks

## Requirements

* 

## Step 1: Create config file(s) from template

To enable processing of a large number of networks, we create config files to assign the parameters for network processing. 

E.g. for DIP:
```
cd Network_Evaluation_Tools/data_configs
cp template.config dip.config
```
Open the file `dip.config` in a text editor and fill out the required arguments:
* `name` - prefix to be assigned to all downstream files for the network
* `input_datafile` - full filepath for downloaded network
* `outpath` - directory path for saving the processed data
* `nodeA_col` - name of the column with first interactor
* `nodeB_col` - name of the column with the second interactor
* `input_id_type` - Identifier type(s) (e.g. UniProt, Symbol)

For example, for DIP we create dip.config with:
```
name=dip
input_datafile=../Data/PathwayCommons12.dip.hgnc.txt.clean.gz
outpath=../Data/example_outputs/
nodeA_col=PARTICIPANT_A
nodeB_col=PARTICIPANT_B
input_id_type="Symbol Uniprot"
```

## Step 2: Run the processing pipeline

The processing pipeline performs the following tasks:
 1. **Data Cleaning**
    a. Binarize complexes
    b. Drop interactions with missing node information
    c. Sort the node pairs
    d. Subset to human only interactions (if species columns are specified)
    e. Process scores to floats, and -log10P if score is p-value (if score column is specified)
    e. Drop duplicate interactions (keeping highest score if score column is specified)
    f. Extract gene identifiers
    g. Drop duplicates after indentifier extraction
    h. Convert all identifiers to strings
 2. **Convert node identifiers to NCBI Gene IDs (Entrez), iterating over all specified indentifier types**
    a. Query input type database to ensure identifier is latest version
    b. Map updated indentifiers to NCBI Gene IDs.
    c. Report the unmapped nodes
    d. Convert all indentifiers in the network, dropping interactions where one or both genes could not be mapped
 3. **Remove duplicates and self interactions after identifier mapping**
 4. **Write the processed data to file, including a score susbet if specified.**

To check the set up of the processing pipeline, first run with the flag `--testmode 1` to process only the first 10000 interactions
```
sh Network_Evaluation_Tools/NP_process_data.sh dip.config 1
```
Then run the full processing
```
sh Network_Evaluation_Tools/NP_process_data.sh dip.config 0
```

## Step 3: Assess the outputs

The network processing pipeline generates the following outputs, where `<pref>` is the name assigned to the network:

* `<pref>_net.txt` - The processed edge list, including score (if available) and index in the original dataset
* `<pref>.nodelist` - A list of all distinct genes present in the processed network
* `<pref>_unmapped_nodes.txt` - A list of gene identifiers that could not be mapped to NCBI Gene IDs.
* `<pref>.stats` - Statistics of interactions and genes removed through the processing pipeline. 

To ensure appropriate processing of a network we recommend checking:

**`<pref>_unmapped_nodes.txt`**. Expect some unmapped nodes due to Withdrawn entries etc. However, a large percentage of unmapped nodes could indicate a processing issue.
1. Check for additional identifier types not specified in the `.config` file
2. Check for prefixes not properly handled and update `node_prefix_separator` in `.config` as needed

**`<pref>.stats`** - This file reports the size of the network after various processing steps. A large decrease at any step should be investigated to ensure that it is consistent with expectations. 

```
                    edges   nodes
raw                 16655.0 3464.0    # initial input
Node prefix         16655.0           # after parsing node prefixes
de-duped1           16655.0           # after first de-duplication step
mapped              16460.0 3433.0    # after mapping gene indentifiers
de-duped_final      16460.0           # after second de-duplication
removed_self_edges  16460.0           # after removing self edges
unmapped                    31.0      # total number of unmapped nodes
```