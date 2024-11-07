# Create Composite Networks

## Pre-requisites
* For the creation of *Ranked* composites, you will need a metric to rank the networks (e.g. from Gene Set Recovery)

All example code below is intended to be run from the base directory of `Network_Evaluation_Tools/`. All python scripts should be installed automatically in bin/ after installaion via PyPi.

## Step 2: Create *Global* composite networks

Inputs:
* `pref_file` - File listing the prefixes to be considered
* `n` - minimum network coverage for an interaction to be included. Can be a series of values
* `nodepref` - Prefix of the node columns in the network files
* `m` - 'parsimonious' for global composites
* `netdir` - Directory containing the network files
* `outdir` - Directory to save the networks
* `name` - Prefix for naming the created networks

Outputs:
* `<name>_composite_min<n>.nodeslist` - Composite network node list
* `<name>_composite_min<n>_net.txt` - Composite network edgelist 

**Usage:**
```
python network_constructor.py -n <n> -d <netdir> -o <outdir> -m parsimonious --name <name> --nodepref <nodepref> <pref_file>
```
E.g. for the example prefixes
```
python network_constructor.py -n 1 2 3 -d Data/example_outputs -o Data/example_outputs/Composites -m parsimonious \
    --name Global --nodepref Entrez_ Data/example_prefix_file.txt
```

## Step 3: Create *Ranked* composite networks

**Rank the networks to be considered**
Create a file listing the network prefixes in ranked order. For example, in the manuscript we ranked the networks based on mean size-adjusted performance across Literature and Genetic gene sets. 

**Create the composites**

This method iterates over the number of top ranked networks (k) (i.e. top-2, then top-3), creating composites where each interaction is present in at least `n` of the top-ranked networks. 

Inputs:
* `pref_file` - File listing the prefixes to be considered in ranked order
* `n` - minimum network coverage for an interaction to be included. Can be a series of values
* `nodepref` - Prefix of the node columns in the network files
* `m` - 'ordered' for ranked composites
* `netdir` - Directory containing the network files
* `outdir` - Directory to save the networks
* `name` - Prefix for naming the created networks

Outputs:
* `<name>_<n>_composite_min<k>.nodeslist` - Composite network node list
* `<name>_<n>_composite_min<k>_net.txt` - Composite network edgelist 

**Usage:**
```
python network_constructor.py -n <n> -d <netdir> -o <outdir> -m ordered --name <name> --nodepref <nodepref> <pref_file>
```
E.g. for the example prefixes (assuming these are in rank order)
```
python network_constructor.py -n 2 -d Data/example_outputs -o Data/example_outputs/Composites -m ordered \
    --name Ranked --nodepref Entrez_ Data/example_prefix_file.txt
```

Composites can also be run from the script:
```
sh ExampleUsage/PC_create_composites.sh Data/example_prefix_file.txt Global global
```