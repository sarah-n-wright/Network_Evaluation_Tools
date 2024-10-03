# Hierarchical Complex Prediction with HiDeF

## Additional Requirements

* [hidef v1.1.5](https://pypi.org/project/hidef/)

## Step 1: Identify the path to the file `hidef_finder.py`

It is usually found in environment bin e.g.  
`.../anaconda3/envs/<env_name>/bin`

## Step 2: Run HiDeF

The provided file `hidef_run.sh` first reformats the network files to be compatible with HiDeF, then models the hierarchy.

Inputs:
* `network_pref` - prefix of the network to be evaluated
* `maxres` - maximum resolution parameter for HiDeF. Higher values will generate more communities.
* `datadir` - directory containing the processed network files
* `hidef_path` - directory containing `hidef_finder.py` 

Outputs:
* `.nodes` file - All communities predicted, including size and community members 
* `.edges` file - Relationships between communities 

**Usage:**
```
sh hidef_run.sh <network_prefix> <maxres> <datadir> <hidef_path>
```
