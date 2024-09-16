# Summarized gene conservation scores README

Additional requirements:
* BEDOPS v2.4.41

## Step 1: Download data
### Downloading phyloP scores

Positional scores were downloaded from 
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP30way/.  

File names with the format `chr{chr#}.phyloP30way.wigFix.gz`, ignoring all `alt` 
and `random` files.  

### Downloading gene position information

Gene annotations were downloaded from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gff3.gz  


## Step 2: Process phyloP wig files

From the directory containing the downloaded `wigFix.gz` files run:
```
sh /path/to/Git/Network_Evaluation_Tools/StateOfTheInteractomes_Notebooks/Supplemental_Code/phyloP_run_wig_bed.sh
```
This will create a compressed bed file for each chromosome.

## Step 3: Process gene annotation information

From the directory containing the gencode annotations `gene.v46.basic.bed` run:
```
sh //path/to/Git/Network_Evaluation_Tools/StateOfTheInteractomes_Notebooks/Supplemental_Code/phyloP_run_pergene.sh
```
This will create a bed file for each gene

## Step 4: Classify each gene by chromosome
From the `perGene` directory created in step 3, create a list of all CDS files:
```
find . -maxdepth 1 -type f -name 'CDS:*' | xargs -n 1 >> all_CDS.txt
```

Then classify all genes based on their chromosome:
```
sh phyloP_classify_chrs.sh all_CDS.txt
```
This will create lists of chromosome files chr{chr#}_files.txt

## Step 5: Summarize the phyloP scores per gene CDS

For each chromosome run:
```
# e.g. for chromosome 22
sh phyloP_run_CDS.sh 22
```
This will loop through every gene assigned to the chromosome and summarize the  
conservation scores per gene using `bedmap`, creating [GENE].bed.summary.bed

## Step 6: Aggregate CDS for each gene

First, from the data directory create a list of all the summarized files:

```
find . -maxdepth 1 -type f -name '*summary.bed' | xargs -n 1 >> all_summary_files.txt
```
Then run the python script `phyloP_aggregation.py` to generate 
`all_summarized_scores.txt`
```
python phyloP_aggregation.py /path/to/data/ all_summary_files.txt
```
