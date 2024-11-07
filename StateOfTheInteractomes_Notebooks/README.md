# State of the Interactomes: Notebooks

This directory contains code and guidelines for reproducing data and figures contained 
in the manuscript Wright SN, et al. **State of the Interactomes: an evaluation of 
molecular networks for generating biological insights.** *Molecular Systems Biology* (2024). 

Due to the computational requirements of the underlying analyses, these notebooks 
leverage pre-computed data and example implementations with small networks. Much 
of the State of the Interactomes pipeline is designed to run in a high-performance 
computing environment. Please see `ExampleUsage` for guidelines on implementing 
each stage of the pipeline.

`Supplemental_Code` contains additional supplementary analyses.  

## Installation Instructions

To run the State Of the Interactomes Notebooks, install the required dependencies:

```
pip install -r requirements_stateoftheinteractomes.txt
```

To run Supplemental Code, see additional dependencies in the associated README files:
* Gene Function Prediction by GBA (`EGAD_README.md`)
* Processing of Gene Conservation Scores (`phyloP_README.md`)

## Inputs/Outputs

* `StateOfTheInteractomes/Notebooks/Data/` contains pre-computed data for visualization
and analysis, including Datasets EV1-6
* Other data neccessary for analysis is contained in `Data/`
* Generated figures are saved to `StateOfTheInteractomes_Notebooks/Figures`
* Generated data is saved to `Data/example_outputs/`
