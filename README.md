# Sample Organization Does Not Affect Experience-Based Risky Choice Regardless of Option Complexity	

[![DOI](https://zenodo.org/badge/1094932157.svg)](https://doi.org/10.5281/zenodo.24109605)

	
Research on the sampling paradigm of decisions from experience has documented a correlational link between switch rates during predecisional search and the final choice.
Yet, the causal direction of this link is ambiguous. 
Here, we present an experiment ($N=199$) in which we examined whether different predecisional search patterns (sampling with infrequent switching, with frequent switching, and parallel sampling) shift people’s processing strategies and their preferences between a short-term and long-term winning option. 
Our predictions were derived from the ecological and cognitive search hypothesis, which posits that the more complex the choice options, the more switch rates affect the costs of implementing different strategies for processing the samples. 
Our results did not support the core prediction that frequent switching and parallel sampling would increase reliance on a roundwise processing strategy that favors short-term winning choices.
Instead, option complexity shifted preferences under all sampling organizations similarly, with more short-term winner choices in complex than in simple problems. 
A computational latent classification analysis suggests that this complexity-related preference shift reflects an increase in guessing under high complexity, rather than a shift in processing strategy.

## How to reproduce the analyses?

1. To ensure full reproducibility, recreate the software setup using the following steps:

-   Install and select `R version 4.3.3` in `RStudio` under `Tools > Global Options > General`
-   Install the [`C++17 toolchain`](https://mc-stan.org/install/index.html) (required for `Stan`)
-   Download the repo via Fork & Clone or ZIP Download (`<> Code` button) and initialize it as R project
-   Run `renv::restore()` to install all required packages/dependencies with the same version

2. All results reported in the respective chapter of the dissertation can be reproduced using the provided [data/](https://github.com/linushof/sample-organization/tree/main/data) and [code/](https://github.com/linushof/sample-organization/tree/main/code) files. 

- Option 1: Start from the raw data and redo the data pre-processing and model fitting (takes some time). 

    - Retrieve raw data from [latest release](https://github.com/linushof/sample-organization/releases) and store subfolders in [data/raw/](https://github.com/linushof/sample-organization/tree/main/data/raw)
    - Source `preprocessing.R` to (over)write data sets in [data/clean/](https://github.com/linushof/sample-organization/tree/main/data/clean)
    - Source `analysis.R` to redo all analyses on clean data, write files in [fits/](https://github.com/linushof/sample-organization/tree/main/fits) and create all tables and figures

- Option 2: Start from the clean data and use provided model fits to access results more directly. 

    - Retrieve model fits from [latest release](https://github.com/linushof/sample-organization/releases) and store in [fits/](https://github.com/linushof/sample-organization/tree/main/fits)
    - Source `analysis.R` to inspect all analyses outputs and create all tables and figures


<!--Browse the [documentation](https://github.com/linushof/sample-organization/tree/main/documentation) for details about the files and variables in the data sets.-->