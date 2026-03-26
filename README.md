# UKB-MDD Proteomics Analysis Pipeline

## Overview
Major Depressive Disorder (MDD) is a heterogeneous condition, yet its biological basis remains fragmented across molecular, physiological, and clinical domains. Here, leveraging plasma proteomics from 52,996 UK Biobank participants, we delineate a systems-level architecture of MDD using a biologically interpretable modeling framework. 
![BioM2](https://github.com/BioTransAI/UKB-MDD/edit/main/UKB_MDD.png)
## Key Features
* **Machine Learning & Identify Biomarkers:** Utilizes the BioM2 frameworks to identify highly predictive protein biomarkers and biological pathways
* **Phenotype Association:** Robust generalized linear models (GLM, LM, POLR) to correlate predicted risk scores with various clinical phenotypes.
* **Bidirectional Mendelian Randomization (MR):** Leverages TwoSampleMR to infer causal relationships between specific proteins and MDD.
* **Mediation Analysis:** Evaluates the indirect and direct effects of upstream/downstream proteins on MDD.
* **Druggable Genome Enrichment:** Overlaps identified proteins with the druggable genome (Tiers 1-3) using Fisher's Exact Tests.



## Prerequisites & Dependencies

### R Packages
Ensure you have the following R packages installed before running the pipeline:

```R
install.packages(c("data.table", "parallel", "caret", "MASS", "dplyr", "tidyr", "ggplot2", "scales", "openxlsx"))
```
Additionally, install specialized bioinformatics and ML packages:
* **mlr3verse:** `install.packages("mlr3verse")`
* **TwoSampleMR:** `remotes::install_github("MRCIEU/TwoSampleMR")`
* **ieugwasr:** `remotes::install_github("MRCIEU/ieugwasr")`
* **mediation:** `install.packages("mediation")`
* **BioM2:** `install.packages("BioM2") / devtools::install_github("BioTransAI/BioM2")`
* **GO.db:** `BiocManager::install("GO.db")`

### External Software
* **Python:** Required to run the GREP drug annotation script (`grep.py`).
