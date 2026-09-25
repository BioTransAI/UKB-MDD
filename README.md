# UK Biobank Proteomics & MDD

**Analysis workflow accompanying the manuscript**

Pathway-based prediction of major depressive disorder (MDD), followed by protein interpretation, polygenic risk scoring, phenotype associations, genetic analyses, and druggability enrichment.

[Analysis script](analysis_workflow.R) · [Workflow](#workflow) · [Inputs](#inputs) · [Use](#use)

## Workflow

| Module | Analysis | Purpose |
| :---: | :--- | :--- |
| 1 | Model-fitting functions | Fit learners and generate predictions. |
| 2 | Feature-selection functions | Select pathway scores and optional unmapped proteins. |
| 3 | Model configuration comparison | Compare two-stage configurations by cross-validation. |
| 4 | Training and test evaluation | Evaluate the selected configuration in training and held-out test data. |
| 5 | Pathway coefficients | Rank pathway coefficients and annotate GO terms. |
| 6 | Protein coefficients | Estimate protein coefficients within retained pathways. |
| 7 | Polygenic risk scores | Calculate MDD scores using PRSice-2. |
| 8 | Phenotype associations | Relate predicted MDD scores to phenotypes, adjusting for covariates. |
| 9 | Forward Mendelian randomization | Evaluate the protein-to-MDD direction. |
| 10 | Reverse Mendelian randomization | Outline the MDD-to-protein direction; comments only. |
| 11 | Mediation analysis | Commented template for upstream protein, MDD status, and downstream protein. |
| 12 | Druggability enrichment | Assess druggability tiers and provide GREP analysis commands. |

The prediction workflow uses BioM2 with `glmnet` learners, 5-fold outer cross-validation, and 100-fold internal pathway-score reconstruction. Analysis parameters are specified in the script.

## Inputs

The script starts with prepared analysis objects; sample preparation and covariate preprocessing are outside its scope.

| Object | Expected structure |
| :--- | :--- |
| `train_data_final` | Training data frame with `label` first, followed by protein features. |
| `test_data_final` | Test data frame with the same outcome and feature structure. |
| `confounder_metadata` | Prepared covariates for downstream adjustment. |

Participant IDs are stored in row names. Outcome labels are `0` for controls and `1` for MDD cases. Additional inputs include pathway/protein annotations, phenotypes, selected protein sets, GWAS and pQTL statistics, genotype and LD resources, and druggability annotations.

Covariates are used in the phenotype analysis and the commented mediation template. Optional prediction-score adjustment is disabled by `adjust_confounders <- FALSE`.

## Use

1. Load the prepared objects and required software.
2. Replace descriptive `<PLACEHOLDERS>` with local input, output, or executable locations.
3. Run the implemented modules in sequence. Keep the reverse-MR outline commented and specify protein pairs before enabling the mediation template.
4. Execute GREP examples in a terminal; PRSice-2 commands are invoked from R.

<details>
<summary><strong>Software dependencies</strong></summary>

R packages referenced include `data.table`, `parallel`, `caret`, `mlr3verse`, `BioM2`, `glmnet`, `ModelMetrics`, `ROCR`, `mlr3measures`, `GO.db`, `AnnotationDbi`, `MASS`, `dplyr`, `TwoSampleMR`, `ggplot2`, `openxlsx`, `ieugwasr`, `genetics.binaRies`, `mediation`, `tidyr`, and `scales`. External tools include PLINK, PRSice-2, and GREP.

</details>

Model and analysis objects are saved as RDS files; some summary tables remain in R for inspection or export. This distribution contains code and documentation without participant-level data or computed results. The analysis was not rerun for these presentation edits.
