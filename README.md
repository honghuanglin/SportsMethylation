Duration of Contact Sports Play and Aberrant DNA Methylation in Human Frontal Cortex

This repository contains the analysis scripts used in the manuscript:

“Duration of contact sports play associated with aberrant DNA methylation in human frontal cortex.”

The scripts in this folder implement preprocessing, statistical modeling, and result generation for genome-wide DNA methylation analysis in postmortem frontal cortex tissue.


The repository includes R scripts for:
Processing smoothed methylation proportion data
Aligning samples with phenotype metadata
Running CpG-wise logistic regression models
Generating summary statistics for downstream analysis



The analyses evaluate whether DNA methylation levels across the genome are associated with:
CTE case status
Duration of contact sports participation
Other relevant covariates and neuropathologic markers (depending on script)
Methylation data consist of smoothed methylation proportions (per CpG, per sample) derived from whole-genome bisulfite sequencing.
Each chromosome-specific methylation matrix is processed separately and combined downstream.

This script performs the following:
Reads phenotype metadata
(Phenotype_CTE.txt)
Loads chromosome-level methylation data
(methy_proportion_smoothed.chr{chr}.txt)

Aligns subjects between phenotype and methylation matrices

Running the Script
Command-line usage
Rscript run_cte_methy_glm.R <chr>


Example (for chromosome 1):

Rscript Perform_association_test.R 1

Input directory structure
../Phenotype_CTE.txt
../methy_proportion_smoothed.chr1.txt
../methy_proportion_smoothed.chr2.txt
...

Output

A file named:
result.chr<chr>.txt

containing per-CpG regression summary statistics.

Dependencies

The scripts require:

R ≥ 4.0

CRAN packages:
data.table

You can install missing packages with:
install.packages("data.table")

Citation
If you use any scripts from this repository, please cite the associated manuscript:

"Duration of contact sports play associated with aberrant DNA methylation in human frontal cortex".

Contact
For questions about the code or manuscript, please contact the corresponding authors.
