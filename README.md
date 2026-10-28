# Latent class analysis code and item coding

Materials for reproducing the latent class analysis reported in:

> Exploratory Symptom Profiles and Clinical Correlates in Adults with Attention-Deficit/Hyperactivity Disorder: A Retrospective Study Using the Wender-Reimherr Adult Attention Deficit Disorder Scale

## Contents

- `lca/`: the master script and helper scripts used for the latent class analysis and its model-fit and stability checks.
- `post_processing/alternative_groupings_S9.R`: calculates the alternative grouping proportions, posterior-weighted estimate, and five-class domain means from existing fitted models; it does not re-estimate any model.
- `item_domain_mapping.csv`: the Excel and SPSS item variable names, WRAADDS domain labels, response codes, and coding note.

No participant-level data, fitted model objects, logs, or analysis outputs are included.

## Software

The latent class analyses were conducted using **R 4.3.3** and **poLCA 1.6.0.1**. The scripts retain version checks. The master script also requires the `readxl` package.

The S9 post-processing script was run using **R 4.6.1** and requires `readxl`. It reads the existing `fit_k5.rds` and `fit_k6.rds` objects produced by the LCA workflow; it does not re-estimate any model.

## Input data and execution

The clinical input workbook is not included because it contains sensitive participant-level information. The workflow expects an authorized Excel workbook with a sheet named `IN_EX_Final`, 749 complete records, the 28 item variables listed in `item_domain_mapping.csv`, and a `Class_neu` column containing the six-class assignment used in the manuscript. The workflow uses that column in model-fit and stability checks. Response codes 1, 2, and 3 correspond to the original WRAADDS ratings 0–2 shifted by +1.

To run the analysis, use R 4.3.3 with `poLCA` 1.6.0.1 and `readxl`. In RStudio, open and source `lca/run_final_lca.R`, then select the authorized Excel workbook when prompted. The workflow performs extensive multi-start model searches and bootstrap procedures and may require many hours.

Output files are written next to the selected workbook in `lca_final_R433_poLCA1601`. They include fitted models and logs; these can contain participant-level item responses, posterior probabilities, or class assignments and should not be shared publicly without appropriate authorization.

The S9 script writes aggregate CSV files and session information to a local output folder; these generated files are not included in this repository.

## Reuse

No license is specified for this code; all rights are reserved by the authors. For reuse, please contact the corresponding author, Nicole Mauche, Department of Psychiatry and Psychotherapy, University of Leipzig Medical Center.
