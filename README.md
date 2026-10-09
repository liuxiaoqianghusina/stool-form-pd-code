# Stool form and incident Parkinson's disease in UK Biobank

Analysis code for the manuscript *"Opposite associations of hard and loose stools with incident Parkinson's disease: a prospective phenome-wide study in UK Biobank"* (submitted).

The study relates stool form (frequency of hard and of loose stools) and bowel movement frequency, reported on the UK Biobank Digestive Health Questionnaire (2017–2020), to incident Parkinson's disease and to 395 other diseases. It adds external data from NHANES 2005–2010, two-sample Mendelian randomization, colocalisation and brain iron imaging (QSM).

## Data

No individual-level data are included in this repository. UK Biobank data are available to approved researchers through https://www.ukbiobank.ac.uk. NHANES 2005–2010 files and the public-use linked mortality files are available from the US Centers for Disease Control and Prevention. GWAS summary statistics were obtained from the GWAS Catalog (GCST90565484, stool frequency; GCST009325, Parkinson's disease) and from FinnGen release 12.

## Setup

Paths are written as placeholders. Before running, replace `/path/to/project/` with the project folder (which should contain `scripts/`, `data/`, `results/`, `figures/` and `manuscript/`) and `/path/to/ukb_data/` with the folder holding the UK Biobank extract files.

Analyses were run with R 4.2.2 (survival 3.4.0, survey 4.4.2, mice 3.16.0, TwoSampleMR 0.7.11, MRPRESSO 1.0, data.table 1.15.4, splines) and Python 3.14 (pandas 2.3.3, numpy 2.4.6, scipy 1.17.1, statsmodels 0.14.6, matplotlib, python-docx, openpyxl, Pillow).

## Pipeline

Scripts are numbered in the order they were run. `prep.R` is sourced by the R survival scripts and defines the analytic sample, exposures, covariates and models.

| Script | Purpose |
|---|---|
| `01_build_dataset.py` | Build the cohort (exposures, covariates, exclusions, censoring) and the outcome table from the UK Biobank extract |
| `02_phewas_cox.R` | Phenome-wide Cox models for 395 outcomes |
| `03_key_outcomes.R` | Key outcomes: model ladder, sensitivity analyses, follow-up windows, Fine–Gray models, E-values |
| `04_dose_response.R` | Dose–response analyses |
| `05_figures.py` | Phenome-wide figure and frequency forest plot |
| `06_nhanes_build.py`, `07_nhanes_analysis.R`, `07b_nhanes_pd_meds.R` | NHANES 2005–2010 data preparation and survey-weighted analyses |
| `08_mr.R`, `11_reverse_mr.R` | Forward and reverse two-sample Mendelian randomization |
| `09_qsm_nigra.py` | Bowel habits and QSM/T2* of subcortical nuclei |
| `10_coloc.py` | Colocalisation (approximate Bayes factors, coloc.abf algorithm) |
| `13_pd_mutual_adjustment.R` | Joint models of stool form and frequency; Schoenfeld tests |
| `14_flowchart.py` | Study flow figure |
| `16_time_origin_comparison.R` | Follow-up from questionnaire completion versus recruitment |
| `17_prediction.R` | Incremental predictive value (C index, cross-validation, net reclassification, enrichment) |
| `18a_extract_extra_fields.py`, `18_tier1.R` | Negative-control exposures, alternative PD definitions, transit spectrum score, subgroups, absolute risk |
| `18b_extract_accel_meds.py`, `19a_med_flags.py`, `19_tier2.R` | Joint analysis with other prodromal markers, accelerometry, medication exclusions, multiple imputation |
| `20_figures_pd.py`, `fig_tiff.py`, `21_collect_figures.sh` | Main and supplementary figures (PDF, PNG and 600-dpi TIFF) |
| `22_loose_stool_causes.R` | Sensitivity analysis excluding causes of loose stools |
| `12_tables.py`, `23_manuscript_tables.py` | Main and supplementary tables, Supplementary Information and Supplementary Data files |
| `15_assemble_manuscript.py`, `24_renumber_supplement_text.py`, `md2docx.py` | Manuscript assembly utilities (reference list from PubMed, supplementary numbering, Markdown to Word) |

## License

MIT License (see `LICENSE`).
