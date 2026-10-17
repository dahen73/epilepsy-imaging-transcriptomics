# Pathology-specific transcriptomic coupling of cortical morphology in drug-resistant focal epilepsy

Analysis code accompanying the manuscript submitted to *Epilepsia* (2026).

The pipeline combines standardised FreeSurfer-derived cortical metrics from the **IDEAS dataset**
(OpenNeuro [ds005602](https://openneuro.org/datasets/ds005602), CC0; 442 drug-resistant focal epilepsy
patients who underwent resective surgery, 100 controls) with region-level gene expression from the
**Allen Human Brain Atlas** (via the [BrainEnrich](https://github.com/ZZH-Finalize/BrainEnrich) R package),
and tests whether inter-individual differences in cortical thickness, surface area, volume and sulcal
depth are coupled to 86 pre-specified epilepsy-relevant gene sets, at cohort level and within
pathology strata (hippocampal sclerosis vs non-HS as primary strata).

## Requirements

- **R** (≥ 4.3; analyses run on R 4.3.2 / 4.5.3) with packages: `BrainEnrich`, `data.table`, `Matrix`,
  `ggplot2`, `patchwork`
- **Python** (≥ 3.10) with `nibabel` (sulcal-depth extraction only)
- IDEAS data release (OpenNeuro ds005602) including the raw FreeSurfer archives for sulcal depth
- AHBA expression digest as shipped with BrainEnrich (`get_geneExp(atlas = "desikan", ...)`)

## Pipeline (scripts run in numeric order)

### Stage 1 — Input preparation and validation
| Script | Purpose |
|---|---|
| `01_prepare_inputs.R` | Build analysis inputs from the IDEAS release: AHBA left-hemisphere expression (34 Desikan regions × 6,513 genes), individual/group thickness matrices, phenotype table |
| `02_smoke_test.R` | 10-case end-to-end smoke test (DisGeNET 4,622 sets) |
| `03_name_hits.R`, `align_check.R`, `09_anno_check.R`, `id_check.R`, `na_check.R` | Region-name alignment, annotation and missing-data checks |

### Stage 2 — Full-library screens
| Script | Purpose |
|---|---|
| `04_full_analysis.R` | Group-level enrichment (13 configurations) + individual-level full-library screens (parametric and spin-permutation, 5,000 permutations) |
| `05_b2_summary.R` – `08_fcd_check.R`, `b1_debug.R`, `corr_diff.R`, `t_diff.R`, `t_diff2.R`, `lmfit_test.R` | Stage summaries and checks |

### Stage 3 — Targeted 86-set panel (main FDR clusters)
| Script | Purpose |
|---|---|
| `10_targeted.R` | 86-set panel (8 mTOR / 13 GABAergic / 65 epilepsy-DisGeNET) × 4 modalities × 6 predictors, individual level, spin 5,000; FCD II (n = 37) extracted from operative memos |
| `11_tg_debug.R` – `14_tg_final.R`, `verify_repro.R` | Targeted-pipeline debugging/finalisation and reproducibility check against the reference implementation (agreement < 1e-6) |

### Stage 4 — Modality extension (surface area / volume / sulcal depth)
| Script | Purpose |
|---|---|
| `15_prepare_area_vol.R`, `16_targeted_area_vol.R` | Area and volume modalities into the pipeline |
| `extract_sulc.sh`, `chain_sulc_prep.sh`, `17_sulc_extract.py`, `18_prepare_sulc.R`, `19_targeted_sulc.R` | Vertex-wise sulcal depth from 542 FreeSurfer tar.gz archives, then into the pipeline |
| `20_cross_modal_summary.R`, `21_region_summary.R` | Cross-modal summary; 34-region age-at-onset gradient |

### Stage 5 — Robustness and sensitivity analyses (S1–S6)
| Script | Purpose |
|---|---|
| `24_robustness.R` | S1 resection-burden covariation; S2 leave-one-out jackknife; S3 80% bootstrap ×200 |
| `25_individual_resectexcl.R` | S4 exclusion of 19 high-resection regions with full spin re-testing |
| `26_robustness_summary.R` | S1–S4 summary + robustness figure |
| `29_prepare_rh.R`, `30_individual_rh.R`, `31_rh_vs_lh.R`, `32_make_fig6.R` | S5 right-hemisphere sensitivity (5,000 right-hemisphere spin permutations) |
| `34_asm_mechanism.R`, `36_make_fig7.R` | S6 full severity adjustment (Freedman–Lane), dose–response, gene-set-family specificity |

### Stage 6 — Pathology-stratified analyses and outputs
| Script | Purpose |
|---|---|
| `41_stratified_pathology.R` | Primary HS vs non-HS strata (within-stratum BH-FDR) and exploratory small strata (DNT, cavernoma, OTHER, dual pathology); stratum-vs-cohort beta concordance |
| `22_make_figures.R`, `gen_strat_outputs.R` | Main figures (Figures 1–3) |
| `37_make_tables.R`, `gen_suptable1.R` | Main and supplementary tables |

Helper scripts used during development (`install_brainenrich*.R`, `install_v3.R`, `install_v4.R`)
install the BrainEnrich package and are included for completeness; they are not part of the analysis
pipeline.

## Data availability

The IDEAS dataset is publicly available at OpenNeuro (ds005602, CC0). Gene-expression data derive
from the Allen Human Brain Atlas via the BrainEnrich package. No new human data were generated.

## Licence

We suggest releasing this code under the MIT licence — please add a `LICENSE` file with the
copyright holders' names before uploading.

## Citation

If you use this code, please cite the manuscript (citation to be added on publication) and the
IDEAS dataset paper: Taylor PN et al., *Epilepsia* 2025;66(2):471–481, doi:10.1111/epi.18192.
