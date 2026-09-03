# Source data manifest

Every table below is an aggregate over studies, conditions, cohort phases or
populations, except the two marked DE-IDENTIFIED, which are per-sample value
vectors stripped of every column that could link a value to a sequencing run.
No table contains an enzyme-to-axis mapping. See the repository README for why.

| file | panel | rows | cols | sha256 (12) |
|---|---|--:|--:|---|
| `cmd/fig1B_3axis_plane.tsv` | Fig. 1B | 37 | 8 | `c70cdcbe29cc` |
| `cmd/cMD_dimensionality_scree.tsv` | Fig. 1B inset | 22 | 4 | `5f4fa316eff8` |
| `cmd/cMD_transdiagnostic_ordering_CORRECTED.tsv` | Fig. 3A | 6 | 4 | `dd3597e30034` |
| `cmd/cMD_disease_channel_matrix_V42_CORRECTED.tsv` | Fig. 3B | 24 | 8 | `c87e32fb0b51` |
| `cmd/channel_evidence_grade_cMD.tsv` | Fig. 3B bubble size | 306 | 7 | `73c74c52ae0e` |
| `neuro/neuro_v42_signatures.tsv` | Fig. 3C | 54 | 7 | `8571b2bc5cfa` |
| `engine/engine_gate_vs_continuous_cMD.tsv` | Fig. 2A | 13 | 9 | `ce5f9d55e546` |
| `engine/engine_gate_trajectory_anchors.tsv` | Fig. 2B | 30 | 7 | `e3d5b0053b9b` |
| `engine/entero_by_gate_cMD.tsv` | Fig. 2C | 4 | 6 | `c963b63bc507` |
| `anchors/fig2_paired_delta_feubioplus.tsv` | Fig. 4A | 12 | 14 | `f75260b5c02c` |
| `engine/fig4B_deltaEUBIO_decomposition_cohort_means.tsv` | Fig. 4B | 10 | 7 | `8e903d479466` |
| `anchors/fig3A_baryoseph_carbapenemase_paired_delta.tsv` | Fig. 4C | 13 | 5 | `1e249e9060c4` |
| `ici/fig4D_ici_adjusted_levels.tsv` | Fig. 4D | 104 | 3 | `44b36ff327b2` |
| `ici/ici_PRJEB70966_channel_response.tsv` | Fig. 4D annotation | 30 | 7 | `a0ac162c6933` |
| `ici/ici_GABA_composite_AUC.tsv` | Fig. 4D annotation | 4 | 6 | `8ecd229be4cc` |
| `populations/S10_Fig5_regenerated_EC200.tsv` | Fig. 5A, 5B | 6 | 9 | `81c896d47e7f` |
| `supp_fig1/suppfig1_calibration_invariance.tsv` | Supplementary Fig. 1 | 720 | 2 | `84c804b0fd2f` |
| `engine/tail_vs_mean_channel_condition_cMD.tsv` | Supplementary Fig. 2 | 175 | 14 | `fe1493827062` |
| `loso/TAB1_LOSO.tsv` | Supplementary Fig. 3 | 136 | 11 | `1c5eded51443` |

## What each file contains

**`cmd/fig1B_3axis_plane.tsv`** — Fig. 1B. Per-condition position on the three axes (covariate-adjusted betas), with the number of studies and cases. One row per cMD condition.

**`cmd/cMD_dimensionality_scree.tsv`** — Fig. 1B inset. Variance explained by each principal component of the axis matrix.

**`cmd/cMD_transdiagnostic_ordering_CORRECTED.tsv`** — Fig. 3A. Per-condition mean level of the damage and protective axes and their net difference.

**`cmd/cMD_disease_channel_matrix_V42_CORRECTED.tsv`** — Fig. 3B. Covariate-adjusted within-study Hedges' g per channel and condition, with channel valence. Same values as Supplementary Table 5.

**`cmd/channel_evidence_grade_cMD.tsv`** — Fig. 3B bubble size. Odds ratio, p, Benjamini-Hochberg q, number of studies and evidence grade per channel and condition. Internal bookkeeping column 'roster' removed.

**`neuro/neuro_v42_signatures.tsv`** — Fig. 3C. Covariate-adjusted channel effects in the Parkinson and Alzheimer cohorts.

**`engine/engine_gate_vs_continuous_cMD.tsv`** — Fig. 2A. Per-condition engine read continuously (odds ratio per SD) and as a threshold (odds ratio for the oxidised niche). Same values as Supplementary Table 3.

**`engine/engine_gate_trajectory_anchors.tsv`** — Fig. 2B. Per-cohort, per-phase mean and standard error of the sink-to-respiration log-ratio, with the fraction of samples in the oxidised niche.

**`engine/entero_by_gate_cMD.tsv`** — Fig. 2C. Enterobacteriaceae relative abundance in the four niche-by-group strata.

**`anchors/fig2_paired_delta_feubioplus.tsv`** — Fig. 4A. Per-cohort paired change in EUBIO under antibiotics or transplantation: n pairs, mean and median delta, fraction in the expected direction, Wilcoxon p.

**`engine/fig4B_deltaEUBIO_decomposition_cohort_means.tsv`** — Fig. 4B. Per-cohort mean decomposition of the within-subject change in EUBIO into an engine and a channel share. Derived from the internal per-subject table by keeping only the cohort means; see deidentify.py.

**`anchors/fig3A_baryoseph_carbapenemase_paired_delta.tsv`** — Fig. 4C. Paired pre/post carbapenemase abundance (EC 3.5.2.6, ppm) by decolonisation outcome, Bar-Yoseph cohort. Subject labels are the cohort's own, not sequencing accessions.

**`ici/fig4D_ici_adjusted_levels.tsv`** — Fig. 4D. Covariate-adjusted EUBIO and GABA-synthesis levels at baseline in the PRIMM checkpoint-immunotherapy cohort, by objective response. DE-IDENTIFIED: sample identifiers removed, rows shuffled, values rounded to 3 decimals; see deidentify.py.

**`ici/ici_PRJEB70966_channel_response.tsv`** — Fig. 4D annotation. Per-channel adjusted beta, p, Cliff's delta, AUC and q against objective response in the same cohort.

**`ici/ici_GABA_composite_AUC.tsv`** — Fig. 4D annotation. AUC and adjusted effect of the GABA-synthesis composite against two endpoints. Rows referencing an enzyme the paper does not disclose were removed; they are not used by the panel.

**`populations/S10_Fig5_regenerated_EC200.tsv`** — Fig. 5A, 5B. Per-population median EUBIO as measured and after adjustment for functional richness, with median enzyme count, Shannon diversity and a beneficial-taxa index. Same values as Supplementary Table 10.

**`supp_fig1/suppfig1_calibration_invariance.tsv`** — Supplementary Fig. 1. EUBIO under the frozen calibration against a cohort-internal standardisation, Parkinson cohort. DE-IDENTIFIED: sample identifiers and sequencing depth removed, rows shuffled, values rounded to 2 decimals; see deidentify.py.

**`engine/tail_vs_mean_channel_condition_cMD.tsv`** — Supplementary Fig. 2. Mean effect and threshold effect for every channel-by-condition association, with the detection class of each.

**`loso/TAB1_LOSO.tsv`** — Supplementary Fig. 3. Random-effects pooled Hedges' g per headline signal and the re-estimate with each contributing study removed in turn. Same values as Supplementary Table 6.

