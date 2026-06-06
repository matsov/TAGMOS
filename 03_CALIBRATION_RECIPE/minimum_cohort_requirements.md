# Minimum cohort requirements for TAGMOS calibration

This document specifies what your cohort must look like to support a reliable
TAGMOS calibration.

---

## 1. Sample size

| Cohort size | Status | Implications |
|---|---|---|
| n ≥ 500 | **recommended** | Robust quartile / tertile cuts, stable mu / sigma, suitable for disease-association analyses with adjusted models. |
| 200 ≤ n < 500 | **acceptable** | Cuts and z-scores are computable; recommend permutation null and bootstrap on disease-association results to confirm stability. |
| 100 ≤ n < 200 | **discouraged** | Cuts unstable; sigma estimates noisy; bootstrap mandatory. Override with `--min-samples 100`. |
| n < 100 | **refused** | Calibration is not statistically meaningful. Use the public Wellmicro calibration only as an external reference for individual sample classification (not implementable in this bundle by design). |

Default refusal threshold: `--min-samples 200`.

---

## 2. EC coverage

Your upstream functional profiler must cover at least 60 % of the 151
bottleneck EC entries of the framework. The calibration script reports the
per-sub-axis coverage; sub-axes with coverage = 0 will produce sigma = 0
and therefore z = 0 for every sample on that sub-axis (effectively
neutralising it).

**Critical EC subsets** — if any of the following are missing, results are
materially affected:

- Engine z eubiotic sinks: at least 1 EC each from H_sink_acetogenesis_WL,
  H_sink_methanogenesis, H_sink_sulfate_reduction. Missing all three of one
  category gives Engine z a 1-pole-only definition.
- Engine z dysbiotic sinks: at least 2 of EC 7.1.1.9 / 7.1.1.7 / 1.7.1.15.
  Missing these collapses the dysbiotic pole.
- SYN_cof: at least 1 EC from each of B12, K2, POL, GABA.
- SYN_carb: BUT sub-axis (EC 2.8.3.8) — if missing, falls back to TAX panel
  if available, otherwise SYN_carb is effectively PRP + LAC only.

The calibrator reports `EC coverage per axis` in stdout. Inspect this before
trusting the output.

---

## 3. Sample homogeneity

The TAGMOS calibration assumes that all samples in the calibration cohort
were processed through the **same upstream profiling pipeline** (same
adapter trimming, same host-read depletion, same EC quantification thresholds).
If your cohort mixes samples from multiple pipelines (e.g. some HUMAnN3 +
some custom DIAMOND), the per-EC mu / sigma will be pipeline-confounded.

Recommended: calibrate one TAGMOS instance per pipeline, then compare
instances at the partition level (cell-3D occupancy) rather than at the
z-score level.

---

## 4. Demographic structure (for downstream disease-association analyses)

The calibration step itself is **demographic-blind** — the EC matrix and
the cohort size are all that enter. However, if you intend to use the
calibration for disease-association analyses with adjusted models, your
cohort should have:

- balanced sex distribution (40-60 %); pure single-sex cohorts can be
  calibrated but lose cross-sex generalisability
- age range covering the relevant clinical age range
- minimum metadata per sample: sample_id, sex, age, BMI, outcome label(s)

If demographic metadata is incomplete, restrict downstream analyses to
unadjusted T4-vs-T1 odds ratio + Cochran-Armitage trend Z + Cliff δ
(`run_disease_OR.py`, `run_CA_trend.py`, `run_cliff_delta.py` in
`04_DISEASE_ASSOCIATION_TEMPLATE/`).

---

## 5. Outcome / disease labels

For disease-association analyses (Step 4 of QUICK_START.md):

- Binary outcomes: 0 / 1 or "control" / "case" (the script will normalise)
- Ordinal outcomes: ordinal integer encoding (e.g. healthy / mild / severe → 0 / 1 / 2)
- Minimum prevalence per outcome: at least 10 cases of each class in T1
  AND in T4 (so the T4-vs-T1 cell of the 2×2 OR table is populated). If
  any cell is < 5, the script falls back to Fisher exact + Haldane–Anscombe
  correction with a warning.

---

## 6. Practical sanity checks

After calibration, before publishing, please run these:

- `verify_calibration.py --calibration my_calibration.json` (in
  `03_CALIBRATION_RECIPE/`) checks internal consistency.
- Inspect `ENG_tier` distribution after classification: it should be
  approximately 25 / 25 / 25 / 25 on the calibration cohort (the cuts are
  quartiles by construction).
- Inspect SYN_cof / SYN_carb / channel class distributions: they should be
  approximately 33 / 33 / 33.
- Inspect 36-cell occupancy: typically 30+ of 36 cells are populated on a
  cohort of n ≥ 500; cohorts of n < 200 may show 10-25 populated cells.

---

*Minimum cohort requirements v10.5 · 2026-06-06*
