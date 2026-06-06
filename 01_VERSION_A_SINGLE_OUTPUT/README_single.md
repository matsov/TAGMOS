# Version A · single-output 14 axes (paper 1 clean)

This is the classifier used in *Soverini et al. bioRxiv 2026*. It implements
the **14 single-output framework axes** of TAGMOS v4.5.13 FORMULA A FULL16
with **no composite indices**.

Use this version when you want to:

- replicate the primary paper's clinical-gradient analyses on your own cohort
- compare your cohort's substrate-functional partition against the paper's
  6,508-Italian RWE / 18,138-subject cMD external validation
- extend or modify the base framework (new axis, new EC, new formula)
- publish results that are directly comparable to the primary TAGMOS paper

For composite indices (MUC_INFLAM, BAI_MCAS, PATHOBIONT_triple,
EUB_complete_resilience, TMA_net), use **Version B** at
`../02_VERSION_B_WITH_COMPOSITE/`.

---

## Files

| File | Purpose |
|---|---|
| `classifier_single.py` | core module · AXIS_EC_DICT + AXIS_FORMULA + compute_axis_* functions |
| `calibrate_local_single.py` | CLI · compute calibration JSON on YOUR cohort (one time) |
| `classify_local_single.py` | CLI · apply calibration to a new EC-count matrix |
| `README_single.md` | this file |

---

## Three-step workflow

### Step 1 · provide an EC-count TSV

EC ids as first column, samples as remaining columns.
Counts as non-negative integers (or floats).

```
EC_id        sample_A   sample_B   sample_C
1.2.7.4      12         0          7
2.3.1.169    8          3          21
2.8.4.1      0          0          15
...
```

### Step 2 · calibrate

```bash
python calibrate_local_single.py \
    --ec-counts your_cohort_EC_counts.tsv \
    --cohort-tag MyStudyName \
    --out my_calibration.json
```

This runs ONCE per cohort and produces `my_calibration.json` with all the
per-EC / per-sub-axis / per-axis mu / sigma and quartile / tertile cuts.

Minimum recommended n = 200 samples; see
`../03_CALIBRATION_RECIPE/minimum_cohort_requirements.md`.

### Step 3 · classify

```bash
python classify_local_single.py \
    --ec-counts your_cohort_EC_counts.tsv \
    --calibration my_calibration.json \
    --out my_classification.tsv
```

Output columns:

- `sample_id`
- For each of the 14 framework axes: `<AXIS>_z` (the axis z-score)
- For Engine z: `ENG_tier` (T1_EUBIOTIC / T2_PRESERVED / T3_ALTERED / T4_DYSBIOTIC)
- For SYN_cof: `SYN_cof_class` (SYN_FRAGILE / SYN_PARTIAL / SYN_RESILIENT)
- For SYN_carb: `SYN_carb_class` (SYN_carb_FRAGILE / SYN_carb_PARTIAL / SYN_carb_RESILIENT)
- For each channel: `<CHANNEL>_class` (<CHANNEL>_LOW / _MID / _HIGH)
- `cell_3D`: 3D state-space cell label
  (e.g. `T1_EUBIOTIC × SYN_RESILIENT × SYN_carb_PARTIAL`)

---

## Determinism + reproducibility

- Every output is deterministic given (EC count matrix, calibration JSON).
- The calibration JSON is SHA-256 hash-locked; the hash is in the JSON itself.
- Two analysts running the same EC matrix through the same calibration JSON
  will get bit-identical output.

---

## Citation

You must cite TAGMOS in any publication. See `../HOW_TO_CITE.md`.

The recommended naming for your calibration instance:
`TAGMOS-<your-cohort-tag>-<n>` (e.g. `TAGMOS-StanfordIBD-850`).

Reference the Wellmicro calibration of the primary paper as
`TAGMOS-WMP-IT-6508`.

---

*Version A README v10.5 · 2026-06-06*
