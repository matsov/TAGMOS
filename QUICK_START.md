# Quick start · 15 minutes

A minimal end-to-end walk-through on the bundled synthetic demo cohort
(`06_EXAMPLES/synthetic_demo_cohort_n300.tsv`, 300 samples, 151 bottleneck EC
entries, synthetic disease labels).

## Step 0 · check your Python env

```bash
python3 -c "import numpy, pandas, scipy, statsmodels; print('OK')"
```

If anything is missing:

```bash
pip install numpy pandas scipy statsmodels
```

## Step 1 · pick a version

You have **two parallel branches**:

- **Version A** · paper 1 clean · 14 single-output axes
- **Version B** · research layer · 14 single-output axes + composite indices

Read `README.md` §"Two versions" to decide. For this walk-through we use
**Version A**.

## Step 2 · calibrate on your cohort

This step is **done once** per cohort. It computes per-EC mean / sigma and
per-axis mean / sigma / quartile / tertile cuts on YOUR samples, and saves
them into a calibration JSON. The classifier will use this JSON, not any
Wellmicro proprietary numerics.

```bash
cd 01_VERSION_A_SINGLE_OUTPUT/
python calibrate_local_single.py \
  --ec-counts ../06_EXAMPLES/synthetic_demo_cohort_n300.tsv \
  --out ../06_EXAMPLES/my_calibration_demo.json
```

You should see something like:

```
[TAGMOS calibrate] reading EC counts:          ../06_EXAMPLES/synthetic_demo_cohort_n300.tsv
[TAGMOS calibrate] n samples:                  300
[TAGMOS calibrate] n EC entries detected:      151
[TAGMOS calibrate] computing per-EC log1p mu/sigma...
[TAGMOS calibrate] aggregating axis-level z scores (14 axes)...
[TAGMOS calibrate] computing quartile cuts on Engine z + tertile cuts on syntrophy and channels...
[TAGMOS calibrate] writing calibration JSON to: ../06_EXAMPLES/my_calibration_demo.json
[TAGMOS calibrate] done.
```

Cohort minimum requirements: see `03_CALIBRATION_RECIPE/minimum_cohort_requirements.md`
(suggested n ≥ 200, balanced sex / age structure if you intend to use the
adjusted models).

## Step 3 · classify

Now apply the calibration you just generated to each sample:

```bash
python classify_local_single.py \
  --ec-counts ../06_EXAMPLES/synthetic_demo_cohort_n300.tsv \
  --calibration ../06_EXAMPLES/my_calibration_demo.json \
  --out ../06_EXAMPLES/my_classification_demo.tsv
```

The output is a TSV with one row per sample and, for each sample, one column
per axis (z-score) plus tier / class labels:

```
sample_id   ENG_z   ENG_tier   SYN_cof_z   SYN_cof_class   SYN_carb_z   SYN_carb_class   BA_z   TRP_z   ...
sample_001  +0.42   T1         +0.18       SYN_PARTIAL     -0.31        SYN_FRAGILE      +1.02  -0.50   ...
sample_002  -1.13   T4         -0.84       SYN_FRAGILE     +0.65        SYN_RESILIENT    -0.41  -0.13   ...
...
```

## Step 4 · disease association

Now ask: are TAGMOS tiers associated with my outcome?

```bash
cd ../04_DISEASE_ASSOCIATION_TEMPLATE/
python run_disease_OR.py \
  --classification ../06_EXAMPLES/my_classification_demo.tsv \
  --outcomes ../06_EXAMPLES/synthetic_outcomes_n300.tsv \
  --out ../06_EXAMPLES/demo_disease_OR.tsv
```

You get a table with:

- Per-outcome prevalence in each tier (T1, T2, T3, T4)
- T4-vs-T1 odds ratio with 95 % CI
- Cochran-Armitage trend Z statistic
- Benjamini–Hochberg q-value across all tested outcomes
- Cliff δ between case and control on continuous Engine z

See `04_DISEASE_ASSOCIATION_TEMPLATE/README_disease_template.md` for the
full output format and how to interpret it.

## Step 5 · (optional) reproduce the β-diversity falsification

If you want to verify on your own cohort that the TAGMOS partition explains
substantially more compositional variance than diversity-based partitions:

```bash
cd ../05_VALIDATION_TOOLKIT/
python beta_diversity_falsification.py \
  --ec-counts ../06_EXAMPLES/synthetic_demo_cohort_n300.tsv \
  --classification ../06_EXAMPLES/my_classification_demo.tsv \
  --outcomes ../06_EXAMPLES/synthetic_outcomes_n300.tsv \
  --out ../06_EXAMPLES/demo_falsification.tsv
```

On the published 6,508-Italian RWE cohort this analysis returns 38× to 94×
more variance explained by the 36-cell partition than by any individual
clinical outcome on the same distance matrix.

## Step 6 · cite TAGMOS

Before publishing any result from the bundle, read `HOW_TO_CITE.md` and
`LICENSE.md` §4. Citation is mandatory under the Public Research Licence.

---

## What's next

- For a deeper walk-through of all 14 axes and their biological interpretation,
  read `00_FRAMEWORK_SPEC/axes_registry.md`.
- For Version B (with composite indices) the workflow is identical — replace
  `01_VERSION_A_SINGLE_OUTPUT/` with `02_VERSION_B_WITH_COMPOSITE/`.
- For frequently asked questions: `07_FAQ.md`.

---

*Quick start v10.5 · 2026-06-06*
