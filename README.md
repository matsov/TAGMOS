![alt text](https://github.com/matsov/TAGMOS/blob/main/TΛGM%CA%98S_Logo.jpg)


# TAGMOS — public research bundle v10.5

**TAGMOS** is a substrate-functional multi-axis framework for stratifying the
human gut microbiome. It reduces a metagenome to a small set of single-output
axes calibrated on a reference cohort, partitions samples into four ordinal
tiers (T1_EUBIOTIC → T4_DYSBIOTIC) and produces an operational definition of
eubiosis and dysbiosis that has been externally validated across 18,138
metagenomes (92 cMD studies) and tested for cross-population transferability
on > 590 traditional / paleofecal samples.

This bundle is the **public release for academic, non-commercial research use**.

> **Read first**: `LICENSE.md` (non-commercial-research-only · mandatory TAGMOS citation) and `QUICK_START.md` (15-minute walk-through).

---

## What's inside

```
TAGMOS_PUBLIC_BUNDLE_v10_5/
├── README.md                      ← this file
├── LICENSE.md                     ← Public Research Licence v1.0
├── QUICK_START.md                 ← 15-minute walk-through
├── HOW_TO_CITE.md                 ← citation template + naming convention
├── 00_FRAMEWORK_SPEC/             ← framework definition (no calibration numerics)
│   ├── schema_tagmos_public_v4513.json
│   ├── axes_registry.md
│   ├── ec_dictionaries.md
│   └── architecture_overview.md
├── 01_VERSION_A_SINGLE_OUTPUT/    ← *paper 1 clean* · 14 single-output axes only
│   ├── README_single.md
│   ├── classifier_single.py
│   ├── calibrate_local_single.py
│   └── classify_local_single.py
├── 02_VERSION_B
├── 03_CALIBRATION_RECIPE/         ← step-by-step end-to-end calibration on YOUR cohort
│   ├── CALIBRATION_PROTOCOL.md
│   ├── minimum_cohort_requirements.md
│   └── verify_calibration.py
├── 04_DISEASE_ASSOCIATION_TEMPLATE/
│   ├── README_disease_template.md
│   ├── run_disease_OR.py          ← T4-vs-T1 odds ratio + 95% CI + q_BH
│   ├── run_CA_trend.py            ← Cochran-Armitage T1→T4 trend Z
│   └── run_cliff_delta.py         ← continuous Cliff δ per axis
├── 05_VALIDATION_TOOLKIT/
│   ├── beta_diversity_falsification.py
│   ├── dimensional_saturation.py
│   └── input_method_routing.py
├── 06_EXAMPLES/
│   ├── synthetic_demo_cohort_n300.tsv     ← synthetic EC count table
│   ├── synthetic_outcomes_n300.tsv        ← synthetic disease outcomes
│   └── walkthrough_results.md
└── 07_FAQ.md
```

---


### Version A · single-output axes only (paper 1 clean)

`01_VERSION_A_SINGLE_OUTPUT/`

- 14 framework single-output axes: **1 PRIMARY** (Engine z) + **2 CO-PRIMARY**
  aggregator axes (SYN_cof, SYN_carb) + **11 independent ancillary channels**
  (BA, TRP, PROT, MUC, HIS, IRON, LPS, ETU, TMA, UREM, HYS).
- One z-score per axis per sample — no composite metrics.
- This is the classifier used in *Soverini et al. bioRxiv 2026*.
- Suitable for: monotonic clinical-gradient testing, cross-cohort replication,
  reductionist disease-association analyses, methodological extensions of the
  base framework.

### Version B · with composite indices (research layer)
 
 Under development

---

## What you need to provide

To use TAGMOS on your own cohort you only need an EC-level functional
count matrix, in the format:

```
EC_id           sample_A   sample_B   sample_C   ...
EC_1.2.7.4      12         0          7
EC_2.3.1.169    8          3          21
EC_2.8.4.1      0          0          15
...
```

Either Level-4 EC numbers from shotgun functional profilers (WMP, HUMAnN3,
Meteor2) or any equivalent EC-level quantification.

Plus, for disease-association analyses, a binary or ordinal outcome table:

```
sample_id    outcome_label    age    sex    BMI
sample_A     control          45     F      24
sample_B     case             52     M      27
...
```

---

## What this bundle does NOT contain (and why)

The numeric values of the Wellmicro Italian RWE calibration anchor —
per-axis mean, standard deviation and quartile / tertile cuts — are NOT
released in this bundle. They are held as proprietary calibration parameters
of Wellmicro S.r.l.

**This does not prevent you from using TAGMOS.** The framework is designed so
that any researcher can **calibrate on their own cohort** following the
recipe in `03_CALIBRATION_RECIPE/`. The classifier scripts operate on the
user's own calibration parameters, never on the Wellmicro frozen anchors.

See `03_CALIBRATION_RECIPE/CALIBRATION_PROTOCOL.md` for the step-by-step
procedure, and `07_FAQ.md` §3 for the rationale.

---

## Dependencies

Python ≥ 3.9 with:

```bash
pip install numpy pandas scipy statsmodels
```

Optional (for visualisation):

```bash
pip install matplotlib seaborn
```

No external network access required at runtime. No proprietary database lookups.

---

## Quick start (15 minutes)

```bash
# 1. activate your Python env
# 2. recalibrate on your cohort (one-time):
cd 03_CALIBRATION_RECIPE
python ../01_VERSION_A_SINGLE_OUTPUT/calibrate_local_single.py \
  --ec-counts your_cohort_EC_counts.tsv \
  --out my_calibration.json

# 3. classify your samples:
python ../01_VERSION_A_SINGLE_OUTPUT/classify_local_single.py \
  --ec-counts your_cohort_EC_counts.tsv \
  --calibration my_calibration.json \
  --out my_classification.tsv

# 4. run disease association on the result:
cd ../04_DISEASE_ASSOCIATION_TEMPLATE
python run_disease_OR.py \
  --classification ../my_classification.tsv \
  --outcomes your_outcomes.tsv \
  --out disease_OR_results.tsv
```

End-to-end on the synthetic 300-sample demo cohort takes < 1 minute on a laptop.

---

## Citation requirement

Use of this bundle for any publication, preprint, presentation, dataset or
software release requires citation of the primary TAGMOS publication in the
form specified in `HOW_TO_CITE.md` and `LICENSE.md` §4.

Failure to cite TAGMOS as specified constitutes a breach of the Public
Research Licence and terminates your right to use the Software.

---

## Help / contact

For methodological questions or to report a bug:

- File an issue on the public GitHub repository at: `https://github.com/<TBD>/TAGMOS`
- Or contact the corresponding author: `andrea.castagnetti@wellmicro.com`

For commercial licensing, see `LICENSE.md` §3 and §10.

---

*TAGMOS v10.5 public research bundle - Wellmicro®*
