![alt text](https://github.com/matsov/TAGMOS/blob/main/TΛGM%CA%98S_Logo.jpg)


# TAGMOS - research · v4.5.11

Research-use-only distribution of the TAGMOS framework for individual-level
classification of gut microbiome ecosystems via a 24-axis architecture of
140 curated bottleneck enzyme commission entries.

> **License:** All components (code, schema, helpers) are released under the
> **TAGMOS Academic Research License v1.0** — academic research use ONLY.
> Commercial use of the Software, of its outputs ("Results"), or of any
> derivative work is prohibited and requires a separate written commercial
> license from the corresponding author (Wellmicro S.r.l.). See `LICENSE`.
> **No clinical or commercial use.** No calibration data are shipped — you
> must build your own calibration on your reference cohort.

## What's in the box

```
TAGMOS-research-v4.5.11-public-2D/
├── README.md                                       # this file
├── CALIBRATION.md                                  # how to build your own calibration
├── LICENSE                                         # TAGMOS Academic Research License v1.0
├── LICENSE-SCHEMA                                  # Schema-specific terms (academic use only)
├── CITATION.cff                                    # how to cite
├── requirements.txt                                # pip dependencies
├── TAGMOS/
│   ├── __init__.py
│   ├── classifier.py                               # core scoring (CLI: -m TAGMOS.classifier)
│   ├── calibrate.py                                # build calibration from controls
│   └── data/
│       ├── TAGMOS_schema_v4511_public_2D.json     # PUBLIC schema: 24 axes paper-aligned 2D
│       └── TAGMOS_calibration_TEMPLATE.json       # placeholder calibration (NOT clinical)
└── examples/
    ├── synthetic_ec.tsv                            # tiny demo input
    └── demo_run.sh                                 # end-to-end example
```

## Quickstart

### 1 · install Python dependencies

```bash
python3 -m pip install -r requirements.txt
```

### 2 · build your calibration

This software ships **without** calibration data. To produce meaningful
classifications you must construct a calibration JSON from a sample × EC
count table of *healthy reference controls* drawn from your own cohort
or institution. Minimum 50 controls (recommended ≥ 200).

```bash
python3 -m TAGMOS.calibrate \
    --controls my_healthy_controls_ec.tsv \
    --output   my_calibration.json \
    --calibration-id "my_cohort_v1"
```

Read **CALIBRATION.md** for the full protocol (sample-pipeline matching,
quality control, axis sanity checks).

### 3 · classify samples

```bash
python3 -m TAGMOS.classifier \
    --ec-tsv      my_samples_ec.tsv \
    --calibration my_calibration.json \
    --output      my_samples_TAGMOS.tsv
```

Output: per-sample TSV with 24 axis z-scores + 18 composite indices +
tier/class assignments. A `*.manifest.json` sidecar records SHA-256
hashes of the input, schema, calibration, and output for provenance.

## Demo

Run the end-to-end example on the synthetic dataset shipped under
`examples/`:

```bash
bash examples/demo_run.sh
```

This will (i) build a (toy) calibration from the synthetic controls,
(ii) classify the synthetic samples, (iii) print the first lines of the
classification output. It uses only the dummy template calibration, so
the resulting scores are **not** clinically interpretable — the demo is
purely for verifying that the software runs end-to-end on your system.

## Input format

Tab-separated values, first column is the sample id, remaining columns
are EC numbers (e.g. `1.1.1.1`, `2.3.1.85`, ...). Counts are non-negative
floats or integers; the classifier applies `log1p` then z-scoring internally.

```
sample_id    1.1.1.1   1.1.1.27   2.3.1.85   ...
SAMP_001     42        0          17         ...
SAMP_002     8         13         1          ...
```

Compatible with outputs from common functional metagenomics pipelines
(WMP, HUMAnN3 after EC aggregation, custom EC profilers).

## What this software does NOT contain

- **L0 calibration data** from the 10,000-subject Italian reference cohort
  (μ, σ, tertile cutoffs from the proprietary calibration). These are
  available from the corresponding author upon reasonable request,
  subject to a confidentiality agreement.
- **The WMP upstream profiling pipeline.** This software starts from a
  sample × EC count matrix and is upstream-agnostic — use any pipeline
  you prefer to produce the EC matrix.

## Citation

If you use TAGMOS-research in your work, please cite:

> Soverini M, et al. *Functional multi-axis decomposition of the human gut microbiome: reveals operational definition of eubiosis and dysbiosis* (2026) BiorXiv

See `CITATION.cff` for machine-readable metadata.

## License & contact

All components (code, schema, examples, helpers) are released under the
**TAGMOS Academic Research License v1.0** (see `LICENSE`).

This License authorizes use of the Software ONLY for non-commercial
academic research conducted within a recognized academic, public-research
or non-profit research institution. Any commercial use of the Software,
of its outputs ("Results"), or of any derivative work is expressly
prohibited and constitutes a material breach of the License. The
non-commercial restriction propagates to and binds every derivative
work, regardless of programming language or implementation medium.

Commercial use, clinical deployment, IVD deployment, internal R&D in a
for-profit organisation, monetization of Results, and any use that
generates revenue requires a separate written commercial license from
the Licensor — contact `andrea.castagnetti@wellmicro.com`.

For calibration access under confidentiality agreement, scientific
collaboration, or bug reports, contact
`andrea.castagnetti@wellmicro.com`.

---

© 2026 Wellmicro S.r.l. All rights reserved.
