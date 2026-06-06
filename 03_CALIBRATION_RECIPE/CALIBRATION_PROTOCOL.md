# Calibration protocol · how to compute YOUR TAGMOS calibration

This is the **central operational protocol** of the public bundle. It describes
how to compute a TAGMOS calibration on your own cohort, what to check, and
how to name and publish the resulting calibration instance.

The recipe is identical for Version A (single-output) and Version B
(with composite indices) — only the script names change.

---

## Why you calibrate locally

The numeric calibration parameters of the Wellmicro Italian RWE 6,508-subject
anchor (per-axis mean, standard deviation, quartile / tertile cuts) are
proprietary calibration parameters of Wellmicro S.r.l. and are NOT shipped
with this bundle.

This is by design. The framework was built so that any researcher can
**calibrate on their own cohort**, producing a per-cohort instance of TAGMOS
that classifies their samples according to *their* cohort's internal
distribution.

Calibration on your own cohort is not a workaround for the missing
Wellmicro numerics — it is the **canonical use mode** of the public bundle.
A cohort-internal calibration ensures that the tier partition (T1–T4) and
the syntrophy / channel class partitions (FRAGILE / PARTIAL / RESILIENT,
LOW / MID / HIGH) are computed against the natural variation of your own
samples, which is what makes the resulting clinical associations
meaningful.

---

## Five-step protocol

### Step 1 · prepare the input

You need an EC-count TSV with EC ids in the first column and samples in the
remaining columns. EC ids are Level-4 numeric strings (e.g. "1.2.7.4").

```
EC_id        sample_001   sample_002   sample_003   ...
1.2.7.4      12           0            7
2.3.1.169    8            3            21
2.8.4.1      0            0            15
...
```

The matrix can come from any upstream profiler that produces Level-4 EC
counts (HUMAnN3 / bioBakery, Meteor2, custom DIAMOND BLASTX pipelines, etc.).

### Step 2 · check the minimum requirements

See `minimum_cohort_requirements.md`. In short:

- n ≥ 200 samples recommended (refused by default below this)
- balanced sex / age structure if you intend to use adjusted models
- minimum EC coverage 60 % of the framework dictionary (151 entries)
- if BUT sub-axis EC coverage is low, expect the BUT z to be unreliable;
  consider providing a complementary species-abundance table for the
  TAX-first BUT fallback (documented in `axes_registry.md` §SYN_carb)

### Step 3 · calibrate

```bash
# Version A · single output
python ../01_VERSION_A_SINGLE_OUTPUT/calibrate_local_single.py \
    --ec-counts your_cohort_EC_counts.tsv \
    --cohort-tag MyStudyName \
    --out my_calibration.json
```

Or, for Version B with composite indices:

```bash
python ../02_VERSION_B_WITH_COMPOSITE/calibrate_local_composite.py \
    --ec-counts your_cohort_EC_counts.tsv \
    --cohort-tag MyStudyName \
    --out my_calibration_composite.json
```

`--cohort-tag` should be a short identifier (max 20 characters, alphanumeric
+ hyphen) that describes the calibration cohort. This becomes part of the
`calibration_id`. Examples:

- `--cohort-tag StanfordIBD` → produces `TAGMOS-StanfordIBD-<n>`
- `--cohort-tag EurMS` → produces `TAGMOS-EurMS-<n>`
- `--cohort-tag NishijimaJP` → produces `TAGMOS-NishijimaJP-<n>`

The output JSON contains:

- `calibration_id` (e.g. `TAGMOS-StanfordIBD-850`)
- `calibration_sha256` (SHA-256 hash-lock of the JSON contents)
- per-sub-axis and per-axis mu / sigma
- Engine z quartile cuts (Q25 / Q50 / Q75)
- tertile cuts for SYN_cof, SYN_carb, and each channel (Q33 / Q66)
- the EC dictionary used (for reproducibility)
- the bundle version + the calibration date

### Step 4 · verify

Run the verification script to confirm the calibration is internally
consistent:

```bash
python verify_calibration.py --calibration my_calibration.json
```

The verifier checks that:

- the SHA-256 hash matches the recomputed hash of the calibration content;
- the Engine z cuts are monotonic (Q25 < Q50 < Q75);
- each tertile cut is monotonic (Q33 < Q66);
- each sub-axis has finite mu and positive sigma;
- the framework axes registry is complete;
- no L0 / Wellmicro fields are present in the JSON (the bundle never
  introduces them).

### Step 5 · publish

Once calibration is verified, treat the calibration JSON as a versioned
research artifact. Recommended publication:

- Deposit the calibration JSON in a public repository (e.g. Zenodo) with a
  DOI;
- Reference the DOI + the calibration `calibration_id` and `calibration_sha256`
  in your Methods;
- Cite TAGMOS as specified in `HOW_TO_CITE.md`.

Example Methods sentence:

> *"We recalibrated TAGMOS (Soverini et al., bioRxiv 2026) on the present
> cohort (n = 850; IBD case-control sub-cohort) following the public
> calibration recipe (TAGMOS_PUBLIC_BUNDLE_v10_5/03_CALIBRATION_RECIPE/) under
> the TAGMOS Public Research Licence v1.0. The resulting calibration instance,
> hereafter referred to as TAGMOS-StanfordIBD-850, is deposited at Zenodo
> doi:[…] with SHA-256 hash a1b2c3…"*

---

## Optional · cross-calibration alignment with the published Italian anchor

If you want to verify that your cohort's substrate-functional partition is
consistent with the published Italian-RWE anchor (the *Soverini et al. 2026*
calibration `TAGMOS-WMP-IT-6508`), you can:

1. Compute the **distribution shape** statistics on your Engine z (skewness,
   kurtosis, IQR / SD ratio) and compare them to those reported in the paper's
   Supplementary Note 5.
2. Run the **β-diversity falsification** toolkit (in
   `05_VALIDATION_TOOLKIT/`) on your cohort. The published Italian benchmark
   is **38× to 94× variance advantage** of the 36-cell partition vs any single
   clinical outcome. Your cohort should produce a similar order-of-magnitude
   ratio if the framework is performing consistently.
3. If you have an outcome that overlaps with the published paper (e.g. T2D,
   Celiac, IBD), compare your T4-vs-T1 odds ratio with the Italian paper's.
   Order-of-magnitude agreement (factor 2–3) is expected; exact numerical
   agreement is not, because cohort architecture and effect size differ.

None of these checks require the Wellmicro numerics. They are
distribution-level and effect-level sanity checks.

---

## Frequently encountered issues

### Insufficient EC coverage

If your upstream profiler covers < 60 % of the framework EC dictionary, the
calibration is still computed but some sub-axes will be weak or near-zero.
Inspect the calibration JSON `per_sub_axis_mu_sigma` for sub-axes with
sigma ≈ 0 — these were ill-defined on your cohort.

### Very small cohorts (n < 200)

The default `--min-samples 200` prevents calibration on cohorts that are too
small to estimate quartile cuts reliably. You can lower this with explicit
awareness (`--min-samples 100`) but expect noisy partitions and weak
disease-association statistics.

### TAX-first BUT sub-axis

The Version B production schema uses a TAX-first species-marker panel for the
BUT sub-axis. The public bundle implements an EC-fallback (EC 2.8.3.8) which
may be less accurate on cohorts where the butyryl-CoA pathway is taxonomy-
specific. If your cohort has species-abundance data, see
`02_VERSION_B_WITH_COMPOSITE/composite_definitions.md` for guidance on
overlaying the TAX panel.

---

*Calibration protocol v10.5 · 2026-06-06*
