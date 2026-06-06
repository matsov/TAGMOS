# Walkthrough results · what to expect on the demo cohort

This is the expected output of the QUICK_START walk-through on the synthetic
`synthetic_demo_cohort_n300.tsv` cohort (300 samples, 83 EC entries from the
framework dictionary).

The synthetic cohort was constructed so that the **T2D-like outcome shows a
strong monotonic gradient** from T1 to T4 (designed effect), **Celiac-like a
moderate one**, and **IBS_null no association** (negative control).

---

## Step 1 · calibrate (Version A)

```
[TAGMOS calibrate] n samples:                  300
[TAGMOS calibrate] n EC entries detected:      83
[TAGMOS calibrate] writing calibration JSON to: 06_EXAMPLES/my_calibration_demo.json
[TAGMOS calibrate] calibration_id: TAGMOS-DemoSynth-300
[TAGMOS calibrate] sha256:         fbe361f3000643b7...
```

## Step 2 · verify

```
✓ SHA-256 hash matches
✓ Engine z quartile cuts monotonic
✓ tertile cuts monotonic for 13 axes
✓ all 36 sub-axes have finite mu and positive sigma
✓ all 14 framework axes present in per_axis_mu_sigma
✓ no forbidden L0 / Wellmicro keys present
✅ PASS
```

## Step 3 · classify

```
Engine tier distribution (% of samples):
  T1_EUBIOTIC     25.0 %
  T2_PRESERVED    25.0 %
  T3_ALTERED      25.0 %
  T4_DYSBIOTIC    25.0 %

36-cell occupancy: 28 populated cells of 36
```

The 25 / 25 / 25 / 25 distribution is **by construction** — the cuts are
the quartiles of the calibration cohort. When you classify a **different**
cohort using THIS calibration, you'll see a non-balanced distribution.

## Step 4 · disease OR

```
 outcome  n_total  n_case  prev_T1  prev_T2  prev_T3  prev_T4  prev_factor  OR    OR_95CI         CA_Z   BH_q
 T2D      300      64      9.3 %    18.7 %   25.3 %   32.0 %   3.43×        4.57  [1.83, 11.4]    +3.53  0.001
 Celiac   300      36      6.7 %    12.0 %   13.3 %   16.0 %   2.4 ×        2.67  [0.89,  7.99]   +1.75  0.121
 IBS_null 300      38     14.7 %    9.3 %    5.3 %    21.3 %   1.45×        1.58  [0.68,  3.67]   +0.93  0.352
```

Reading the table:

- **T2D shows the convergent four-way signature** of a real association:
  CA Z > +3, prev factor > 3, OR with 95 % CI excluding 1, BH q < 0.01.
- **Celiac is marginal** — direction is right but the OR 95 % CI crosses 1
  and BH q > 0.10. With n = 300 and 12 % prevalence we don't reach
  conclusive significance.
- **IBS_null is correctly null** — no monotonic gradient, OR ~ 1, BH q ~ 0.35.

This is a clean demonstration that:
1. The framework correctly identifies the designed association (T2D)
2. The framework correctly does NOT over-call a random outcome (IBS_null)
3. Small effect sizes need larger cohorts to detect (Celiac)

## Step 5 · validation toolkit · input-method routing

```
Per-sub-axis EC coverage and routing recommendation:
  most sub-axes:                 100 % EC coverage → use EC (rule R9 pass)
  H_sink_acetogenesis_WL:        100 %
  H_sink_methanogenesis:         100 %
  H_sink_sulfate_reduction:      100 %
  ...
```

On a real cohort with HUMAnN3 output you may see lower coverage in
specialised mucin sub-axes (MUC_specialised: 60-80 %) and selected B12
sub-axes (B12_DMB_ligand: 40-60 %). Inspect `demo_R9_routing.tsv` for the
per-sub-axis numbers and decide whether to flag any caveats in your
downstream interpretation.

---

## Differences when using Version B (with composite indices)

```
PATHOBIONT_triple positive:    9 / 300 samples
EUB_complete_resilience pos:  19 / 300 samples
```

The composite indices add columns to the classification TSV
(MUC_INFLAM, BAI_MCAS_score, PATHOBIONT_triple, PATHOBIONT_continuous,
TMA_net, EUB_complete_resilience). The single-output axes are bit-identical
between Version A and Version B.

---

## What this demo shows you

After this walk-through you should be confident that:

1. The bundle scripts execute correctly on a typical input.
2. The framework produces the expected tier distribution on the calibration
   cohort.
3. The disease-association template returns the four-way uniform metrics
   ready for publication tables.
4. The calibration JSON is L0-clean (no Wellmicro proprietary numerics).
5. Both Version A and Version B work end-to-end on the same EC matrix.

You are ready to swap the synthetic cohort for your real EC count table and
run the same workflow on your data.

---

*Walkthrough v10.5 · 2026-06-06*
