# FAQ

## 1. Why don't you ship the Wellmicro Italian calibration numerics?

The numeric calibration parameters of the Wellmicro 6,508-subject Italian RWE
anchor are held as proprietary calibration parameters of Wellmicro S.r.l.
The framework was however designed so that **any researcher can calibrate
locally** on their own cohort. The classifier scripts operate exclusively on
your own calibration; they never read the Wellmicro frozen anchor.

This stance is consistent with the way most production-grade frameworks
release: schema + classifier are public, the production calibration of the
issuing company is held. The scientific value is in the framework + the
calibration recipe, both of which are fully released here.

## 2. Can I publish results computed with this bundle?

Yes, in academic non-commercial research venues. Cite TAGMOS as required by
`HOW_TO_CITE.md`. Failure to cite constitutes a breach of the Public Research
Licence (`LICENSE.md` §4).

## 3. Can I use this bundle in a commercial product or service?

No. The Public Research Licence explicitly prohibits commercial use. Contact
`andrea.castagnetti@wellmicro.com` for a separate commercial licence.

## 4. What's the difference between Version A and Version B?

- **Version A** ships the 14 single-output framework axes used in the
  primary TAGMOS paper. One z-score per sample per axis. No composites.
- **Version B** adds research-layer **composite indices** (MUC_INFLAM,
  BAI_MCAS_score, PATHOBIONT_triple, PATHOBIONT_continuous, TMA_net,
  EUB_complete_resilience) defined in `02_VERSION_B_WITH_COMPOSITE/composite_definitions.md`.

If you're doing methodologically-rigorous work where each axis must have a
clear mechanistic interpretation, use Version A. If you're doing exploratory
sub-phenotype discovery, Version B may give better discrimination but at
some interpretability cost.

The single-output axes are bit-identical between Version A and Version B.

## 5. How do I cite a recalibrated instance?

Use the naming convention `TAGMOS-<cohort-tag>-<n>` (see `HOW_TO_CITE.md`
§4). Example: `TAGMOS-StanfordIBD-850`.

The original Wellmicro calibration is referenced as `TAGMOS-WMP-IT-6508`.

## 6. How do I interpret the four uniform metrics (CA Z / prev factor / OR / BH q)?

The primary paper reports the same four metrics for every disease association.
Each tells you a different thing about the same association:

| Metric | Question it answers |
|---|---|
| **CA Z** (Cochran-Armitage trend) | Is there a monotonic gradient from T1 to T4? |
| **prev factor** (T4 prevalence / T1 prevalence) | Is the outcome more common in T4 than in T1, and by what factor? |
| **OR T4-vs-T1** (Fisher exact odds ratio with 95 % CI) | What's the strength of the T4-vs-T1 association? |
| **BH q-value** | After correcting for the multiple outcomes I tested, is this still significant? |

A solid association shows **convergence** across all four:
- CA Z strong and positive (typically > +3)
- prev factor > 2
- OR with 95 % CI excluding 1
- BH q < 0.10 (paper-defined threshold)

A finding showing CA Z significant but OR ~ 1 or BH q > 0.10 is suspicious —
either the effect size is small (worth flagging but not headline-worthy) or
the trend is driven by a single tier rather than a true gradient.

## 7. My cohort is small (n < 200). Can I still use this?

The default minimum is `--min-samples 200`. Below that, quartile / tertile
cuts become noisy. You can override (`--min-samples 100`) but expect:
- noisy partitions (cell-3D occupancy < 15 cells of 36)
- wider OR confidence intervals (often crossing 1 even for true associations)
- BH q corrections may eliminate all signals

Recommended workflow for small cohorts: bootstrap the calibration and
disease-association statistics 100-200× to obtain CI on every estimate.

## 8. My upstream pipeline isn't HUMAnN3 or WMP. Can I still use this?

Yes, **as long as your pipeline produces Level-4 EC counts**. The bundle is
profiler-agnostic. We have validated it on WMP and HUMAnN3 in the primary
paper; users have applied it to Meteor2 and custom DIAMOND BLASTX pipelines
without issue.

The only pipeline-specific caveat is **EC coverage**: pipelines with restrictive
gene-family granularity may miss some bottleneck ECs. Run
`05_VALIDATION_TOOLKIT/input_method_routing.py` after calibration to check
per-sub-axis EC coverage on your cohort. If coverage < 60 % on multiple
sub-axes, results should be reported with explicit pipeline-conditional
caveat.

## 9. What if my cohort and the published paper disagree on a result?

This is expected and informative. Disease-association effect sizes depend on:
(a) cohort architecture (RWE vs case-control, age structure, comorbidity)
(b) sequencing depth and pipeline
(c) sample size
(d) outcome definition (self-report vs ICD vs clinical interview)

If your effect size differs from the published paper by < 2× and the direction
is preserved, results are consistent. If direction is opposite, investigate
(typically: pipeline-conditional EC coverage gap, or outcome-definition
mismatch).

The paper's robustness battery (Supp Note 5) addresses sources of variation
systematically. Refer to it for the expected range of variation.

## 10. Can I propose changes to the framework?

Yes — the framework is open to extension via derivative work under the
Public Research Licence. Common extensions:

- New EC entries (e.g. for novel bottleneck enzymes)
- New channel axes (e.g. for substrates not yet covered)
- Reweighting of composite indices for new clinical contexts

Document your modification clearly in your paper's Methods and reference
TAGMOS as the base framework. See `HOW_TO_CITE.md` §5 for the exact wording.

For changes you'd like to see merged into the official TAGMOS schema,
contact `andrea.castagnetti@wellmicro.com` with your rationale + data.

## 11. Where's the official online repo?

The official TAGMOS public bundle is hosted at:
- GitHub: `https://github.com/<TBD>/TAGMOS`
- Zenodo (versioned snapshots): TBD

Both are linked from the primary paper's Data and code availability section.

## 12. How do I report a bug?

File an issue on the public GitHub repository with:
- the exact command you ran
- the full stdout / error
- the bundle version (`v10.5`)
- (if relevant) a minimal reproducible EC count subset

For security issues, contact `andrea.castagnetti@wellmicro.com` directly.

## 13. Why is "MUC_INFLAM" the name of the composite and not a single axis?

In Version A the framework axis is plain `MUC` — it aggregates
MUC_specialised + MUC_general. In Version B the **composite** index
`MUC_INFLAM` is different: it z-means MUC_specialised with LPS and HIS,
providing a multi-axis inflammatory signature.

`MUC` (axis) and `MUC_INFLAM` (composite) are NOT the same metric. Use
`MUC` in mechanistic disease-association work; use `MUC_INFLAM` in
inflammatory phenotype discrimination.

## 14. The Engine z tier in my cohort is 50 % T1 and 30 % T4. Is that normal?

The calibration step always produces approximately 25 / 25 / 25 / 25 on
the **calibration cohort itself** (by construction — the cuts are quartiles
of the cohort). If you classify a different cohort using the calibration of
the first one, the distribution can be skewed.

If you classified your cohort using ITS OWN calibration and got a skewed
distribution, something is wrong: check that you used the same EC matrix
for both calibration and classification, and that no sample-id mismatch
caused the calibration to be computed on a different sub-cohort than the
one you classified.

## 15. The clinical OR I'm getting doesn't replicate the paper. Should I worry?

Read FAQ #9 carefully. Most discrepancies are explained by cohort
architecture differences. If after due diligence (matched outcome definition,
matched pipeline, matched sample size, matched age structure) you still
cannot replicate, contact the corresponding author with details. This is
how the field improves.

---

*FAQ v10.5 · 2026-06-06*
