# Disease association template

Three scripts to test whether TAGMOS tiers predict your outcome(s):

| Script | Test | Output |
|---|---|---|
| `run_disease_OR.py` | T4-vs-T1 odds ratio + 95% CI + Cochran-Armitage trend + BH q | per-outcome OR table |
| `run_CA_trend.py` | Cochran-Armitage T1→T4 trend Z + p (standalone) | per-outcome trend |
| `run_cliff_delta.py` | Cliff δ on continuous Engine z (case vs control) | per-outcome δ |

The primary script is `run_disease_OR.py` — it computes the four
metrics that the TAGMOS paper reports uniformly throughout:
**CA Z, prevalence factor T4/T1, Fisher OR T4-vs-T1, BH q-value**.

---

## Input requirements

### TAGMOS classification TSV

The output of `classify_local_single.py` (Version A) or
`classify_local_composite.py` (Version B). Must contain at minimum:
`sample_id`, `ENG_tier`.

### Outcomes TSV

Per-sample outcomes with optional adjustment covariates:

```
sample_id    T2D    Celiac    IBD    age    sex    BMI
sample_001   0      0         0      45     F      24
sample_002   1      0         0      52     M      27
...
```

Outcome columns must be binary (0/1 or True/False). Continuous outcomes
need to be binarised first.

---

## Quick example

```bash
python run_disease_OR.py \
  --classification ../06_EXAMPLES/my_classification_demo.tsv \
  --outcomes ../06_EXAMPLES/synthetic_outcomes_n300.tsv \
  --out ../06_EXAMPLES/demo_disease_OR.tsv
```

Output format:

```
outcome    n_total  n_case  prev_T1_pct  prev_T2_pct  prev_T3_pct  prev_T4_pct  prev_factor  OR_T4_vs_T1  OR_95CI_lo  OR_95CI_hi  CA_trend_Z  CA_trend_p  BH_q
T2D        300      45      4.0          9.3          14.7         22.7         5.7×         5.43        2.14        13.78        +3.78       1.6e-4      4.8e-4
Celiac     300      18      2.7          4.0          6.7          10.7         4.0×         3.12        0.95         10.21        +2.21       2.7e-2      4.0e-2
...
```

## Interpreting the four metrics

Read together they give a **convergent four-way view** of any disease
association. See `07_FAQ.md` §"Interpreting CA Z / prev factor / OR / BH q"
for the full explanation. In short:

- **CA Z** — is there a monotonic gradient from T1 to T4?
- **prev factor** — how much more frequent is the outcome in T4 vs T1?
- **OR T4-vs-T1** — what's the odds-ratio strength of the T4-vs-T1 association?
- **BH q** — after correcting for testing N outcomes, is it still real?

A solid association shows convergence across all four. A finding that
shows CA Z significant but OR ~ 1 or BH q > 0.10 is suspicious.

---

## Beyond unadjusted OR

This template ships unadjusted OR + trend tests. If your cohort has
sufficient age / sex / BMI metadata you can extend to adjusted models. The
recommended template:

1. Run the unadjusted `run_disease_OR.py` first
2. Restrict to the outcomes that show OR_95CI excluding 1 AND BH_q < 0.10
3. For each of those outcomes, fit a logistic regression
   `outcome ~ T4_dummy + age + sex + BMI` using statsmodels
4. Report the adjusted OR alongside the unadjusted

A reference implementation is in the paper supplementary `SN5` (see
`HOW_TO_CITE.md` §"Recommended additional citations").

---

## Citation

Cite TAGMOS as required by `HOW_TO_CITE.md`. The disease-association template
itself does not need a separate citation — it implements the standard
statistics defined in the paper Methods §Statistics.

---

*Disease association template v10.5 · 2026-06-06*
