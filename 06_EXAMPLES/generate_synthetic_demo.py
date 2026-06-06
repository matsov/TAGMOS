"""
Generate a synthetic demo cohort for walking through the TAGMOS bundle.

Produces:
  - synthetic_demo_cohort_n300.tsv  : EC count matrix · 151 EC × 300 samples
  - synthetic_outcomes_n300.tsv     : sample_id + binary outcomes + age/sex/BMI

The synthetic data is designed to *exhibit* TAGMOS-typical behaviour:
  - A T2D-like binary outcome with elevated prevalence in the dysbiotic tail
  - A Celiac-like binary outcome with monotonic gradient
  - A "null" outcome (random, no association) for negative control

It is NOT intended as a benchmark — only as a tutorial cohort to verify
that the bundle scripts run end-to-end on a typical input.
"""

import sys, math
from pathlib import Path
import numpy as np
import pandas as pd

OUT_DIR = Path(__file__).resolve().parent

# A representative subset of the 151 bottleneck EC entries — keep the same
# IDs as in classifier_single.AXIS_EC_DICT so the example actually exercises
# all framework axes.
EC_LIST = [
    # Engine WL
    "1.2.7.4", "2.3.1.169", "6.3.4.3", "1.5.1.5", "3.5.4.9", "1.5.1.20",
    # Engine methano
    "2.8.4.1", "1.5.98.1", "1.5.98.2", "1.5.98.3", "1.8.7.3",
    # Engine sulfate
    "1.8.99.5", "1.8.4.8",
    # Engine dysbiotic
    "7.1.1.9", "7.1.1.7", "1.7.1.15", "1.7.2.2", "1.7.5.1", "1.7.2.1",
    "1.7.1.4", "1.9.3.1",
    # B12
    "2.1.1.151", "4.99.1.3", "2.5.1.17", "6.6.1.2", "2.4.2.21",
    "5.4.99.2", "2.1.1.13",
    # K2
    "4.1.1.71", "6.2.1.26", "4.1.2.62", "4.1.3.36", "3.1.2.28", "2.1.1.163",
    # Polyamine + GABA
    "4.1.1.17", "4.1.1.19", "3.5.3.11", "4.1.1.15",
    # SYN_carb
    "2.8.3.8", "4.2.1.28", "1.1.1.202", "5.1.2.1",
    # BA
    "3.5.1.24", "1.3.1.115", "1.17.99.1", "1.1.1.391",
    "1.1.1.50", "1.1.1.52", "6.2.1.7", "2.3.1.65",
    # TRP
    "4.1.99.1", "1.4.3.2", "4.1.1.74", "1.13.11.11", "1.13.11.52", "3.5.1.9",
    # PROT
    "4.1.1.25", "4.1.1.18",
    # MUC
    "3.2.1.49", "3.2.1.97", "3.2.1.96", "3.2.1.51", "3.2.1.18",
    "3.2.1.22", "3.2.1.4",
    # HIS
    "4.1.1.22", "4.3.1.3", "2.3.1.5",
    # TMA
    "4.3.99.4", "4.3.99.3", "1.7.99.1", "4.3.99.1",
    "1.14.13.148", "2.1.1.157",
    # UREM, IRON, LPS, ETU, HYS
    "4.1.99.5", "1.13.11.27",
    "6.3.2.39", "6.3.2.14",
    "2.3.1.129", "2.3.1.241",
    "4.3.1.7",
    "4.2.1.22", "2.5.1.48",
]


def main():
    rng = np.random.default_rng(seed=42)
    n_samples = 300
    n_ec = len(EC_LIST)

    sample_ids = [f"sample_{i:03d}" for i in range(1, n_samples + 1)]

    # Latent "dysbiosis score" per sample, used to modulate eubiotic vs
    # dysbiotic EC means. Range roughly [-1.5, +1.5] with normal shape.
    dysbiosis = rng.normal(0, 1, size=n_samples)

    # Build EC matrix:
    # eubiotic ECs (Engine WL/methano/sulfate + SYN aggregator inputs) are
    # higher in low-dysbiosis samples; dysbiotic ECs (Engine dysbiotic + LPS)
    # are higher in high-dysbiosis samples. Channel ECs are independent.
    eubiotic_ecs = set(EC_LIST[:13] + EC_LIST[21:38] + EC_LIST[42:50])
    dysbiotic_ecs = set(EC_LIST[13:21] + ["2.3.1.129", "2.3.1.241"])

    matrix = np.zeros((n_ec, n_samples), dtype=float)
    for i, ec in enumerate(EC_LIST):
        if ec in eubiotic_ecs:
            mu = 25 - 8 * dysbiosis
        elif ec in dysbiotic_ecs:
            mu = 25 + 10 * dysbiosis
        else:
            mu = 25 + rng.normal(0, 4, size=n_samples)
        # negative-binomial-like spread
        scale = 1.0 + 0.2 * rng.normal(size=n_samples)
        vals = np.maximum(0, mu * scale + rng.normal(0, 6, size=n_samples)).astype(int)
        matrix[i] = vals

    ec_df = pd.DataFrame(matrix, index=EC_LIST, columns=sample_ids)
    ec_df.index.name = "EC_id"
    ec_path = OUT_DIR / "synthetic_demo_cohort_n300.tsv"
    ec_df.to_csv(ec_path, sep="\t")
    print(f"[demo] EC matrix written:  {ec_path}  ({n_ec} EC × {n_samples} samples)")

    # Outcomes:
    # T2D-like: probability ~ logistic(0.6 * dysbiosis - 1.5)
    # Celiac-like: probability ~ logistic(0.45 * dysbiosis - 2.0)
    # IBS_null: random Bernoulli 10 %
    def _logistic(x):
        return 1 / (1 + np.exp(-x))

    p_t2d = _logistic(0.65 * dysbiosis - 1.4)
    p_celiac = _logistic(0.45 * dysbiosis - 2.0)
    t2d = rng.binomial(1, p_t2d)
    celiac = rng.binomial(1, p_celiac)
    ibs_null = rng.binomial(1, 0.10, size=n_samples)

    # demographics
    age = rng.integers(20, 80, size=n_samples)
    sex = rng.choice(["F", "M"], size=n_samples)
    bmi = np.round(rng.normal(25, 3, size=n_samples), 1)

    out_df = pd.DataFrame({
        "sample_id": sample_ids,
        "T2D": t2d,
        "Celiac": celiac,
        "IBS_null": ibs_null,
        "age": age,
        "sex": sex,
        "BMI": bmi,
    })
    out_path = OUT_DIR / "synthetic_outcomes_n300.tsv"
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"[demo] outcomes written:    {out_path}")
    print(f"[demo] expected prevalences:")
    print(f"        T2D       {100*t2d.mean():.1f} %")
    print(f"        Celiac    {100*celiac.mean():.1f} %")
    print(f"        IBS_null  {100*ibs_null.mean():.1f} %  (random)")


if __name__ == "__main__":
    main()
