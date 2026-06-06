"""
TAGMOS calibrate · Version A · single-output 14 axes only
=========================================================

Compute the cohort-local calibration JSON from an EC-count table. This is the
"set-up step": it is run ONCE per cohort and produces all the per-sub-axis
mu / sigma, per-axis mu / sigma, quartile / tertile cuts that the classifier
will subsequently apply to every sample of the cohort.

The output JSON is SHA-256 hash-locked. It contains:
    - calibration_id  (with cohort tag + n + SHA-256)
    - schema version  (TAGMOS public v4513 FORMULA_A_FULL16)
    - per_sub_axis_mu_sigma   (one entry per sub-axis)
    - per_axis_mu_sigma       (one entry per framework axis)
    - engine_quartile_cuts    (Q25 / Q50 / Q75 of Engine z)
    - tertile_cuts            (Q33 / Q66 of each tertile-partitioned axis)
    - bundle_version
    - calibration_date_utc

Usage
-----
    python calibrate_local_single.py \\
        --ec-counts your_cohort_EC_counts.tsv \\
        --cohort-tag StanfordIBD \\
        --out my_calibration.json

Required input format
---------------------
EC count TSV with EC ids as the first column ("EC_id" or "ec") and samples
as the remaining columns. EC ids must be Level-4 numeric strings, e.g.
"1.2.7.4". Counts must be non-negative integers (or floats; log1p is applied
upstream so non-integers are tolerated).

Example::

    EC_id        sample_001  sample_002  sample_003  ...
    1.2.7.4      12          0           7
    2.3.1.169    8           3           21
    2.8.4.1      0           0           15
    ...

Output
------
A self-contained JSON file at the path given by --out. The file is the
"calibration instance" of TAGMOS on YOUR cohort. Name it consistently with
the convention `TAGMOS-<cohort-tag>-<n>` (see HOW_TO_CITE.md §4) when
referencing it in publications.

License
-------
TAGMOS Public Research License v1.0 — see ../LICENSE.md. Mandatory citation.
"""

import argparse, json, hashlib, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

# allow `python calibrate_local_single.py` from within the version dir
sys.path.insert(0, str(Path(__file__).resolve().parent))
from classifier_single import (
    AXIS_EC_DICT, AXIS_FORMULA, FRAMEWORK_AXES_14,
    compute_sub_axis_log1p, compute_axis_raw,
)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--ec-counts", required=True,
                    help="TSV with EC ids as first column, samples as columns")
    ap.add_argument("--cohort-tag", default="UserCohort",
                    help="Short identifier of the calibrated cohort, "
                         "alphanumeric + hyphen (max 20 chars). "
                         "Becomes part of the calibration_id.")
    ap.add_argument("--out", required=True,
                    help="Output calibration JSON path")
    ap.add_argument("--min-samples", type=int, default=200,
                    help="Refuse to calibrate on fewer than this many samples "
                         "(default 200; see minimum_cohort_requirements.md)")
    args = ap.parse_args()

    print(f"[TAGMOS calibrate] reading EC counts:          {args.ec_counts}")
    df = pd.read_csv(args.ec_counts, sep="\t")
    if df.columns[0] not in ("EC_id", "ec", "EC", "ec_id"):
        # fall back to first column as index
        ec_col = df.columns[0]
        print(f"[TAGMOS calibrate] using first column as EC id: '{ec_col}'")
    else:
        ec_col = df.columns[0]
    df = df.set_index(ec_col)
    df.index = df.index.astype(str)
    print(f"[TAGMOS calibrate] n samples:                  {df.shape[1]}")
    print(f"[TAGMOS calibrate] n EC entries detected:      {df.shape[0]}")

    if df.shape[1] < args.min_samples:
        sys.exit(f"[TAGMOS calibrate] ERROR: only {df.shape[1]} samples (< "
                 f"min_samples={args.min_samples}). Calibration on small cohorts "
                 f"is unreliable; either provide more samples or lower "
                 f"--min-samples explicitly with awareness of the implications "
                 f"(see minimum_cohort_requirements.md).")

    # ---- sub-axis level scores per sample
    print(f"[TAGMOS calibrate] computing per-sub-axis log1p sums...")
    sub_axis_matrix = compute_sub_axis_log1p(df)

    # ---- per-sub-axis mu / sigma on the calibration cohort
    sub_axis_mu_sigma = {}
    for sub_axis in sub_axis_matrix.index:
        vals = sub_axis_matrix.loc[sub_axis].dropna().values
        mu = float(np.mean(vals))
        sigma = float(np.std(vals, ddof=1)) if len(vals) > 1 else 1.0
        if sigma <= 0:
            sigma = 1.0
        sub_axis_mu_sigma[sub_axis] = (mu, sigma)

    # ---- compute axis raw scores from sub-axes
    print(f"[TAGMOS calibrate] aggregating axis-level scores ({len(FRAMEWORK_AXES_14)} framework axes)...")
    axis_raw = compute_axis_raw(sub_axis_matrix, sub_axis_mu_sigma)

    # ---- per-axis mu / sigma on the calibration cohort
    axis_mu_sigma = {}
    for axis in axis_raw.index:
        vals = axis_raw.loc[axis].dropna().values
        mu = float(np.mean(vals))
        sigma = float(np.std(vals, ddof=1)) if len(vals) > 1 else 1.0
        if sigma <= 0:
            sigma = 1.0
        axis_mu_sigma[axis] = (mu, sigma)

    # ---- z-score axis_raw using axis_mu_sigma to get final axis z
    axis_z = {}
    for axis in axis_raw.index:
        mu, sigma = axis_mu_sigma[axis]
        axis_z[axis] = ((axis_raw.loc[axis] - mu) / sigma).values

    # ---- compute quartile cuts on Engine z + tertile cuts on the others
    print(f"[TAGMOS calibrate] computing quartile cuts on Engine z + tertile cuts on syntrophy + channels...")
    engine_q = np.percentile(axis_z["ENG"], [25, 50, 75]).tolist()
    tertile_cuts = {}
    for axis in FRAMEWORK_AXES_14:
        if axis == "ENG":
            continue
        cuts = np.percentile(axis_z[axis], [33.333, 66.667]).tolist()
        tertile_cuts[axis] = cuts

    # ---- assemble calibration JSON
    cal = {
        "calibration_id": f"TAGMOS-{args.cohort_tag}-{df.shape[1]}",
        "schema_version": "v4513.public.1.0",
        "formula_layer":  "FORMULA_A_FULL16",
        "bundle_version": "v10.5.public",
        "n_samples_calibration": int(df.shape[1]),
        "cohort_tag": args.cohort_tag,
        "calibration_date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "per_sub_axis_mu_sigma": {k: list(v) for k, v in sub_axis_mu_sigma.items()},
        "per_axis_mu_sigma": {k: list(v) for k, v in axis_mu_sigma.items()},
        "engine_quartile_cuts": {
            "Q25": engine_q[0],
            "Q50": engine_q[1],
            "Q75": engine_q[2],
        },
        "tertile_cuts": {
            axis: {"Q33": cuts[0], "Q66": cuts[1]}
            for axis, cuts in tertile_cuts.items()
        },
        "tier_labels": {
            "ENG":      ["T4_DYSBIOTIC", "T3_ALTERED", "T2_PRESERVED", "T1_EUBIOTIC"],
            "SYN_cof":  ["SYN_FRAGILE", "SYN_PARTIAL", "SYN_RESILIENT"],
            "SYN_carb": ["SYN_carb_FRAGILE", "SYN_carb_PARTIAL", "SYN_carb_RESILIENT"],
        },
        "ec_axis_dictionary_used": AXIS_EC_DICT,
        "license_note": "TAGMOS Public Research License v1.0 · non-commercial research · mandatory citation",
        "primary_citation": "Soverini M., et al. bioRxiv 2026, doi:[TBD]",
    }

    payload_for_hash = json.dumps(cal, sort_keys=True).encode("utf-8")
    cal["calibration_sha256"] = hashlib.sha256(payload_for_hash).hexdigest()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(cal, f, indent=2)

    print(f"[TAGMOS calibrate] writing calibration JSON to: {out_path}")
    print(f"[TAGMOS calibrate] calibration_id: {cal['calibration_id']}")
    print(f"[TAGMOS calibrate] sha256:        {cal['calibration_sha256'][:16]}...")
    print(f"[TAGMOS calibrate] done.")


if __name__ == "__main__":
    main()
