"""
TAGMOS classify · Version A · single-output 14 axes only
========================================================

Apply a previously-computed calibration JSON to a new EC-count matrix.
Output: per-sample z-score for each of the 14 single-output framework axes,
plus tier label for Engine z and tertile class label for SYN_cof / SYN_carb /
channels.

Usage
-----
    python classify_local_single.py \\
        --ec-counts your_samples_EC_counts.tsv \\
        --calibration my_calibration.json \\
        --out my_classification.tsv

License: TAGMOS Public Research License v1.0 — see ../LICENSE.md
Mandatory citation: see ../HOW_TO_CITE.md
"""

import argparse, json, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from classifier_single import (
    AXIS_EC_DICT, AXIS_FORMULA, FRAMEWORK_AXES_14,
    compute_sub_axis_log1p, compute_axis_raw,
    assign_quartile_tier, assign_tertile_class,
)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--ec-counts", required=True)
    ap.add_argument("--calibration", required=True,
                    help="Calibration JSON from calibrate_local_single.py")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    print(f"[TAGMOS classify] reading calibration:         {args.calibration}")
    with open(args.calibration) as f:
        cal = json.load(f)
    print(f"[TAGMOS classify] calibration_id:              {cal['calibration_id']}")
    print(f"[TAGMOS classify] calibration sha256:          {cal.get('calibration_sha256','?')[:16]}...")

    print(f"[TAGMOS classify] reading EC counts:           {args.ec_counts}")
    df = pd.read_csv(args.ec_counts, sep="\t")
    ec_col = df.columns[0]
    df = df.set_index(ec_col)
    df.index = df.index.astype(str)
    print(f"[TAGMOS classify] n samples:                   {df.shape[1]}")

    # --- recompute sub-axis log1p sums ---
    sub_axis_matrix = compute_sub_axis_log1p(df)

    # --- apply calibration mu/sigma ---
    sub_axis_mu_sigma = {k: tuple(v) for k, v in cal["per_sub_axis_mu_sigma"].items()}
    axis_raw = compute_axis_raw(sub_axis_matrix, sub_axis_mu_sigma)

    axis_mu_sigma = {k: tuple(v) for k, v in cal["per_axis_mu_sigma"].items()}
    axis_z = {}
    for axis in axis_raw.index:
        mu, sigma = axis_mu_sigma.get(axis, (0.0, 1.0))
        if sigma <= 0:
            sigma = 1.0
        axis_z[axis] = (axis_raw.loc[axis] - mu) / sigma
    axis_z_df = pd.DataFrame(axis_z).T  # axis × sample

    # --- assign tiers / classes ---
    out = pd.DataFrame(index=df.columns)
    out.index.name = "sample_id"
    for axis in FRAMEWORK_AXES_14:
        out[f"{axis}_z"] = axis_z_df.loc[axis].values

    # Engine z tier (quartile)
    eng_cuts = [cal["engine_quartile_cuts"]["Q25"],
                cal["engine_quartile_cuts"]["Q50"],
                cal["engine_quartile_cuts"]["Q75"]]
    out["ENG_tier"] = assign_quartile_tier(
        axis_z_df.loc["ENG"],
        cuts=eng_cuts,
        labels=cal["tier_labels"]["ENG"],
    ).values

    # SYN_cof / SYN_carb / channels tertile
    tertile_labels = {
        "SYN_cof":  cal["tier_labels"]["SYN_cof"],
        "SYN_carb": cal["tier_labels"]["SYN_carb"],
    }
    for axis in FRAMEWORK_AXES_14:
        if axis == "ENG":
            continue
        if axis not in cal["tertile_cuts"]:
            continue
        cuts = [cal["tertile_cuts"][axis]["Q33"],
                cal["tertile_cuts"][axis]["Q66"]]
        if axis in tertile_labels:
            labels = tertile_labels[axis]
        else:
            labels = [f"{axis}_LOW", f"{axis}_MID", f"{axis}_HIGH"]
        out[f"{axis}_class"] = assign_tertile_class(
            axis_z_df.loc[axis], cuts=cuts, labels=labels,
        ).values

    # 36-cell label (3D state-space) — derived
    out["cell_3D"] = (out["ENG_tier"] + " × "
                     + out["SYN_cof_class"] + " × "
                     + out["SYN_carb_class"])

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.reset_index().to_csv(out_path, sep="\t", index=False)

    # --- summary
    print(f"[TAGMOS classify] writing classification to:   {out_path}")
    print()
    print(f"[TAGMOS classify] Engine tier distribution (% of samples):")
    eng_dist = out["ENG_tier"].value_counts(normalize=True).sort_index()
    for tier, pct in eng_dist.items():
        print(f"    {tier:<14} {100*pct:5.1f} %")
    print()
    print(f"[TAGMOS classify] 36-cell occupancy: "
          f"{out['cell_3D'].nunique()} populated cells of 36")
    print(f"[TAGMOS classify] done. Cite TAGMOS — see ../HOW_TO_CITE.md")


if __name__ == "__main__":
    main()
