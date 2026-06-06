"""
TAGMOS input-method routing audit (Rule R9)

Inspect the EC-coverage of each framework axis on YOUR cohort and recommend
input-method routing decisions (EC vs TAX fallback).

Rule R9 (from paper 1 Methods §Per-axis input-method routing):
  - If an axis has EC coverage ≥ 60 % of its dictionary, use EC input method.
  - If coverage < 60 %, route to TAX fallback (if available).
  - The BUT sub-axis is always TAX-first by design.

This script reports per-axis EC coverage so you can verify whether your
cohort meets the routing thresholds.
"""

import argparse, sys
from pathlib import Path
import json
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--ec-counts", required=True)
    ap.add_argument("--calibration", required=True,
                     help="Calibration JSON (contains ec_axis_dictionary_used)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.calibration) as f:
        cal = json.load(f)
    ec_dict = cal["ec_axis_dictionary_used"]

    ec_df = pd.read_csv(args.ec_counts, sep="\t")
    ec_col = ec_df.columns[0]
    ec_df = ec_df.set_index(ec_col)
    ec_df.index = ec_df.index.astype(str)
    present = set(ec_df.index)

    rows = []
    for sub_axis, ec_list in ec_dict.items():
        ec_list = list(ec_list)
        n_total = len(ec_list)
        n_present = sum(1 for e in ec_list if e in present)
        # also check per-sample non-zero coverage
        if n_present > 0:
            sub_matrix = ec_df.loc[[e for e in ec_list if e in present]]
            n_samples_above_zero = int((sub_matrix.sum(axis=0) > 0).sum())
            sample_coverage_pct = 100 * n_samples_above_zero / sub_matrix.shape[1]
        else:
            sample_coverage_pct = 0
        coverage_pct = 100 * n_present / n_total if n_total > 0 else 0
        rule_R9_input = "EC" if coverage_pct >= 60 else "TAX_fallback_needed"
        if sub_axis.startswith("BUT_"):
            rule_R9_input = "TAX_first_by_design"
        rows.append({
            "sub_axis": sub_axis,
            "n_EC_in_dict": n_total,
            "n_EC_present_in_cohort": n_present,
            "EC_coverage_pct": round(coverage_pct, 1),
            "sample_coverage_pct": round(sample_coverage_pct, 1),
            "rule_R9_recommendation": rule_R9_input,
        })

    res = pd.DataFrame(rows)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_path, sep="\t", index=False)

    print(f"[R9] per-sub-axis EC coverage and routing recommendation:")
    print()
    print(res.to_string(index=False))
    print()
    print(f"[R9] saved to: {out_path}")
    weak = res[res["EC_coverage_pct"] < 60]
    if not weak.empty:
        print(f"[R9] WARNING: {len(weak)} sub-axes have EC coverage < 60 % "
              f"and should use TAX fallback (or be excluded from clinical claims).")


if __name__ == "__main__":
    main()
