"""
TAGMOS Cliff's δ on continuous Engine z (case vs control)

For each binary outcome, compute Cliff's δ between cases and controls
on the continuous Engine z axis. Negative δ means cases are shifted toward
dysbiotic (lower Engine z) relative to controls.

|δ| interpretation: < 0.15 trivial · 0.15-0.33 small · 0.33-0.47 medium · > 0.47 large.
"""

import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats


def cliff_delta(a, b):
    """Vectorised Cliff δ between two 1-D arrays."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    m, n = len(a), len(b)
    if m == 0 or n == 0:
        return np.nan
    diff = a[:, None] - b[None, :]
    return float(np.sign(diff).sum() / (m * n))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--classification", required=True)
    ap.add_argument("--outcomes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--axis-col", default="ENG_z",
                     help="Continuous axis column to use (default: ENG_z)")
    ap.add_argument("--outcome-cols", nargs="+", default=None)
    args = ap.parse_args()

    clf = pd.read_csv(args.classification, sep="\t")
    if args.axis_col not in clf.columns:
        sys.exit(f"ERROR: column '{args.axis_col}' not found in classification")
    out_df = pd.read_csv(args.outcomes, sep="\t")
    df = clf.merge(out_df, on="sample_id", how="inner")

    if args.outcome_cols is None:
        tagmos_cols = set(clf.columns)
        outcome_cols = [c for c in out_df.columns
                        if c != "sample_id" and c not in tagmos_cols
                        and set(out_df[c].dropna().unique()).issubset({0, 1, True, False})]
    else:
        outcome_cols = args.outcome_cols

    rows = []
    for oc in outcome_cols:
        sub = df[[args.axis_col, oc]].dropna()
        case_vals = sub.loc[sub[oc] == 1, args.axis_col].values
        ctrl_vals = sub.loc[sub[oc] == 0, args.axis_col].values
        delta = cliff_delta(case_vals, ctrl_vals)
        if len(case_vals) > 0 and len(ctrl_vals) > 0:
            try:
                u_stat, mw_p = stats.mannwhitneyu(case_vals, ctrl_vals,
                                                    alternative="two-sided")
            except Exception:
                u_stat, mw_p = np.nan, np.nan
        else:
            u_stat, mw_p = np.nan, np.nan
        rows.append({
            "outcome": oc,
            "n_case": len(case_vals), "n_ctrl": len(ctrl_vals),
            "median_axis_case": float(np.median(case_vals)) if len(case_vals) > 0 else np.nan,
            "median_axis_ctrl": float(np.median(ctrl_vals)) if len(ctrl_vals) > 0 else np.nan,
            "cliff_delta": delta,
            "mannwhitney_p": mw_p,
        })

    res = pd.DataFrame(rows)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_path, sep="\t", index=False)
    print(res.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print(f"[TAGMOS Cliff-delta] saved to: {out_path}")


if __name__ == "__main__":
    main()
