"""
TAGMOS dimensional saturation test

Test on YOUR cohort whether 3 axes (Engine × SYN_cof × SYN_carb) is the
dimensional sweet spot, as reported in paper 1 Methods §Architectural
saturation testing.

Reports the bias-corrected Cramér's V* (Bergsma 2013) for each dimensionality
k = 1, 2, 3 on the binary outcomes provided. Per the published finding,
k = 3 should be the saturation point (gain at k = 4 minimal; per-cell
sample density drops below the n ≥ 30 robust floor).
"""

import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd


def cramers_v_bias_corrected(contingency: np.ndarray) -> float:
    """Bergsma 2013 bias-corrected Cramér's V*."""
    if contingency.sum() == 0:
        return np.nan
    chi2 = 0
    n = contingency.sum()
    row_sums = contingency.sum(axis=1)
    col_sums = contingency.sum(axis=0)
    for i in range(contingency.shape[0]):
        for j in range(contingency.shape[1]):
            exp = row_sums[i] * col_sums[j] / n if n > 0 else 0
            if exp > 0:
                chi2 += (contingency[i, j] - exp) ** 2 / exp
    phi2 = chi2 / n if n > 0 else 0
    r, k = contingency.shape
    phi2_corr = max(0, phi2 - (r - 1) * (k - 1) / (n - 1))
    r_corr = r - (r - 1) ** 2 / (n - 1)
    k_corr = k - (k - 1) ** 2 / (n - 1)
    denom = min(r_corr - 1, k_corr - 1)
    if denom <= 0:
        return np.nan
    return float(np.sqrt(phi2_corr / denom))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--classification", required=True)
    ap.add_argument("--outcomes", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    clf = pd.read_csv(args.classification, sep="\t")
    out_df = pd.read_csv(args.outcomes, sep="\t")
    df = clf.merge(out_df, on="sample_id", how="inner")

    tagmos_cols = set(clf.columns)
    outcome_cols = [c for c in out_df.columns
                    if c != "sample_id" and c not in tagmos_cols
                    and set(out_df[c].dropna().unique()).issubset({0, 1, True, False})]

    # build groupings at increasing dimensionality
    d1 = df["ENG_tier"]
    d2 = df["ENG_tier"] + "|" + df["SYN_cof_class"]
    d3 = df["cell_3D"]

    rows = []
    for k_label, grouping in [("1D (Engine)", d1),
                                ("2D (Engine × SYN_cof)", d2),
                                ("3D (full 36-cell)", d3)]:
        for oc in outcome_cols:
            cont = pd.crosstab(grouping, df[oc]).values
            v_star = cramers_v_bias_corrected(cont)
            rows.append({"dimensionality": k_label, "outcome": oc,
                         "V_star": v_star,
                         "n_groups_populated": pd.Series(grouping).nunique()})

    res = pd.DataFrame(rows)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_path, sep="\t", index=False)

    print(f"[dim-sat] dimensionality vs Cramér V*")
    print()
    pivoted = res.pivot(index="outcome", columns="dimensionality", values="V_star")
    print(pivoted.to_string(float_format=lambda x: f"{x:.3g}"))
    print()
    print(f"[dim-sat] saved to: {out_path}")
    print(f"[dim-sat] Published Italian-RWE benchmark: 3D Cramér V* peak at saturation.")


if __name__ == "__main__":
    main()
