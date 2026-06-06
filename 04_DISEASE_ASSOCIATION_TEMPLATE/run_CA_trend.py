"""
TAGMOS Cochran-Armitage trend test · standalone version
"""

import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


TIERS = ["T1_EUBIOTIC", "T2_PRESERVED", "T3_ALTERED", "T4_DYSBIOTIC"]


def cochran_armitage(cases_per_tier, total_per_tier):
    scores = np.array([1, 2, 3, 4], dtype=float)
    cases = np.array(cases_per_tier, dtype=float)
    totals = np.array(total_per_tier, dtype=float)
    N = totals.sum()
    if N <= 0:
        return np.nan, np.nan
    p = cases.sum() / N
    T = np.sum(scores * (cases - totals * p))
    var = p * (1 - p) * (np.sum(totals * scores ** 2)
                         - (np.sum(totals * scores)) ** 2 / N)
    if var <= 0:
        return np.nan, np.nan
    Z = T / np.sqrt(var)
    p_val = 2 * (1 - stats.norm.cdf(abs(Z)))
    return float(Z), float(p_val)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--classification", required=True)
    ap.add_argument("--outcomes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--outcome-cols", nargs="+", default=None)
    args = ap.parse_args()

    clf = pd.read_csv(args.classification, sep="\t")
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
    raw_p = []
    for oc in outcome_cols:
        df_oc = df[["ENG_tier", oc]].dropna()
        cases = [int(df_oc[df_oc["ENG_tier"] == t][oc].sum()) for t in TIERS]
        totals = [int((df_oc["ENG_tier"] == t).sum()) for t in TIERS]
        Z, p = cochran_armitage(cases, totals)
        rows.append({
            "outcome": oc,
            "n_T1": totals[0], "n_T2": totals[1], "n_T3": totals[2], "n_T4": totals[3],
            "cases_T1": cases[0], "cases_T2": cases[1], "cases_T3": cases[2], "cases_T4": cases[3],
            "CA_Z": Z, "CA_p": p,
        })
        raw_p.append(p if not np.isnan(p) else 1.0)

    res = pd.DataFrame(rows)
    if len(raw_p) > 1:
        _, qvals, _, _ = multipletests(raw_p, method="fdr_bh")
        res["BH_q"] = qvals
    else:
        res["BH_q"] = raw_p

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_path, sep="\t", index=False)
    print(res.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print(f"[TAGMOS CA-trend] saved to: {out_path}")


if __name__ == "__main__":
    main()
