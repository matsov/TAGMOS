"""
TAGMOS disease association · T4-vs-T1 odds ratio with 95% CI + BH q-value

Given a TAGMOS classification TSV (output of classify_local_*.py) and an
outcomes TSV (sample_id + binary outcome column(s)), compute for each outcome:
  - per-tier prevalence (T1 / T2 / T3 / T4)
  - T4-vs-T1 odds ratio + 95% Wald CI (Haldane-Anscombe correction if any cell < 5)
  - Cochran-Armitage T1->T4 trend Z statistic + p-value
  - Benjamini-Hochberg q-value across all outcomes tested jointly

Usage:
    python run_disease_OR.py \\
        --classification my_classification.tsv \\
        --outcomes your_outcomes.tsv \\
        --out disease_OR_results.tsv \\
        --outcome-cols T2D Celiac IBD     (optional · default: all non-id columns)

Outcomes TSV format:
    sample_id   T2D   Celiac   IBD   age   sex   BMI
    sample_001  0     0        0     45    F     24
    sample_002  1     0        0     52    M     27
    ...
"""

import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


TIERS = ["T1_EUBIOTIC", "T2_PRESERVED", "T3_ALTERED", "T4_DYSBIOTIC"]


def or_with_95ci(case_T4, ctrl_T4, case_T1, ctrl_T1):
    """Standard 2x2 odds ratio + Wald 95% CI. Haldane-Anscombe correction
    (+ 0.5 to all cells) applied if any cell < 5."""
    cells = [case_T4, ctrl_T4, case_T1, ctrl_T1]
    if min(cells) < 5:
        case_T4 += 0.5; ctrl_T4 += 0.5; case_T1 += 0.5; ctrl_T1 += 0.5
    odds_case = case_T4 / case_T1 if case_T1 > 0 else np.nan
    odds_ctrl = ctrl_T4 / ctrl_T1 if ctrl_T1 > 0 else np.nan
    if np.isnan(odds_case) or np.isnan(odds_ctrl) or odds_ctrl == 0:
        return np.nan, np.nan, np.nan
    OR = odds_case / odds_ctrl
    if min(case_T4, ctrl_T4, case_T1, ctrl_T1) <= 0:
        return OR, np.nan, np.nan
    se_log_or = np.sqrt(1/case_T4 + 1/ctrl_T4 + 1/case_T1 + 1/ctrl_T1)
    ci_lo = OR * np.exp(-1.96 * se_log_or)
    ci_hi = OR * np.exp(+1.96 * se_log_or)
    return OR, ci_lo, ci_hi


def cochran_armitage_trend(cases_per_tier, total_per_tier):
    """Cochran-Armitage trend test Z statistic across an ordered tier sequence.
    Scores (1, 2, 3, 4) assigned to T1, T2, T3, T4."""
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
    ap.add_argument("--outcome-cols", nargs="+", default=None,
                     help="Subset of outcome column names. "
                          "Default: all numeric (0/1) columns of the outcomes TSV.")
    args = ap.parse_args()

    print(f"[TAGMOS disease-OR] reading classification: {args.classification}")
    clf = pd.read_csv(args.classification, sep="\t")
    print(f"[TAGMOS disease-OR] reading outcomes:       {args.outcomes}")
    out_df = pd.read_csv(args.outcomes, sep="\t")
    df = clf.merge(out_df, on="sample_id", how="inner")
    print(f"[TAGMOS disease-OR] joined samples:         {len(df)}")
    if len(df) == 0:
        sys.exit("ERROR: no overlap between classification and outcomes on sample_id")

    # determine outcome columns
    if args.outcome_cols is None:
        # any 0/1 binary column that is NOT a TAGMOS column and not metadata
        tagmos_cols = {c for c in clf.columns}
        candidates = []
        for c in out_df.columns:
            if c == "sample_id": continue
            if c in tagmos_cols: continue
            vals = set(out_df[c].dropna().unique())
            if vals.issubset({0, 1}) or vals.issubset({0, 1, True, False}):
                candidates.append(c)
        outcome_cols = candidates
    else:
        outcome_cols = args.outcome_cols
    print(f"[TAGMOS disease-OR] testing {len(outcome_cols)} outcomes: {outcome_cols}")

    rows = []
    raw_p = []
    for oc in outcome_cols:
        df_oc = df[["ENG_tier", oc]].dropna()
        n_total = len(df_oc)
        n_case = int(df_oc[oc].sum())

        cases_per_tier = []
        totals_per_tier = []
        prev_per_tier = []
        for t in TIERS:
            mask = (df_oc["ENG_tier"] == t)
            tot = int(mask.sum())
            cas = int(df_oc.loc[mask, oc].sum())
            totals_per_tier.append(tot)
            cases_per_tier.append(cas)
            prev_per_tier.append(100 * cas / tot if tot > 0 else np.nan)

        case_T1, case_T4 = cases_per_tier[0], cases_per_tier[3]
        tot_T1, tot_T4 = totals_per_tier[0], totals_per_tier[3]
        ctrl_T1 = tot_T1 - case_T1
        ctrl_T4 = tot_T4 - case_T4
        OR, ci_lo, ci_hi = or_with_95ci(case_T4, ctrl_T4, case_T1, ctrl_T1)
        ca_Z, ca_p = cochran_armitage_trend(cases_per_tier, totals_per_tier)
        prev_factor = (prev_per_tier[3] / prev_per_tier[0]) if prev_per_tier[0] > 0 else np.nan

        rows.append({
            "outcome": oc,
            "n_total": n_total, "n_case": n_case,
            "prev_T1_pct": prev_per_tier[0], "prev_T2_pct": prev_per_tier[1],
            "prev_T3_pct": prev_per_tier[2], "prev_T4_pct": prev_per_tier[3],
            "prevalence_factor_T4_vs_T1": prev_factor,
            "OR_T4_vs_T1": OR, "OR_95CI_lo": ci_lo, "OR_95CI_hi": ci_hi,
            "CA_trend_Z": ca_Z, "CA_trend_p": ca_p,
        })
        raw_p.append(ca_p if not np.isnan(ca_p) else 1.0)

    res = pd.DataFrame(rows)
    if len(raw_p) > 1:
        _, qvals, _, _ = multipletests(raw_p, method="fdr_bh")
        res["BH_q"] = qvals
    else:
        res["BH_q"] = raw_p

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_path, sep="\t", index=False)

    print(f"[TAGMOS disease-OR] writing results to:    {out_path}")
    print()
    print(res.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print()
    print(f"[TAGMOS disease-OR] done. Cite TAGMOS — see ../HOW_TO_CITE.md")


if __name__ == "__main__":
    main()
