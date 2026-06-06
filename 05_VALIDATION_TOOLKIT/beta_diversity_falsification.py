"""
TAGMOS beta-diversity falsification toolkit

Reproduce on your own cohort the falsification finding of the primary paper:
the 36-cell TAGMOS partition explains substantially more compositional
variance than any individual binary outcome on the same Bray-Curtis matrix.

For each binary outcome and for the 36-cell TAGMOS partition, this script
computes the PERMANOVA pseudo-F and the R² (variance explained) on the
Bray-Curtis distance matrix built from the EC count table. The 36-cell
partition's R² is then compared against each outcome's R² and the
ratio is reported.

Published Italian-RWE benchmark: 38× to 94× variance advantage of the
36-cell partition over each individual clinical outcome.
"""

import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd


def bray_curtis_matrix(X: pd.DataFrame) -> np.ndarray:
    """Bray-Curtis dissimilarity. X = sample × feature, row-sum-normalised."""
    X = X.div(X.sum(axis=1).replace(0, np.nan), axis=0).fillna(0).values
    n = X.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            num = np.abs(X[i] - X[j]).sum()
            den = (X[i] + X[j]).sum()
            d = num / den if den > 0 else 0.0
            D[i, j] = D[j, i] = d
    return D


def permanova_r2(D: np.ndarray, groups: np.ndarray, n_perm: int = 99) -> dict:
    """Simplified PERMANOVA: F + R² on a 1-factor grouping. n_perm permutation
    permutation null. groups must be 1-D categorical array, same length as D."""
    n = D.shape[0]
    groups = np.asarray(groups)
    valid = ~pd.isna(groups)
    D = D[np.ix_(valid, valid)]
    groups = groups[valid]
    n = len(groups)
    if n < 4:
        return {"F": np.nan, "R2": np.nan, "p": np.nan}

    def _ss(D, groups):
        SST = (D ** 2).sum() / (2 * n)
        SSW = 0.0
        for g in np.unique(groups):
            idx = np.where(groups == g)[0]
            n_g = len(idx)
            if n_g < 2: continue
            sub = D[np.ix_(idx, idx)]
            SSW += (sub ** 2).sum() / (2 * n_g)
        SSA = SST - SSW
        return SST, SSA, SSW

    SST, SSA, SSW = _ss(D, groups)
    a = len(np.unique(groups))
    if SSW <= 0 or a < 2:
        return {"F": np.nan, "R2": np.nan, "p": np.nan}
    F = (SSA / (a - 1)) / (SSW / (n - a))
    R2 = SSA / SST if SST > 0 else np.nan

    # permutation null
    n_ge = 1
    rng = np.random.default_rng(seed=42)
    for _ in range(n_perm):
        perm = rng.permutation(groups)
        _, SSA_p, SSW_p = _ss(D, perm)
        if SSW_p <= 0: continue
        F_p = (SSA_p / (a - 1)) / (SSW_p / (n - a))
        if F_p >= F:
            n_ge += 1
    p = n_ge / (n_perm + 1)
    return {"F": float(F), "R2": float(R2), "p": float(p)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--ec-counts", required=True)
    ap.add_argument("--classification", required=True)
    ap.add_argument("--outcomes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-perm", type=int, default=99)
    args = ap.parse_args()

    print(f"[β-falsif] loading EC counts...")
    ec = pd.read_csv(args.ec_counts, sep="\t").set_index(
        pd.read_csv(args.ec_counts, sep="\t").columns[0]
    )
    ec = ec.T  # sample × EC

    clf = pd.read_csv(args.classification, sep="\t")
    out_df = pd.read_csv(args.outcomes, sep="\t")
    df = clf.merge(out_df, on="sample_id", how="inner")

    common = sorted(set(ec.index) & set(df["sample_id"]))
    print(f"[β-falsif] aligning to {len(common)} samples")
    ec = ec.loc[common]
    df = df.set_index("sample_id").loc[common].reset_index()

    print(f"[β-falsif] computing Bray-Curtis matrix...")
    D = bray_curtis_matrix(ec)

    print(f"[β-falsif] PERMANOVA on 36-cell partition...")
    cell3d = df["cell_3D"].values
    cell_res = permanova_r2(D, cell3d, n_perm=args.n_perm)
    print(f"[β-falsif] 36-cell partition: F={cell_res['F']:.3f}, "
          f"R²={cell_res['R2']:.4f}, p={cell_res['p']:.3f}")

    tagmos_cols = set(clf.columns)
    outcome_cols = [c for c in out_df.columns
                    if c != "sample_id" and c not in tagmos_cols
                    and set(out_df[c].dropna().unique()).issubset({0, 1, True, False})]

    rows = [{"grouping": "TAGMOS_36cell", "F": cell_res["F"],
             "R2": cell_res["R2"], "p": cell_res["p"],
             "advantage_vs_36cell": 1.0}]
    for oc in outcome_cols:
        print(f"[β-falsif] PERMANOVA on outcome {oc}...")
        groups = df[oc].astype(float).values
        out_res = permanova_r2(D, groups, n_perm=args.n_perm)
        adv = (cell_res["R2"] / out_res["R2"]) if out_res["R2"] and out_res["R2"] > 0 else np.nan
        rows.append({"grouping": oc, "F": out_res["F"], "R2": out_res["R2"],
                     "p": out_res["p"], "advantage_vs_36cell": adv})

    res = pd.DataFrame(rows)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_path, sep="\t", index=False)
    print()
    print(res.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print()
    print(f"[β-falsif] saved to: {out_path}")
    print(f"[β-falsif] Published Italian-RWE benchmark: 38× to 94× advantage of 36-cell partition.")


if __name__ == "__main__":
    main()
