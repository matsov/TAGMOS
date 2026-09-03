#!/usr/bin/env python3
"""Convert curatedMetagenomicData gene-family profiles into EC coordinates.

curatedMetagenomicData distributes HUMAnN3 UniRef90 gene-family abundances in
reads per kilobase (RPK) rather than ready-made EC profiles. This script builds
the EC matrix that the scripts in ../analysis/ expect, by the procedure stated
in the paper's Methods:

  * HUMAnN's official level-4 EC-to-UniRef90 table (`map_level4ec_uniref90`) is
    inverted into a UniRef90-to-EC dictionary;
  * each UniRef90 family contributes its abundance to every EC it represents,
    multi-functional families being kept as such rather than forced to a single
    assignment;
  * abundances are summed over HUMAnN's taxonomic stratification, so that only
    the community-level total of each family is used;
  * per-study matrices are merged on the union of ECs, absent enzymes being
    zero.

Inputs
------
--map      HUMAnN's `map_level4ec_uniref90.txt.gz`, distributed with HUMAnN's
           utility mapping files. One line per EC: the EC identifier, then the
           UniRef90 families that carry it, tab-separated.
--genes    one or more HUMAnN3 gene-family tables, one per study. Each is
           tab-separated with a first column of feature names and one column
           per sample. Feature names are UniRef90 identifiers, optionally
           stratified as `UNIREF|g__Genus.s__Species`; stratified rows are
           summed into their unstratified family.
--out      path of the EC matrix to write: first column `sample_id`, one
           further column per EC, one row per sample.

Example
-------
    python build_ec_matrix.py \\
        --map map_level4ec_uniref90.txt.gz \\
        --genes cmd_genefamilies/*.tsv \\
        --out ec_matrix.tsv

Note
----
This implements the conversion as the Methods describe it. It is not the
in-house script that produced the matrix analysed in the paper, which is part
of the Wellmicro pipeline and is not distributed; it is provided so that the
same matrix can be rebuilt from public inputs. Run it once against your own
data before relying on it.
"""
import argparse
import gzip
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd


def load_uniref_to_ec(path):
    """Invert HUMAnN's EC-to-UniRef90 table into UniRef90 -> [EC, ...]."""
    opener = gzip.open if path.endswith(".gz") else open
    u2e = defaultdict(list)
    n_ec = 0
    with opener(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            ec, families = parts[0], parts[1:]
            n_ec += 1
            for fam in families:
                if fam:
                    u2e[fam].append(ec)
    if not u2e:
        sys.exit(f"{path}: no mappings read; is this HUMAnN's map_level4ec_uniref90 file?")
    multi = sum(1 for v in u2e.values() if len(v) > 1)
    print(f"[map] {n_ec} ECs, {len(u2e)} UniRef90 families "
          f"({multi} of them multi-functional, kept as such)", flush=True)
    return u2e


def genes_to_ec(path, u2e):
    """One study's gene-family table -> DataFrame of samples x EC."""
    df = pd.read_csv(path, sep="\t", low_memory=False)
    feat = df.columns[0]
    df[feat] = df[feat].astype(str)
    # sum HUMAnN's taxonomic stratification into the community-level family
    df[feat] = df[feat].str.split("|").str[0]
    df = df.groupby(feat, sort=False).sum(numeric_only=True)

    samples = df.columns.tolist()
    acc = defaultdict(lambda: np.zeros(len(samples)))
    hit = 0
    for fam, row in zip(df.index, df.to_numpy()):
        ecs = u2e.get(fam)
        if not ecs:
            continue
        hit += 1
        for ec in ecs:
            acc[ec] += row
    if not acc:
        print(f"[{os.path.basename(path)}] no family mapped to an EC; skipped",
              file=sys.stderr, flush=True)
        return None
    out = pd.DataFrame(acc, index=samples)
    print(f"[{os.path.basename(path)}] {len(samples)} samples, "
          f"{hit}/{len(df)} families mapped, {out.shape[1]} ECs", flush=True)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--map", required=True, help="HUMAnN map_level4ec_uniref90 file")
    p.add_argument("--genes", required=True, nargs="+", help="HUMAnN3 gene-family tables")
    p.add_argument("--out", required=True, help="EC matrix to write")
    a = p.parse_args()

    u2e = load_uniref_to_ec(a.map)
    parts = [x for x in (genes_to_ec(g, u2e) for g in a.genes) if x is not None]
    if not parts:
        sys.exit("no study produced an EC matrix")

    # merge on the union of ECs; an enzyme absent from a study is zero there
    ec = pd.concat(parts, axis=0, sort=True).fillna(0.0)
    ec.index.name = "sample_id"
    ec = ec.reset_index()
    ec.to_csv(a.out, sep="\t", index=False)
    print(f"[out] {a.out}: {ec.shape[0]} samples x {ec.shape[1] - 1} ECs", flush=True)

    nec = pd.DataFrame({"sample_id": ec["sample_id"],
                        "n_ec": (ec.drop(columns=["sample_id"]) > 0).sum(axis=1)})
    rich = os.path.join(os.path.dirname(a.out) or ".", "richness.tsv")
    nec.to_csv(rich, sep="\t", index=False)
    print(f"[out] {rich}: per-sample enzyme counts, pass to the analysis scripts "
          f"with --richness", flush=True)


if __name__ == "__main__":
    main()
