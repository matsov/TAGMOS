"""Shared input handling for the TAGMOS analysis scripts.

Every script in this directory takes its inputs from the command line rather
than from hard-coded paths, and takes the definition of the metabolic axes from
a panel file that the user supplies (see ../panel/README.md). None of the
scripts contains, or needs, the enzyme lists used in the paper: they operate on
whatever panel they are given.

Expected inputs
---------------
--ec        EC abundance matrix, tab-separated. First column `sample_id`, one
            further column per EC number, one row per sample. This is the
            matrix produced by ../ec_matrix/build_ec_matrix.py, or any matrix
            in the same shape.
--species   species abundance matrix, tab-separated. First column `feature`
            (the species name), one column per sample. Column names may be
            either `study|sample_id` or a bare `sample_id`.
--meta      sample metadata, tab-separated, one row per sample, with at least
            the columns `sample_id`, `study_name`, `study_condition`,
            `disease_subtype`, `body_site`, `age`, `age_category`, `gender`,
            `BMI`, `country`, `number_reads`, `non_westernized`,
            `antibiotics_current_use`, `subject_id` and `disease`. The `sampleMetadata` table of
            curatedMetagenomicData has this shape; `disease_subtype` is what
            separates Crohn's disease from ulcerative colitis within IBD.
--panel     axis panel, JSON (see ../panel/README.md).
--richness  optional, tab-separated with columns `sample_id` and `n_ec`: the
            number of enzymes detected per sample, used as a technical
            covariate throughout. If omitted it is computed from --ec, which is
            correct only when --ec holds the full enzyme matrix rather than a
            panel-restricted subset.
--out       directory to write results into.
"""
import argparse
import json
import os
import sys

import pandas as pd


def cli(description, ec=False, species=False, panel=False, richness=False, guild=False):
    p = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    if ec:
        p.add_argument("--ec", required=True, help="EC abundance matrix (sample_id x EC)")
    if species:
        p.add_argument("--species", required=True, help="species abundance matrix (feature x sample)")
    p.add_argument("--meta", required=True, help="sample metadata")
    if panel:
        p.add_argument("--panel", required=True, help="axis panel JSON")
    if guild:
        p.add_argument("--condition", default="CRC",
                       help="study_condition to contrast against controls (default CRC)")
        p.add_argument("--guild", required=True,
                       help="species guild panel: one species name per line, first column")
    if richness:
        p.add_argument("--richness", default=None, help="per-sample enzyme count (sample_id, n_ec)")
    p.add_argument("--out", required=True, help="output directory")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    return a


def load_panel(path):
    """Read an axis panel and check it against the schema.

    Returns the panel as a dict. Keys are axis names; each value carries at
    least a valence, and either an `ecs` list (an axis scored as the mean of
    its z-scored enzymes) or a `num`/`den` pair (an axis scored as a
    log-ratio). Axes with neither are ignored by the enzyme-level scripts.
    """
    with open(path) as fh:
        ax = json.load(fh)
    if not isinstance(ax, dict) or not ax:
        sys.exit(f"{path}: expected a non-empty JSON object of axis definitions")
    ax = {n: v for n, v in ax.items() if not n.startswith("_")}   # keys starting with _ are comments
    for name, v in ax.items():
        if not isinstance(v, dict):
            sys.exit(f"{path}: axis {name!r} must be an object")
        if not (v.get("ecs") or (v.get("num") and v.get("den"))):
            continue
        if v.get("valence") not in ("protective", "danger", "context", None):
            sys.exit(f"{path}: axis {name!r} has an unrecognised valence {v.get('valence')!r}")
    usable = [n for n, v in ax.items() if v.get("ecs")]
    if not usable:
        sys.exit(f"{path}: no axis carries an `ecs` list, nothing to score")
    print(f"[panel] {len(ax)} axes, {len(usable)} with an enzyme list", flush=True)
    return ax


def load_ec(path, keep=None):
    """Read the EC matrix, optionally only the columns in `keep`."""
    if keep is not None:
        hdr = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
        use = [e for e in keep if e in hdr]
        missing = [e for e in keep if e not in hdr]
        if missing:
            print(f"[ec] {len(missing)} of {len(keep)} panel enzymes are absent from the matrix",
                  file=sys.stderr, flush=True)
        df = pd.read_csv(path, sep="\t", usecols=["sample_id"] + use, low_memory=False)
    else:
        df = pd.read_csv(path, sep="\t", low_memory=False)
    df = df.drop_duplicates("sample_id").set_index("sample_id")
    print(f"[ec] {df.shape[0]} samples x {df.shape[1]} enzymes", flush=True)
    return df


def load_meta(path):
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    need = {"sample_id", "study_name", "study_condition"}
    missing = need - set(df.columns)
    if missing:
        sys.exit(f"{path}: metadata is missing required column(s): {sorted(missing)}")
    df = df.drop_duplicates("sample_id").set_index("sample_id")
    print(f"[meta] {df.shape[0]} samples, {df['study_name'].nunique()} studies", flush=True)
    return df


def load_richness(path, ec_df=None):
    """Per-sample enzyme count, read from file or derived from the EC matrix."""
    if path:
        s = pd.read_csv(path, sep="\t").drop_duplicates("sample_id").set_index("sample_id")["n_ec"]
        return s
    if ec_df is None:
        sys.exit("--richness is required when no EC matrix is given")
    print("[richness] --richness not given; counting detected enzymes in --ec. "
          "This is the intended covariate only if --ec is the full enzyme matrix.",
          file=sys.stderr, flush=True)
    return (ec_df > 0).sum(axis=1).rename("n_ec")
