"""
TAGMOS classifier · Version A · single-output 14 axes only
==========================================================

This module defines:
  - AXIS_EC_DICT       : EC dictionaries per sub-axis (the framework's single
                         source of truth at runtime)
  - AXIS_FORMULA       : per-axis aggregation formula
  - TIER_PARTITION     : per-axis partition method (quartile / tertile)
  - compute_axis_raw   : reduce an EC-count matrix to per-sample axis raw scores
  - z_score_apply      : apply pre-calibrated mu/sigma to z-score each axis
  - tier_assign        : apply pre-calibrated quartile / tertile cuts to assign
                         tier or class labels

The module is imported by:
  - calibrate_local_single.py   (compute calibration JSON on your cohort)
  - classify_local_single.py    (apply a calibration JSON to a new EC matrix)

No L0 / Wellmicro proprietary numerics are present anywhere in this file or
elsewhere in this bundle. All calibration parameters are user-side.

License: TAGMOS Public Research License v1.0 — see LICENSE.md
Mandatory citation: see HOW_TO_CITE.md
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple


# =============================================================================
# AXIS_EC_DICT · per-sub-axis EC list
# -----------------------------------------------------------------------------
# Format:  sub_axis_name -> list of Level-4 EC numbers (strings).
# Each sub-axis is computed as the log1p-mean across its EC list.
# An aggregator axis (SYN_cof, SYN_carb) is the z-mean of its sub-axes.
# A directional channel (HIS, TMA) is the difference of its production and
# degradation sub-axes (see AXIS_FORMULA below).
# =============================================================================

AXIS_EC_DICT: Dict[str, List[str]] = {

    # ----- ENGINE (PRIMARY) -----
    "H_sink_acetogenesis_WL":   ["1.2.7.4", "2.3.1.169", "6.3.4.3",
                                  "1.5.1.5", "3.5.4.9", "1.5.1.20"],
    "H_sink_methanogenesis":    ["2.8.4.1", "1.5.98.1", "1.5.98.2",
                                  "1.5.98.3", "1.8.7.3"],
    "H_sink_sulfate_reduction": ["1.8.99.5", "1.8.4.8"],
    "H_sink_dysbiotic":         ["7.1.1.9", "7.1.1.7",
                                  "1.7.1.15", "1.7.2.2", "1.7.5.1", "1.7.2.1",
                                  "1.7.1.4", "1.9.3.1"],

    # ----- SYN_cof CO-PRIMARY sub-axes -----
    "B12_biosynthesis_de_novo":  ["2.1.1.151", "4.99.1.3"],
    "B12_last_step_adenosyl":    ["2.5.1.17"],
    "B12_cobalt_insertion":      ["6.6.1.2"],
    "B12_DMB_ligand":            ["2.4.2.21"],
    "B12_consumer":              ["5.4.99.2", "2.1.1.13"],
    "K2_pathway_early":          ["4.1.1.71", "6.2.1.26", "4.1.2.62"],
    "K2_pathway_late":           ["4.1.3.36", "3.1.2.28"],
    "K2_last_step":              ["2.1.1.163"],
    "POL_putrescine":            ["4.1.1.17", "4.1.1.19", "3.5.3.11"],
    "GABA":                      ["4.1.1.15"],

    # ----- SYN_carb CO-PRIMARY sub-axes -----
    "BUT_butyryl_CoA":           ["2.8.3.8"],   # complemented by TAX panel
    "PRP_propanediol":           ["4.2.1.28", "1.1.1.202"],
    "LAC_racemase":              ["5.1.2.1"],

    # ----- BA channel sub-axes -----
    "BA_BSH":                    ["3.5.1.24"],
    "BA_dehydroxylation":        ["1.3.1.115", "1.17.99.1", "1.1.1.391"],
    "BA_epimerization":          ["1.1.1.50", "1.1.1.52"],
    "BA_conjugation":            ["6.2.1.7", "2.3.1.65"],

    # ----- TRP channel sub-axes -----
    "Trp_indole":                ["4.1.99.1", "1.4.3.2", "4.1.1.74"],
    "Trp_kynurenine":            ["1.13.11.11", "1.13.11.52", "3.5.1.9"],

    # ----- PROT channel (aggregated) -----
    "PROT":                      ["4.1.1.25", "4.1.1.18", "4.1.1.17"],

    # ----- MUC channel sub-axes -----
    "Mucin_specialised":         ["3.2.1.49", "3.2.1.97", "3.2.1.96",
                                  "3.2.1.51", "3.2.1.18"],
    "Mucin_general":             ["3.2.1.22", "3.2.1.4"],

    # ----- HIS channel sub-axes -----
    "Histamine_HDC":             ["4.1.1.22"],
    "Histamine_HAL":             ["4.3.1.3"],
    "Histamine_NAT":             ["2.3.1.5"],

    # ----- TMA channel sub-axes -----
    "TMA_production":            ["4.3.99.4", "4.3.99.3", "1.7.99.1",
                                  "4.3.99.1"],
    "TMA_degradation":           ["1.14.13.148", "2.1.1.157"],

    # ----- UREM channel -----
    "UREM":                      ["4.1.99.5", "1.13.11.27"],

    # ----- IRON channel (siderophore biosynthesis) -----
    "IRON":                      ["6.3.2.39", "6.3.2.14"],

    # ----- LPS channel (hexa-acyl) -----
    "LPS":                       ["2.3.1.129", "2.3.1.241"],

    # ----- ETU channel (ethanolamine utilisation) -----
    "ETU":                       ["4.3.1.7"],

    # ----- HYS channel (H2S) -----
    "HYS":                       ["4.2.1.22", "2.5.1.48"],
}


# =============================================================================
# AXIS_FORMULA · how each framework axis is built from sub-axes
# =============================================================================

AXIS_FORMULA: Dict[str, Dict] = {

    "ENG": {
        "type": "engine_balance",
        "eubiotic": ["H_sink_acetogenesis_WL",
                     "H_sink_methanogenesis",
                     "H_sink_sulfate_reduction"],
        "dysbiotic": ["H_sink_dysbiotic"],
        "tier_partition": "quartile",
        "tier_labels": ["T4_DYSBIOTIC", "T3_ALTERED",
                        "T2_PRESERVED", "T1_EUBIOTIC"],  # low→high z
    },

    "SYN_cof": {
        "type": "z_mean_of_sub_axes",
        "sub_axes": ["B12_biosynthesis_de_novo", "B12_last_step_adenosyl",
                     "B12_cobalt_insertion", "B12_DMB_ligand",
                     "B12_consumer",
                     "K2_pathway_early", "K2_pathway_late", "K2_last_step",
                     "POL_putrescine",
                     "GABA"],
        "_group_mapping_to_4_aggregator_subs": {
            "B12": ["B12_biosynthesis_de_novo", "B12_last_step_adenosyl",
                    "B12_cobalt_insertion", "B12_DMB_ligand", "B12_consumer"],
            "K2":  ["K2_pathway_early", "K2_pathway_late", "K2_last_step"],
            "POL": ["POL_putrescine"],
            "GABA": ["GABA"],
        },
        "tier_partition": "tertile",
        "tier_labels": ["SYN_FRAGILE", "SYN_PARTIAL", "SYN_RESILIENT"],
    },

    "SYN_carb": {
        "type": "z_mean_of_sub_axes",
        "sub_axes": ["BUT_butyryl_CoA", "PRP_propanediol", "LAC_racemase"],
        "_note": ("In production the BUT sub-axis is TAX-first via a curated "
                   "species panel; in this public bundle we use the EC marker "
                   "EC 2.8.3.8 as a portable fallback. Users with a species "
                   "abundance table should overlay the TAX panel as per the "
                   "calibration recipe §BUT-TAX-fallback."),
        "tier_partition": "tertile",
        "tier_labels": ["SYN_carb_FRAGILE", "SYN_carb_PARTIAL",
                        "SYN_carb_RESILIENT"],
    },

    "BA": {
        "type": "mean_of_sub_axes",
        "sub_axes": ["BA_BSH", "BA_dehydroxylation",
                     "BA_epimerization", "BA_conjugation"],
        "tier_partition": "tertile",
    },
    "TRP": {
        "type": "mean_of_sub_axes",
        "sub_axes": ["Trp_indole", "Trp_kynurenine"],
        "tier_partition": "tertile",
    },
    "PROT": {
        "type": "log1p_z_of_sub_axis",
        "sub_axis": "PROT",
        "tier_partition": "tertile",
    },
    "MUC": {
        "type": "mean_of_sub_axes",
        "sub_axes": ["Mucin_specialised", "Mucin_general"],
        "tier_partition": "tertile",
    },
    "HIS": {
        "type": "directional",
        "positive": ["Histamine_HDC"],
        "negative": ["Histamine_HAL", "Histamine_NAT"],
        "tier_partition": "tertile",
    },
    "TMA": {
        "type": "directional",
        "positive": ["TMA_production"],
        "negative": ["TMA_degradation"],
        "tier_partition": "tertile",
    },
    "UREM": {
        "type": "log1p_z_of_sub_axis",
        "sub_axis": "UREM",
        "tier_partition": "tertile",
    },
    "IRON": {
        "type": "log1p_z_of_sub_axis",
        "sub_axis": "IRON",
        "tier_partition": "tertile",
    },
    "LPS": {
        "type": "log1p_z_of_sub_axis",
        "sub_axis": "LPS",
        "tier_partition": "tertile",
    },
    "ETU": {
        "type": "log1p_z_of_sub_axis",
        "sub_axis": "ETU",
        "tier_partition": "tertile",
    },
    "HYS": {
        "type": "log1p_z_of_sub_axis",
        "sub_axis": "HYS",
        "tier_partition": "tertile",
    },
}

FRAMEWORK_AXES_14 = ["ENG", "SYN_cof", "SYN_carb",
                     "BA", "TRP", "PROT", "MUC", "HIS", "TMA",
                     "UREM", "IRON", "LPS", "ETU", "HYS"]


# =============================================================================
# Core math
# =============================================================================

def _safe_log1p_sum(ec_counts: pd.DataFrame, ec_list: List[str]) -> pd.Series:
    """Sum log1p(counts) over the requested EC list. EC absent in the matrix
    contributes 0. Returns a per-sample Series."""
    present = [e for e in ec_list if e in ec_counts.index]
    if not present:
        return pd.Series(0.0, index=ec_counts.columns)
    sub = ec_counts.loc[present]
    return np.log1p(sub).sum(axis=0)


def compute_sub_axis_log1p(ec_counts: pd.DataFrame) -> pd.DataFrame:
    """Return a sub-axis × sample matrix of log1p-sum scores.
    ec_counts must have EC ids as index, samples as columns."""
    rows = {}
    for sub_axis, ec_list in AXIS_EC_DICT.items():
        rows[sub_axis] = _safe_log1p_sum(ec_counts, ec_list)
    return pd.DataFrame(rows).T  # sub_axes × samples


def z_score_against_mu_sigma(values: pd.Series,
                              mu: float,
                              sigma: float) -> pd.Series:
    """Z-score against pre-computed mu/sigma. sigma ≤ 0 falls back to 1."""
    if sigma is None or sigma <= 0 or np.isnan(sigma):
        sigma = 1.0
    return (values - mu) / sigma


def compute_axis_raw(sub_axis_matrix: pd.DataFrame,
                     sub_axis_mu_sigma: Dict[str, Tuple[float, float]]
                     ) -> pd.DataFrame:
    """Build the 14 framework axis raw scores per sample from sub-axis scores
    using AXIS_FORMULA. Returns axis × sample matrix.

    Notes
    -----
    For aggregator axes the sub-axis scores are first z-scored against the
    sub-axis-level mu/sigma (provided by the calibration JSON) before being
    averaged. For directional axes the production-side and degradation-side
    z-scored sub-axes are subtracted. For single-component channels the
    log1p-z of the relevant sub-axis is used directly.
    """
    out = {}
    samples = sub_axis_matrix.columns

    def _z(sub_axis_name: str) -> pd.Series:
        mu, sigma = sub_axis_mu_sigma.get(sub_axis_name, (0.0, 1.0))
        return z_score_against_mu_sigma(sub_axis_matrix.loc[sub_axis_name],
                                         mu, sigma)

    for axis in FRAMEWORK_AXES_14:
        rule = AXIS_FORMULA[axis]

        if rule["type"] == "engine_balance":
            eub_z = pd.concat([_z(s) for s in rule["eubiotic"]], axis=1).mean(axis=1)
            dys_z = pd.concat([_z(s) for s in rule["dysbiotic"]], axis=1).mean(axis=1)
            out[axis] = eub_z - dys_z

        elif rule["type"] == "z_mean_of_sub_axes":
            zs = pd.concat([_z(s) for s in rule["sub_axes"]], axis=1)
            out[axis] = zs.mean(axis=1)

        elif rule["type"] == "mean_of_sub_axes":
            zs = pd.concat([_z(s) for s in rule["sub_axes"]], axis=1)
            out[axis] = zs.mean(axis=1)

        elif rule["type"] == "log1p_z_of_sub_axis":
            out[axis] = _z(rule["sub_axis"])

        elif rule["type"] == "directional":
            pos = pd.concat([_z(s) for s in rule["positive"]], axis=1).mean(axis=1)
            neg = pd.concat([_z(s) for s in rule["negative"]], axis=1).mean(axis=1)
            out[axis] = pos - neg

        else:
            raise ValueError(f"Unknown formula type for axis {axis}: {rule['type']}")

    return pd.DataFrame(out).T  # axis × sample


def compute_axis_zscored(axis_raw: pd.DataFrame,
                          axis_mu_sigma: Dict[str, Tuple[float, float]]
                          ) -> pd.DataFrame:
    """Final z-score per framework axis using axis-level mu/sigma from the
    calibration JSON."""
    out = {}
    for axis in axis_raw.index:
        mu, sigma = axis_mu_sigma.get(axis, (0.0, 1.0))
        out[axis] = z_score_against_mu_sigma(axis_raw.loc[axis], mu, sigma)
    return pd.DataFrame(out).T


# =============================================================================
# Tier / class assignment
# =============================================================================

def assign_quartile_tier(z_series: pd.Series, cuts: List[float],
                          labels: List[str]) -> pd.Series:
    """Assign quartile tier given pre-computed Q25 / Q50 / Q75 cuts.
    `labels` must be ordered low→high, i.e. [T4, T3, T2, T1]."""
    bins = [-np.inf] + sorted(cuts) + [np.inf]
    out = pd.cut(z_series, bins=bins, labels=labels, include_lowest=True)
    return out.astype(str)


def assign_tertile_class(z_series: pd.Series, cuts: List[float],
                          labels: List[str]) -> pd.Series:
    """Assign tertile class given pre-computed Q33 / Q66 cuts.
    `labels` ordered low→high, e.g. [FRAGILE, PARTIAL, RESILIENT]."""
    bins = [-np.inf] + sorted(cuts) + [np.inf]
    out = pd.cut(z_series, bins=bins, labels=labels, include_lowest=True)
    return out.astype(str)


# =============================================================================
# Public API
# =============================================================================

__all__ = [
    "AXIS_EC_DICT",
    "AXIS_FORMULA",
    "FRAMEWORK_AXES_14",
    "compute_sub_axis_log1p",
    "compute_axis_raw",
    "compute_axis_zscored",
    "assign_quartile_tier",
    "assign_tertile_class",
]
