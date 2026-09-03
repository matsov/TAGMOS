"""Provenance record for the two per-sample source-data files in this repository.

Every other file under source_data/ is a study-, condition- or population-level
aggregate and is deposited verbatim. Two figure panels are per-sample by
construction and cannot be redrawn from aggregates:

    Fig. 4D              covariate-adjusted EUBIO and GABA levels, PRIMM (n = 104)
    Supplementary Fig. 1 frozen vs cohort-internal EUBIO, Parkinson cohort (n = 720)

For those two, the values are deposited without any key that would link a score
to a public sequencing run. This script is the exact transformation that was
applied to the internal tables; it is included so that the derivation is on the
record, and is not needed to reproduce the figures.

What is removed and why
-----------------------
* sample_id / run / subject : a (accession, score) pair lets a third party
  regress the published score against the public enzyme-abundance matrix and
  recover the axis definition. Without the key the values cannot be joined to
  any sample.
* logdepth (Supplementary Fig. 1) : not used by the panel, and a quasi-identifier
  -- per-sample sequencing depth is recomputable from the public FASTQ files, so
  it would re-establish the join that dropping sample_id removes.
* roster (channel evidence grades) : an internal bookkeeping tag, not plotted.

Row order is shuffled with a fixed seed and values are rounded, so no residual
ordering or precision artefact carries information about the original samples.
Within-row pairing is preserved, which is all the two panels need.
"""
import pandas as pd

SEED = 20260903

# ---- Fig. 4D : ICI (PRIMM) covariate-adjusted residuals -------------------
ici = pd.read_csv("ici_PRJEB70966_boxplot_values.tsv", sep="\t")
ici = ici.drop(columns=["sample_id"])
ici = ici.sample(frac=1.0, random_state=SEED).reset_index(drop=True)
ici[["EUBIO_resid", "GABA_resid"]] = ici[["EUBIO_resid", "GABA_resid"]].round(3)
ici.to_csv("ici/fig4D_ici_adjusted_levels.tsv", sep="\t", index=False)

# ---- Supplementary Fig. 1 : calibration invariance, PD cohort -------------
pd_ = pd.read_csv("edFig_scaleinvariance_frozen_vs_internalz_PD.csv")
pd_ = pd_.drop(columns=["sample_id", "logdepth"])
pd_ = pd_.sample(frac=1.0, random_state=SEED).reset_index(drop=True)
pd_ = pd_.round(2)
pd_.to_csv("supp_fig1/suppfig1_calibration_invariance.tsv", sep="\t", index=False)

# ---- Fig. 4B : cohort-level means only ------------------------------------
# The internal table also carries per-subject rows; the panel plots only the
# per-cohort means (the rows tagged __MEAN__), so only those are deposited.
dec = pd.read_csv("deltaEUBIO_decomposition_anchors.tsv", sep="\t")
dec = dec[dec["subject"] == "__MEAN__"].drop(columns=["subject"])
dec.to_csv("engine/fig4B_deltaEUBIO_decomposition_cohort_means.tsv", sep="\t", index=False)
