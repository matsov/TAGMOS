# Validation toolkit

Three scripts to validate your cohort and your calibration:

| Script | Purpose | Published benchmark |
|---|---|---|
| `beta_diversity_falsification.py` | Quantify TAGMOS 36-cell partition's variance advantage vs each outcome | Italian-RWE: 38× to 94× advantage |
| `dimensional_saturation.py` | Cramér V* curve over 1D / 2D / 3D | Italian-RWE: 3D = saturation point |
| `input_method_routing.py` | Per-sub-axis EC coverage + Rule R9 routing | per-axis ≥ 60 % EC coverage recommended |

Run them after calibrating + classifying. Their outputs let you verify that
your cohort is producing TAGMOS-typical results, before drawing
disease-association conclusions.

---

## Citation

These tools are templates for reproducing the validation analyses described
in paper 1 Methods §Architectural saturation testing, §Beta-diversity
comparison and §Per-axis input-method routing. Cite TAGMOS as required by
`../HOW_TO_CITE.md`.

---

*Validation toolkit v10.5 · 2026-06-06*
