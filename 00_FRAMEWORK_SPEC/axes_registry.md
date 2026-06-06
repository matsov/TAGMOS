# TAGMOS axes registry · 14 single-output framework axes

This is the formal definition of the 14 framework axes shipped with the
public release. Each axis is a single output (one z-score per sample) and
contributes either to the 36-cell state-space partition (PRIMARY +
CO-PRIMARY) or as an independent host-interface channel.

---

## Hierarchy summary

| # | Layer | Axis | Composition | Scoring | Direction of concern |
|---:|---|---|---|---|---|
| 1 | PRIMARY | **Engine z** | terminal hydrogen sinks: 3 eubiotic (WL, methano, sulfate) vs 3 dysbiotic (aerobic O₂, denitrification, DNRA) | EC | negative pole = dysbiotic |
| 2 | CO-PRIMARY | **SYN_cof** | aggregator of B12 + K2 + POL + GABA cooperative biosynthesis | EC | negative pole = cofactor fragility |
| 3 | CO-PRIMARY | **SYN_carb** | aggregator of BUT + PRP + LAC saccharolytic cooperative metabolism | MIXED (TAX-first BUT) | negative pole = saccharolytic fragility |
| 4 | CHANNEL | **BA** | bile-acid: BSH + 7α-dehydroxylation + epimerization + conjugation | EC | positive pole = elevated BA |
| 5 | CHANNEL | **TRP** | tryptophan-derived: indole route + kynurenine route + tryptamine + IAA | EC | positive pole = elevated route |
| 6 | CHANNEL | **PROT** | protein-derived metabolism | EC | positive pole = elevated |
| 7 | CHANNEL | **MUC** | mucin-utilisation: specialised + general | EC | positive pole = elevated MUC |
| 8 | CHANNEL | **HIS** | histamine net: HDC − (HAL + NAT) | EC | positive pole = net production |
| 9 | CHANNEL | **TMA** | methylamine routes: TMA_production − TMA_degradation | EC | positive pole = net production |
| 10 | CHANNEL | **UREM** | uremic / phenolic precursors (p-cresol, indoxyl) | EC | positive pole = elevated precursors |
| 11 | CHANNEL | **IRON** | iron-uptake / siderophore biosynthesis | EC | positive pole = elevated siderophore |
| 12 | CHANNEL | **LPS** | lipopolysaccharide signalling (hexa-acylated LpxLM) | EC | positive pole = pro-inflammatory LPS |
| 13 | CHANNEL | **ETU** | ethanolamine utilisation | EC | positive pole = elevated ETU |
| 14 | CHANNEL | **HYS** | hydrogen-sulfide (H₂S) handling | EC | positive pole = elevated H₂S |

---

## How each axis is computed

The general recipe is identical for all axes:

1. **Per-EC log1p z-score** — for each EC in the axis dictionary, compute
   `log1p(count_per_sample)` then z-score against the per-EC mean / sigma
   computed on the calibration cohort.
2. **Axis raw score** — sum or mean of the per-EC z-scores within the axis
   dictionary. For directional axes (e.g. HIS net production = HDC − HAL − NAT),
   the formula is implemented as a directional combination of the relevant
   sub-dictionaries.
3. **Axis z-score** — z-score of the axis raw score against the axis-level
   mean / sigma computed on the calibration cohort.
4. **Axis class / tier** — quartile partition for Engine z (T1_EUBIOTIC →
   T4_DYSBIOTIC) or tertile partition for syntrophy aggregators and channels
   (FRAGILE / PARTIAL / RESILIENT for SYN axes; LOW / MID / HIGH for channels).

All four steps are deterministic given the calibration JSON. The calibration
JSON is computed once on the user's cohort by
`calibrate_local_single.py` (Version A) or `calibrate_local_composite.py`
(Version B) and re-used for every subsequent classification.

---

## Engine z · the PRIMARY axis

### Biological rationale

Engine z is the normalised balance between counts assigned to **eubiotic
terminal hydrogen sinks** and **dysbiotic alternative-electron-acceptor sinks**.
Hydrogen is the conserved redox currency of the gut microbial consortium;
its terminal disposal balance therefore is a natural candidate for the
primary organising axis. The formal expression is:

```
H_eubiotic_pct  =  100 × Σ_eubiotic / (Σ_eubiotic + Σ_dysbiotic)
Engine z         =  (H_eubiotic_pct − μ_calib) / σ_calib
```

Where Σ_eubiotic and Σ_dysbiotic are the cohort-level log1p-aggregates of
the eubiotic and dysbiotic terminal-sink EC dictionaries respectively.

### Eubiotic terminal sinks · 3 sub-axes

| Sub-axis | Diagnostic ECs | Pathway |
|---|---|---|
| Wood-Ljungdahl acetogenesis (WL) | EC 1.2.7.4 (Ni-Fe CO-dehydrogenase) + EC 2.3.1.169 (CO-methylating acetyl-CoA synthase) | reductive acetyl-CoA, autotrophic acetogenesis |
| Hydrogenotrophic methanogenesis | EC 2.8.4.1 (methyl-CoM reductase) + EC 1.5.98.1/2/3 + EC 1.8.7.3 | archaeal methanogenesis |
| Sulfate reduction | EC 1.8.99.5 + EC 1.8.4.8 | dissimilatory sulfate reduction |

### Dysbiotic alternative-electron-acceptor sinks · 1 sub-axis aggregating 3 pathways

| Sub-pathway | Diagnostic ECs |
|---|---|
| Aerobic respiration | EC 7.1.1.9 (cytochrome c oxidase) + EC 7.1.1.7 |
| Denitrification | EC 1.7.1.15 + EC 1.7.2.2 + EC 1.7.5.1 + EC 1.7.2.1 |
| DNRA (dissimilatory nitrate reduction to ammonium) | EC 1.7.1.4 + EC 1.9.3.1 |

### Tier partition

Engine z is partitioned into quartiles on the calibration cohort:

| Tier | Range | Calibration prevalence | Interpretation |
|---|---|---|---|
| T1_EUBIOTIC | Engine z > Q75 | ~ 25 % | Dominant eubiotic terminal-sink activity |
| T2_PRESERVED | Q50 < z ≤ Q75 | ~ 25 % | Eubiotic-prevalent intermediate |
| T3_ALTERED | Q25 < z ≤ Q50 | ~ 25 % | Mixed eubiotic / dysbiotic activity |
| T4_DYSBIOTIC | z ≤ Q25 | ~ 25 % | Dysbiotic terminal-sink activity dominant |

The Q25 / Q50 / Q75 numerical values are computed on the user's calibration
cohort. They are not constants of the framework.

---

## SYN_cof · the cooperative-cofactor CO-PRIMARY axis

### Sub-axes (4)

| Sub-axis | Description |
|---|---|
| B12 | cobalamin biosynthesis / consumption / transport |
| K2 | menaquinone biosynthesis |
| POL | polyamine biosynthesis (putrescine, cadaverine, agmatine chains) |
| GABA | γ-aminobutyrate production via glutamate decarboxylase (GAD) |

### Formula

```
SYN_cof_z  =  z-mean of (B12_z, K2_z, POL_z, GABA_z)
```

### Tertile partition

| Class | Range | Interpretation |
|---|---|---|
| SYN_FRAGILE | SYN_cof_z ≤ Q33 | depressed cooperative-cofactor activity |
| SYN_PARTIAL | Q33 < z ≤ Q66 | intermediate |
| SYN_RESILIENT | z > Q66 | elevated cooperative-cofactor activity |

---

## SYN_carb · the saccharolytic CO-PRIMARY axis

### Sub-axes (3)

| Sub-axis | Description | Input method |
|---|---|---|
| BUT | butyrogenic activity | TAX (curated species panel + 4 racemase / propanediol ECs) |
| PRP | propanediol pathway (P. copri P-K route) | EC |
| LAC | lactate racemase / D-L stereo-balance | EC |

### BUT taxonomy-first panel

The BUT sub-axis is taxonomy-first by design because pilot EC-only butyrogenic
scoring through the butyryl-CoA pathway was unstable across functional
profilers (WMP, HUMAnN3, Meteor2) due to pipeline-specific gene-family
granularity. The species-marker panel:

- *Faecalibacterium prausnitzii* cluster (M21/2, SL3/3, L2-6 strains)
- *Roseburia intestinalis* group
- *Agathobacter rectalis* (= *Eubacterium rectale*)
- *Eubacterium ramulus*
- *Dorea longicatena*
- *Coprococcus* spp.

### Formula

```
SYN_carb_z  =  z-mean of (BUT_z, PRP_z, LAC_z)
```

### Tertile partition

| Class | Range | Interpretation |
|---|---|---|
| SYN_carb_FRAGILE | SYN_carb_z ≤ Q33 | depressed saccharolytic activity |
| SYN_carb_PARTIAL | Q33 < z ≤ Q66 | intermediate |
| SYN_carb_RESILIENT | z > Q66 | elevated saccharolytic activity |

---

## 11 channels — at a glance

Each channel axis is computed as a single z-score per sample and partitioned
into tertiles using cuts computed on the calibration cohort. The biological
direction of concern is given in the hierarchy summary table above and in
`schema_tagmos_public_v4513.json`.

For directional channels (HIS, TMA), the formula is a **directional
combination** of the relevant sub-dictionaries:

```
HIS_net_z   =   HDC_z − (HAL_z + NAT_z)
TMA_net_z   =   TMA_production_z − TMA_degradation_z
```

For aggregating channels (BA, MUC, TRP), the formula is the mean of relevant
sub-component z-scores.

For single-component channels (PROT, UREM, IRON, LPS, ETU, HYS), the formula
is the direct log1p-z-score of the dictionary aggregate.

---

## 3 schema-architecture nominees (not in the 14 framework axes)

Three additional axes are nominated in the schema architecture for future
calibration releases but are NOT counted as framework axes in this release:

- **NAD+** · longevity-relevant NAD+ biosynthesis
- **RIB** · riboflavin biosynthesis
- **FOL** · folate biosynthesis

They are not in the 36-cell state-space partition and are not used in the
present calibration recipe. They are documented for transparency and to
anchor future extensions of the framework.

---

## Schema-architecture sub-axes vs framework axes — terminology note

The TAGMOS schema architecture enumerates **24 entries** (14 framework axes +
7 aggregator sub-axes inside SYN_cof and SYN_carb + 3 nominees). The
**14 framework axes** are the *single-output* axes that emit one z-score per
sample. The 7 sub-axes inside aggregator axes (B12, K2, POL, GABA, BUT, PRP,
LAC) are not counted as separate framework axes because they enter the
analysis only as inputs to the SYN_cof / SYN_carb z-mean.

When citing TAGMOS as a 14-axis framework, you are correctly counting the
single-output axes that contribute to the 36-cell state-space partition or
host-interface channel readouts.

---

*Axes registry v10.5 · 2026-06-06*
