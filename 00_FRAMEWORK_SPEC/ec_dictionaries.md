# EC dictionaries · 151 bottleneck EC entries

This file lists the Level-4 EC numbers used in the 14 framework axes. The
EC dictionaries are the same in **Version A** (single-output) and
**Version B** (with composite indices).

## Selection criteria

Each EC entry was selected against four criteria:

1. **Biochemical bottleneck** — the EC catalyses a step that cannot be
   by-passed within its pathway.
2. **Evolutionary conservation** — the EC corresponds to chemistry
   inherited from anaerobic ancestors and conserved across phylogenetically
   distant gut taxa.
3. **Direct mappability** to Level-4 EC numbers from shotgun functional
   profilers (WMP, HUMAnN3, Meteor2) with controlled pipeline-conditional
   variance.
4. **Pathway specificity and clinical interpretability**.

---

## H_sink_eubiotic group

### H_sink_acetogenesis_WL · Wood-Ljungdahl reductive acetyl-CoA pathway

| EC | Enzyme | Role |
|---|---|---|
| 1.2.7.4 | anaerobic CO dehydrogenase (CODH, Ni-Fe) | C1 reduction in the carbonyl branch (patent-stringent diagnostic) |
| 2.3.1.169 | CO-methylating acetyl-CoA synthase (ACS, the bifunctional CODH/ACS) | last committed step of acetate synthesis (patent-stringent diagnostic) |
| 6.3.4.3 | formate-tetrahydrofolate ligase (FtfL) | methyl branch initiation |
| 1.5.1.5 | methylene-THF dehydrogenase | methyl branch |
| 3.5.4.9 | methenyl-THF cyclohydrolase | methyl branch |
| 1.5.1.20 | methylene-THF reductase (MetF) | methyl branch closure |

### H_sink_methanogenesis · hydrogenotrophic methanogenesis (archaea)

| EC | Enzyme | Role |
|---|---|---|
| 2.8.4.1 | methyl-coenzyme M reductase (Mcr) | final methanogenesis step |
| 1.5.98.1 | F420-dependent methylene-H4MPT dehydrogenase | C1 reduction in methanogenic route |
| 1.5.98.2 | F420-dependent methylene-H4MPT reductase | C1 reduction |
| 1.5.98.3 | methylene-H4MPT reductase | C1 reduction |
| 1.8.7.3 | F420-reducing hydrogenase | H₂ activation in methanogens |

### H_sink_sulfate_reduction · dissimilatory sulfate reduction

| EC | Enzyme | Role |
|---|---|---|
| 1.8.99.5 | dissimilatory sulfite reductase (DsrAB) | sulfite → sulfide reduction |
| 1.8.4.8 | adenylyl-sulfate reductase (Apr) | APS → AMP + sulfite |

---

## H_sink_dysbiotic group

### Aerobic respiration (O₂ as terminal electron acceptor)

| EC | Enzyme | Role |
|---|---|---|
| 7.1.1.9 | cytochrome c oxidase | O₂ reduction |
| 7.1.1.7 | quinol oxidase (cytochrome bd) | O₂ reduction |

### Denitrification (NO₃⁻ → N₂)

| EC | Enzyme | Role |
|---|---|---|
| 1.7.1.15 | nitrate reductase (NarG/NapA) | NO₃⁻ → NO₂⁻ |
| 1.7.2.2 | nitrite reductase (NirS/NirK) | NO₂⁻ → NO |
| 1.7.5.1 | nitric oxide reductase (Nor) | NO → N₂O |
| 1.7.2.1 | nitrite reductase NO-forming | NO₂⁻ → NO |

### DNRA (dissimilatory nitrate reduction to ammonium)

| EC | Enzyme | Role |
|---|---|---|
| 1.7.1.4 | nitrite reductase (Nrf, NH₄⁺-forming) | NO₂⁻ → NH₄⁺ |
| 1.9.3.1 | cytochrome c nitrite reductase | NO₂⁻ → NH₄⁺ |

---

## BA · bile-acid channel

### BA_BSH · bile-salt hydrolase (deconjugation)

| EC | Enzyme |
|---|---|
| 3.5.1.24 | bile-salt hydrolase (BSH) |

### BA_dehydroxylation · 7α-dehydroxylation

| EC | Enzyme |
|---|---|
| 1.3.1.115 | 7α-hydroxysteroid dehydrogenase (BaiCD) |
| 1.17.99.1 | 7α-dehydroxylase associated steps |
| 1.1.1.391 | 12α-hydroxysteroid dehydrogenase |

### BA_epimerization · primary→secondary epimerization

| EC | Enzyme |
|---|---|
| 1.1.1.50 | 3α-hydroxysteroid dehydrogenase |
| 1.1.1.52 | 3α-hydroxysteroid dehydrogenase (NADH) |

### BA_conjugation · re-acylation

| EC | Enzyme |
|---|---|
| 6.2.1.7 | cholate-CoA ligase |
| 2.3.1.65 | bile acid CoA:amino acid N-acyltransferase |

---

## TRP · tryptophan channel

### Trp_indole · indole / IAA route

| EC | Enzyme |
|---|---|
| 4.1.99.1 | tryptophanase (indole-producing) |
| 1.4.3.2 | aromatic amino acid oxidase |
| 4.1.1.74 | indolepyruvate decarboxylase |

### Trp_kynurenine · kynurenine route

| EC | Enzyme |
|---|---|
| 1.13.11.11 | tryptophan 2,3-dioxygenase |
| 1.13.11.52 | indoleamine 2,3-dioxygenase |
| 3.5.1.9 | kynureninase |

---

## PROT · protein-derived channel

| EC | Enzyme |
|---|---|
| 3.4.21.- | trypsin-like proteases |
| 3.4.24.- | metalloproteases |
| 4.1.1.25 | tyrosine decarboxylase |
| 4.1.1.18 | lysine decarboxylase |
| 4.1.1.17 | ornithine decarboxylase |

---

## MUC · mucin-utilisation channel

### Mucin_specialised · MUC-glycan-targeted

| EC | Enzyme |
|---|---|
| 3.2.1.49 | α-N-acetylgalactosaminidase |
| 3.2.1.97 | endo-α-N-acetylgalactosaminidase |
| 3.2.1.96 | mannosyl-glycoprotein endo-β-N-acetylglucosaminidase |
| 3.2.1.51 | α-L-fucosidase |
| 3.2.1.18 | sialidase (neuraminidase) |

### Mucin_general · cross-substrate

| EC | Enzyme |
|---|---|
| 3.2.1.- | broader glycoside hydrolase family |

---

## HIS · histamine channel

### Histamine_HDC · production

| EC | Enzyme |
|---|---|
| 4.1.1.22 | histidine decarboxylase (HDC) |

### Histamine_HAL + Histamine_NAT · degradation

| EC | Enzyme |
|---|---|
| 4.3.1.3 | histidine ammonia-lyase (HAL) |
| 2.3.1.5 | arylamine N-acetyltransferase (NAT) |

---

## B12 · cobalamin sub-axis

### B12_biosynthesis_de_novo

| EC | Enzyme |
|---|---|
| 2.1.1.151 | cobalt-precorrin methylase (CobI) |
| 4.99.1.3 | sirohydrochlorin cobaltochelatase (CbiK) |

### B12_last_step_adenosyl

| EC | Enzyme |
|---|---|
| 2.5.1.17 | cob(I)alamin adenosyltransferase (CobO, BtuR) |

### B12_cobalt_insertion

| EC | Enzyme |
|---|---|
| 6.6.1.2 | cobaltochelatase (CobN/CobS/CobT) |

### B12_DMB_ligand

| EC | Enzyme |
|---|---|
| 2.4.2.21 | nicotinate-nucleotide DMB phosphoribosyltransferase (CobT) |

### B12_transporter + B12_consumer + B12_indep_alt

| EC | Function |
|---|---|
| BtuB / BtuF | cobalamin uptake |
| 5.4.99.2 | methylmalonyl-CoA mutase (consumer) |
| 2.1.1.13 | methionine synthase (consumer, B12-dependent) |
| 2.1.1.10 | homocysteine methyltransferase (B12-independent alternative) |

---

## K2 · menaquinone sub-axis

### K2_pathway_early

| EC | Enzyme |
|---|---|
| 4.1.1.71 | isochorismate synthase (MenF) |
| 6.2.1.26 | O-succinylbenzoate-CoA ligase (MenE) |
| 4.1.2.62 | SHCHC synthase (MenD) |

### K2_pathway_late

| EC | Enzyme |
|---|---|
| 4.1.3.36 | DHNA-CoA synthase (MenB) |
| 3.1.2.28 | DHNA-CoA thioesterase (MenI) |

### K2_last_step

| EC | Enzyme |
|---|---|
| 2.1.1.163 | demethylmenaquinone methyltransferase (MenG) |

### K2_utilization

| EC | Enzyme |
|---|---|
| 1.6.5.- | NADH dehydrogenase (K2 reduction) |

---

## TMA · methylamine channel

### TMA_production

| EC | Enzyme |
|---|---|
| 4.3.99.4 | carnitine TMA-lyase (CntA) |
| 4.3.99.3 | choline TMA-lyase (CutC) |
| 1.7.99.1 | TMAO reductase (TorC) |
| 4.3.99.1 | glycine betaine TMA-lyase (YeaY) |

### TMA_degradation

| EC | Enzyme |
|---|---|
| 1.14.13.148 | TMA monooxygenase (TMM) |
| 2.1.1.157 | trimethylamine methyltransferase (MttB) |

---

## Polyamine_putrescine_branch (POL sub-axis input)

| EC | Enzyme |
|---|---|
| 4.1.1.17 | ornithine decarboxylase (OdcD) |
| 4.1.1.19 | arginine decarboxylase (AdiA / SpeA) |
| 3.5.3.11 | agmatinase (SpeB) |

---

## Other sub-axes — short list

| Sub-axis | Key ECs |
|---|---|
| **GABA** | EC 4.1.1.15 (glutamate decarboxylase, GAD) |
| **BUT** | EC 2.8.3.8 (butyryl-CoA transferase) + species-marker panel (TAX-first) |
| **PRP** | EC 4.2.1.28 (propanediol dehydratase) + EC 1.1.1.202 |
| **LAC** | EC 5.1.2.1 (lactate racemase, LarA) + EC 5.1.2.- (LDH stereo) |
| **UREM** | EC 4.1.99.- (tyrosine ammonia-lyase derivatives, p-cresol pathway) |
| **IRON** | EC 6.3.2.- (siderophore biosynthesis ligases) |
| **LPS** | EC 2.3.1.- (LpxA, LpxD, LpxM hexa-acylation) |
| **ETU** | EC 4.3.1.7 (ethanolamine ammonia-lyase, EutBC) |
| **HYS** | EC 4.2.1.22 (cystathionine β-synthase) + EC 2.8.1.- (sulfide-yielding routes) |

---

## Total inventory

| Group | n EC entries |
|---|---:|
| H_sink_eubiotic (WL + methano + sulfate) | 13 |
| H_sink_dysbiotic | 8 |
| BA (4 sub-axes) | ~ 12 |
| TRP (2 sub-axes) | ~ 8 |
| PROT | ~ 6 |
| MUC (2 sub-axes) | ~ 10 |
| HIS (3 sub-axes) | ~ 4 |
| B12 (7 sub-axes) | ~ 15 |
| K2 (4 sub-axes) | ~ 10 |
| POL | ~ 4 |
| GABA | 1 |
| TMA (2 sub-axes) | ~ 6 |
| BUT + PRP + LAC | ~ 8 |
| UREM, IRON, LPS, ETU, HYS | ~ 20 |
| **Total bottleneck EC** | **~ 125 high-priority + ~ 26 sub-priority = 151** |

The exact per-axis EC list is implemented in `classifier_single.py` (and
`classifier_composite.py`) as a Python dictionary `AXIS_EC_DICT` and is
the single source of truth at runtime. The list above is the human-readable
documentation; the runtime list is the authoritative one.

---

## Citing the EC list

When you publish results based on TAGMOS, you do not need to list each EC
individually in your Methods. It is sufficient to cite TAGMOS as required by
`HOW_TO_CITE.md` and reference this document (and the runtime
`AXIS_EC_DICT`) by version (v10.5).

---

*EC dictionaries v10.5 · 2026-06-06*
