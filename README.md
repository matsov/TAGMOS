<p align="center">
  <img src="T%CE%9BGM%CA%98S_Logo.jpg" alt="TAGMOS" width="360">
</p>

# TAGMOS

Analysis code and source data for ***Disease crosses a thermodynamic threshold in the human gut microbiome***
(Soverini *et al.*, submitted).

This repository turns public metagenomic data into the tables the paper
reports, and holds those tables. It contains no plotting code: the tables under
`source_data/` are the numbers behind every figure panel, and are yours to draw
however you like.

```
ec_matrix/     public HUMAnN3 profiles  ->  enzyme (EC) coordinates
analysis/      axis scores              ->  the statistics the paper reports
panel/         the axis-definition format, and a placeholder panel
example/       a small synthetic dataset, so every script can be run as-is
source_data/   the aggregated tables behind each figure panel
```

## Quick start

```bash
git clone https://github.com/matsov/TAGMOS
cd TAGMOS
pip install -r requirements.txt

cd analysis
python threshold_vs_mean.py \
    --ec ../example/ec_matrix.tsv \
    --meta ../example/metadata.tsv \
    --panel ../panel/example_panel.json \
    --richness ../example/richness.tsv \
    --out /tmp/out
```

Python 3.10 or later; numpy, pandas and scipy. Nothing is downloaded at run
time. The example data are synthetic and carry no biological meaning: they are
there so that every script can be run and inspected without any private input.

### 1 · Public profiles to enzyme coordinates

`ec_matrix/build_ec_matrix.py`

curatedMetagenomicData is distributed as HUMAnN3 UniRef90 gene families, not as
EC profiles. This script inverts HUMAnN's official level-4 EC-to-UniRef90 table
into a UniRef90-to-EC dictionary and applies it, summing over HUMAnN's taxonomic
stratification and merging studies on the union of ECs. Public inputs only.

```bash
python ec_matrix/build_ec_matrix.py \
    --map map_level4ec_uniref90.txt.gz \
    --genes cmd_genefamilies/*.tsv \
    --out ec_matrix.tsv
```

### 2 · Enzyme coordinates to axis scores — *not distributed*

This is the framework: the enzyme dictionary, the assignment of enzymes to
metabolic axes, and the frozen calibration constants. It is proprietary to
Wellmicro S.r.l. and is the subject of pending patent applications. The exact
per-axis enzyme lists, and the per-sample axis scores for the samples analysed
in the paper, are available in full to editors and referees on request, as the
paper's Data availability section states.

The reason for the boundary is specific rather than defensive. Enzyme
abundances are public; per-sample axis scores paired with the accessions they
came from are not, because the two together let a third party regress one
against the other and recover the axis definitions from a few hundred samples.
Withholding this stage is what keeps that from being possible.

### 3 · Axis scores to the reported statistics

`analysis/`

Everything the paper claims methodologically lives here. None of it needs to
know how a score was computed, so all of it is in this repository, and all of it
runs on any panel you supply.

| script | function |
|---|---|
| `threshold_vs_mean.py` | for every axis and condition, the mean reading (OR per SD) against the threshold reading (OR for the dysbiotic tail), and the class of each association |
| `threshold_vs_mean_sensitivity.py` | the same across gate thresholds and FDR cutoffs, showing the threshold-only fraction is not an artefact of either |
| `channel_evidence_grade.py` | per axis and condition: within-study OR, random-effects pool, number of independent studies, evidence grade |
| `channel_condition_matrix.py` | the axis-by-condition effect matrix and the transdiagnostic ordering of conditions |
| `taxonomic_negative_control.py` | a species guild scored against 150 prevalence-matched random species panels |
| `gmhi_benchmark.py` | the Gut Microbiome Health Index put through the identical evaluation, with its own matched-random-panel null |

`gmhi_benchmark.py` needs no panel at all: it uses the published 50-species
index of Gupta *et al.* (2020) and runs end to end on public data alone.

Run any script with `--help` for its inputs. `analysis/tagmos_io.py` documents
the metadata columns each one expects; the `sampleMetadata` table of
curatedMetagenomicData has that shape.

## The panel file

Every script takes its axis definitions from a JSON file passed with `--panel`,
not from its own source, so the same code runs on any set of enzyme groupings.
`panel/README.md` documents the format. `panel/example_panel.json` is a
placeholder built from generic central-metabolism enzymes and has no relation to
the panel used in the paper.

### Reproducing the random-panel control

Supplementary Table 2 reports the curated architecture against 300
prevalence-matched random panels. Because the panel is an input rather than a
constant, that comparison is reproducible here on any curated panel: build your
panels matched on size and detection prevalence, write each as a panel file, run
`threshold_vs_mean.py` or `channel_evidence_grade.py` over them, and the null
distribution is the distribution of those runs.

## The deposited tables

`source_data/` holds the aggregated tables the paper reports, one per figure
panel, described file by file with row counts and checksums in
[`source_data/MANIFEST.md`](source_data/MANIFEST.md).

Seventeen of the nineteen are aggregates over studies, conditions, cohort phases
or populations. Two panels of the paper are per-sample by construction, and for
those the values are deposited with every identifying column removed, rows
shuffled and values rounded, so that they redraw the panel but cannot be joined
to any sample. `source_data/deidentify.py` is the exact transformation, kept as
a record.

The only enzyme identifiers anywhere in this repository are ones the paper names
in full: EC 4.1.1.15 and EC 2.6.1.19, the two enzymes of the GABA channel, and
EC 3.5.2.6, the carbapenemase. No file maps an enzyme to an axis.

The primary metagenomes are not redistributed. Every cohort is public and its
accession is listed in Supplementary Table 13; curatedMetagenomicData is
available through its Bioconductor package.

## Licence

Code: MIT ([`LICENSE`](LICENSE)).
Data in `source_data/`: CC BY 4.0 ([`LICENSE-DATA.md`](LICENSE-DATA.md)).

The licences cover this repository only. They do not extend to the TAGMOS
framework, the enzyme dictionary, the per-axis enzyme lists or the calibration
constants, none of which are contained here.

## Citation

Soverini M., Lotfollahzadeh A., di Rito L., Padella A., Santacroce B.,
Marcante A., Monaldi C., Velichevskaya A., Viciani E., Castagnetti A.
*Disease crosses a thermodynamic threshold in the human gut microbiome.*
Submitted, 2026.

Machine-readable metadata in [`CITATION.cff`](CITATION.cff).

## Contact

Correspondence: Andrea Castagnetti, andrea.castagnetti@wellmicro.com
Requests for academic collaboration can be directed to the corresponding author.
