# How to cite TAGMOS

Use of this bundle for any publication, preprint, presentation, dataset or
software release requires citation of the primary TAGMOS publication.

This page lists the citation forms required by the Public Research Licence
(see `LICENSE.md` §4).

---

## 1. Primary citation · mandatory

> Soverini M., Lotfollahdzadeh A., di Rito L., Viciani E., Padella A.,
> Santacroce B., Marcante A., Monaldi C., Velichevskaya A., Castagnetti A.
> *Functional multi-axis decomposition of the human gut microbiome: an
> operational definition of eubiosis and dysbiosis.* bioRxiv 2026.
> doi: [TBD — populate with the bioRxiv DOI on publication]

### BibTeX

```bibtex
@article{tagmos2026,
  title   = {Functional multi-axis decomposition of the human gut microbiome:
             an operational definition of eubiosis and dysbiosis},
  author  = {Soverini, Matteo and Lotfollahdzadeh, Ashkan and di Rito, Laura
             and Viciani, Elisa and Padella, Antonella and Santacroce, Barbara
             and Marcante, Andrea and Monaldi, Cecilia and Velichevskaya, Alena
             and Castagnetti, Andrea},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {[TBD]}
}
```

## 2. Recommended additional citations · use as applicable

If you have used a specific extension of the framework, also cite the
corresponding Supplementary Note from the primary paper:

| Extension used | Cite |
|---|---|
| Cross-pipeline robustness on cMD 18,138 metagenomes | SN1 + main text §K3 |
| Crohn-vs-UC subtype discrimination | SN1 + main text §K2 |
| Wood-Ljungdahl ancestral preservation (paleofecal / Iceman / Sardinia centenarian) | SN2 + main text §K5 |
| Beta- / alpha-diversity falsification | SN4 + main text §3 |
| Adjusted models + 5-test robustness battery | SN5 + main text §Robustness |

## 3. Mandatory attribution language in Methods

You are required to include in the Methods section of the publication a
sentence equivalent to:

> *"Substrate-functional decomposition was computed using the TAGMOS
> framework (Soverini et al., bioRxiv 2026), recalibrated on the present
> cohort following the public calibration recipe (TAGMOS_PUBLIC_BUNDLE_v10_5)
> under the TAGMOS Public Research Licence v1.0."*

## 4. Naming convention for recalibrated instances

Because TAGMOS calibrates on each user's cohort, the same framework yields
numerically different cuts depending on the anchoring cohort. To avoid
confusion in the literature, please use the following naming convention.

### General template

```
TAGMOS-<cohort-tag>-<n>
```

Where:

- `<cohort-tag>` is a short identifier of the anchoring cohort
  (max 20 characters, alphanumeric + hyphen)
- `<n>` is the cohort size used for calibration

### Examples

- A recalibration on Stanford-IBD case-controls n = 850
  → `TAGMOS-StanfordIBD-850`
- A recalibration on European Multiple Sclerosis cohort n = 412
  → `TAGMOS-EurMS-412`
- A recalibration on Japanese healthy adults from Nishijima et al. n = 400
  → `TAGMOS-NishijimaJP-400`

### How to declare in your Methods

> *"We recalibrated TAGMOS on the present cohort (n = 850; IBD case-control
> sub-cohort of the Stanford cohort), producing the calibration instance
> hereafter referred to as TAGMOS-StanfordIBD-850. Frozen calibration cuts of
> this instance, together with the SHA-256 hash of the calibration JSON,
> are reported in Supplementary Table S[X]."*

When comparing your results against published findings of the primary TAGMOS
paper, refer to the Wellmicro calibration as `TAGMOS-WMP-IT-6508` (the
6,508-subject Italian RWE anchor of *Soverini et al. 2026*).

## 5. Citation requirement for derivative work

If you publish a *modification* of the framework — new axis, new EC entry,
re-engineered formula or analogous extension — you must include the
following sentence:

> *"This work extends/modifies the TAGMOS framework (Soverini et al.,
> bioRxiv 2026). The base architecture, EC dictionary and calibration recipe
> derive from the public release of TAGMOS v4.5.13 under the Public Research
> Licence v1.0. Modifications introduced in the present work are described
> in §[…] below."*

## 6. Data availability statement

If you publish a paper that uses TAGMOS, please include a Data Availability
statement equivalent to:

> *"Substrate-functional classification of the cohort was performed with
> the TAGMOS framework v4.5.13 public release bundle, available at
> https://github.com/<TBD>/TAGMOS under the TAGMOS Public Research Licence v1.0.
> The user-side calibration JSON and per-sample classification output
> generated for the present analysis are deposited at [your repository link]."*

## 7. Failure to cite

Failure to provide the citation as specified above constitutes a breach
of the Public Research Licence (LICENSE.md §4) and terminates your right
to use the Software under §8.

---

*Citation guide v10.5 · 2026-06-06*
