# The panel file

The analysis scripts in `../analysis/` do not contain the definition of any
metabolic axis. They read it from a JSON file passed with `--panel`, so that the
same code can be run on any set of enzyme groupings.

The panel used in the paper is not distributed. It is the object of pending
patent applications and is available to editors and referees on request, as
stated in the paper's Data availability section. `example_panel.json` in this
directory is a placeholder that has nothing to do with it: it groups a handful
of well-known housekeeping and central-metabolism enzymes, chosen only so that
the pipeline has something to run on.

## Format

A JSON object. Each key is an axis name; each value describes how that axis is
scored and how it should be read.

```json
{
  "MY_AXIS": {
    "ecs": ["1.1.1.1", "1.2.1.3"],
    "valence": "protective"
  },
  "MY_RATIO_AXIS": {
    "num": ["2.7.1.11"],
    "den": ["4.1.2.13"],
    "valence": "danger"
  }
}
```

`ecs`
: enzymes scored together. Each enzyme is z-scored across samples and the axis
  is the mean of those z-scores. This is the form the analysis scripts consume.

`num` / `den`
: an axis expressed as a log-ratio of two enzyme groups, the form used for a
  balance between two competing routes. Axes given only as `num`/`den` are
  skipped by the scripts here, which expect an `ecs` list.

`valence`
: one of `protective`, `danger` or `context`. It sets the direction in which
  the dysbiotic tail is defined: for a `danger` axis the tail is the upper one
  (z ≥ +1 against the study's own controls), for a `protective` axis the lower
  one. `context` axes are reported but not used in the transdiagnostic
  ordering. Omitting it is allowed and treated as `context`.

Enzyme identifiers are EC numbers as they appear in the columns of the EC
matrix. Any enzyme not present in the matrix is dropped, with a count reported
on stderr; an axis that loses all of its enzymes is skipped.

## Building your own

Nothing in the code assumes the panel is ours. To reproduce the random-panel
control of Supplementary Table 2 on your own terms, generate panels matched to
your curated one on size and detection prevalence, write each as a panel file,
and run `analysis/channel_evidence_grade.py` or `analysis/threshold_vs_mean.py`
over them: the null distribution the paper reports is the distribution of those
runs.
