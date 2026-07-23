# Changelog

## 2.0.0 - 2026-07-23

- Add an integrated results dashboard for native `.z` and `.eqt` output with
  chromosome LR/LOD curves, empirical thresholds, and sortable peak summaries.
- Add native `Zmapqtl` permutation controls and calculate thresholds from the
  observed distribution of experiment-wide maximum LR statistics.
- Export CSV, PNG, SVG, and self-contained HTML analysis reports.
- Add guided new-project and data-import validation for maps, crosses,
  genotypes, phenotypes, markers, chromosomes, missing values, and sample size.
- Add reusable project save/open with version, selected arguments, queue state,
  elapsed time, and results inventory.
- Add typed numeric and mapped option controls with native help ranges.
- Add workflow progress, per-stage status, elapsed time, failure recovery, and
  retry-from-failed-stage behavior.
- Add filterable project files, file-purpose descriptions, and in-app previews.
- Add overwrite protection, detailed failure remedies, keyboard shortcuts,
  accessible labels, scalable fonts, and high-DPI support.
- Split results parsing, reports, project persistence, guided import, and
  process execution out of `MainForm.cs`.

## 1.1.0 - 2026-07-23

- Preserve multi-chromosome cross imports by carrying the active project stem and
  linkage-map filename between `Rmap` and `Rcross`.
- Warn before `Rcross` uses its legacy one-chromosome fallback when a selected
  input dataset has no available linkage map.
- Synchronize the GUI filename stem and Rmap/Rcross file defaults from
  `qtlcart.rc`, with compatibility fallback to the latest QTL Cartographer log.
- Add a five-chromosome regression workflow using the original Nuzhdin example.
- Add explicit Windows application and package version metadata.

## 1.0.0 - 2026-07-22

- Initial Windows GUI and x64 port of the complete QTL Cartographer 1.17j suite.
