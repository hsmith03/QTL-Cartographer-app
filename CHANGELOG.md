# Changelog

## Unreleased

- Keep empirical permutation thresholds isolated to the active project stem.
- Preserve native marker positions in exported map CSV files.
- Validate quoted CSV fields, row widths, duplicate IDs, and matching R/qtl2
  genotype/phenotype sample IDs.
- Require explicit PLINK genetic-position units and report VCF variants skipped
  because they are multiallelic.
- Preview privacy-sensitive reproducibility bundle contents and exclude
  unrelated file types.
- Require the supported .NET Framework 4.8 runtime, pin CI actions by commit,
  derive package versions from application metadata, and publish SHA-256 files.

## 3.0.0 - 2026-07-23

- Repair clipped banner, subtitle, navigation categories, and tab labels at
  standard and high-DPI display scales; retain the version only in the title bar.
- Add 1.5-LOD, 2-LOD, and native-bootstrap peak confidence intervals.
- Add peak-selected genotype/phenotype distributions, mean-effect plots,
  uncertainty bars, group sample sizes, and additive/dominance summaries.
- Add chart zoom, chromosome selection, marker labels, confidence bands,
  accessible colors, and multi-run comparison.
- Add experiment-wide and chromosome-wide permutation thresholds with clear
  multiple-testing explanations.
- Move permutation tests to recoverable background jobs with deterministic
  seeds, pause/resume, cancellation, progress estimates, and persisted state.
- Add CSV templates and CSV, R/qtl, R/qtl2, PLINK, and VCF interoperability.
- Add missingness, allele-frequency, segregation-distortion, and phenotype
  diagnostics plus model-aware covariate/environment/interaction guidance.
- Add manuscript methods text and reproducibility bundles with settings,
  inputs, commands, logs, versions, random seeds, checksums, and rerun scripts.
- Add golden native benchmarks, independent R/qtl2 CSV comparison, cross-type
  validation, malformed/large/multi-trait cases, and parser fuzz testing.
- Make project saves atomic, retain backup copies, and detect changed inputs
  through SHA-256 hashes.
- Add an interactive example analysis tutorial.

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
