# QTL Cartographer for Windows

Current Windows GUI version: **3.0.0** (native analysis engine: 1.17j).

QTL Cartographer for Windows combines the complete QTL Cartographer 1.17 C
analysis suite with a native desktop GUI. The statistical engine and file
formats remain compatible with the Unix release; the GUI adds discoverable
option forms, project folders, live output, cancellable runs, workflow queues,
sample data, and generated-file management.

## Run the packaged application

1. Download and extract `QTL-Cartographer-v3.0.0-Windows-x64.zip`.
2. Run `QTL-Cartographer.exe`.
3. Choose an analysis tool, select a working directory, enable the desired
   options, and click **Run analysis**.

The **Workflow queue** can execute multiple tools in sequence. **Load example
analysis** creates and configures the original nine-stage sample analysis.
Every original command-line executable is also included under `tools`.

## GUI capabilities

- All 16 analysis and data-preparation programs
- Results dashboard with chromosome LR/LOD curves and highlighted threshold
- 1.5-LOD, 2-LOD, and bootstrap confidence intervals around QTL peaks
- Peak-linked genotype/phenotype distributions, effect plots, uncertainty bars,
  additive/dominance effects, and genotype-group sample sizes
- Chromosome selection, zoom, marker labels, confidence bands, accessible
  color palettes, and side-by-side run/model/trait comparison
- Experiment-wide, chromosome-wide, and manual multiple-testing choices with
  interpretation guidance
- Background permutation jobs with deterministic seeds, progress estimates,
  pause/resume, cancellation, and native partial-result recovery
- Sortable QTL peak table with positions, markers, effects, and significance
- Native empirical permutation workflow with configurable count and alpha
- CSV peak tables, publication-quality PNG/SVG charts, and self-contained HTML reports
- Guided new-project/import wizard for example or user data
- CSV templates and import/export support for validated wide genotype matrices,
  simplified R/qtl2-style CSV matrices, PLINK PED/MAP with explicit genetic-map
  units, and biallelic VCF `GT` data
- Missingness heatmaps, allele frequencies, segregation-distortion flags, and
  phenotype histograms before analysis
- Model-aware covariate, environment, genotype-by-environment, multi-trait, and
  population-structure design guidance for native models that support them
- Preflight checks for chromosomes, markers, sample size, phenotypes, cross type,
  missing values, and chromosome assignments
- Reusable `.qtlproject` files containing inputs, selected arguments, workflow
  state, version, and result inventory
- Typed numeric/dropdown option forms generated from native help ranges
- Input/output file pickers and command preview
- Automatic and non-verbose execution modes
- Live stdout/stderr console with cancellation, elapsed time, and log export
- Recoverable workflow queue with pending/running/succeeded/failed states,
  progress, preserved completed stages, and retry from failure
- Filterable project file browser with file-purpose explanations and text preview
- Overwrite warnings and detailed command/error/remedy dialogs
- Keyboard shortcuts, high-DPI scaling, and accessible control names
- Manuscript-ready methods paragraphs and reproducibility ZIP bundles containing
  inputs, settings, seeds, versions, SHA-256 checksums, logs, and rerun scripts
- Golden-result benchmarks for bundled examples plus independent reference-peak
  CSV comparison support
- Atomic project saves, backup copies, and changed-input hash detection
- Interactive example-data tutorial
- Original CLI interfaces for scripts and reproducible pipelines
- Automatic filename-stem synchronization from project resource/log files
- Missing-map warning that prevents accidental one-chromosome cross imports

## Build

Requirements:

- Windows 11, or a Windows 10 installation covered by LTSC or ESU support
- Visual Studio 2019 or later Build Tools
- **Desktop development with C++** workload
- .NET Framework 4.8 runtime and reference assemblies

From PowerShell:

```powershell
.\build.ps1
.\test.ps1
```

The application is created under `build\app`; the distributable ZIP is
`build\QTL-Cartographer-v3.0.0-Windows-x64.zip`.

## Scientific interpretation

The dashboard reads the native `*.z` likelihood-ratio profile and `*.eqt` peak
estimates. LOD values are calculated as `LR / (2 ln 10)`. When a permutation
test is available, the displayed empirical threshold is the selected quantile
of the global maximum LR distribution. Peaks should be interpreted against a
validated threshold and an appropriate experimental model.

Support intervals are calculated by tracing the native profile down 1.5 or 2
LOD units from each peak. Bootstrap percentile intervals are shown when native
bootstrap output is present. Interoperability imports preserve the source data,
but phenotype and genetic-map files remain required when formats such as VCF
provide genotypes alone.

R/qtl's native combined/rotated CSV layouts and R/qtl2 control-file semantics
are not silently treated as generic matrices. Convert them to the documented
wide CSV templates first so genotype encodings and sample identifiers can be
validated explicitly.

## License and attribution

This repository is licensed under the GNU General Public License v3. The QTL
Cartographer analysis engine was created by C. J. Basten, B. S. Weir, and
Z.-B. Zeng. Its upstream source is available at
[cbasten/qtlcart](https://github.com/cbasten/qtlcart).
