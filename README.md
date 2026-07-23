# QTL Cartographer for Windows

Current Windows GUI version: **2.0.0** (native analysis engine: 1.17j).

QTL Cartographer for Windows combines the complete QTL Cartographer 1.17 C
analysis suite with a native desktop GUI. The statistical engine and file
formats remain compatible with the Unix release; the GUI adds discoverable
option forms, project folders, live output, cancellable runs, workflow queues,
sample data, and generated-file management.

## Run the packaged application

1. Download and extract `QTL-Cartographer-v2.0.0-Windows-x64.zip`.
2. Run `QTL-Cartographer.exe`.
3. Choose an analysis tool, select a working directory, enable the desired
   options, and click **Run analysis**.

The **Workflow queue** can execute multiple tools in sequence. **Load sample
workflow** creates and configures the original nine-stage sample analysis.
Every original command-line executable is also included under `tools`.

## GUI capabilities

- All 16 analysis and data-preparation programs
- Results dashboard with chromosome LR/LOD curves and highlighted threshold
- Sortable QTL peak table with positions, markers, effects, and significance
- Native empirical permutation workflow with configurable count and alpha
- CSV peak tables, publication-quality PNG/SVG charts, and self-contained HTML reports
- Guided new-project/import wizard for example or user data
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
- Original CLI interfaces for scripts and reproducible pipelines
- Automatic filename-stem synchronization from project resource/log files
- Missing-map warning that prevents accidental one-chromosome cross imports

## Build

Requirements:

- Windows 10 or later
- Visual Studio 2019 or later Build Tools
- **Desktop development with C++** workload
- .NET Framework 4 reference assemblies

From PowerShell:

```powershell
.\build.ps1
.\test.ps1
```

The application is created under `build\app`; the distributable ZIP is
`build\QTL-Cartographer-v2.0.0-Windows-x64.zip`.

## Scientific interpretation

The dashboard reads the native `*.z` likelihood-ratio profile and `*.eqt` peak
estimates. LOD values are calculated as `LR / (2 ln 10)`. When a permutation
test is available, the displayed empirical threshold is the selected quantile
of the global maximum LR distribution. Peaks should be interpreted against a
validated threshold and an appropriate experimental model.

## License and attribution

This repository is licensed under the GNU General Public License v3. The QTL
Cartographer analysis engine was created by C. J. Basten, B. S. Weir, and
Z.-B. Zeng. Its upstream source is available at
[cbasten/qtlcart](https://github.com/cbasten/qtlcart).
