# QTL Cartographer for Windows

Current Windows GUI version: **1.1.0** (native analysis engine: 1.17j).

QTL Cartographer for Windows combines the complete QTL Cartographer 1.17 C
analysis suite with a native desktop GUI. The statistical engine and file
formats remain compatible with the Unix release; the GUI adds discoverable
option forms, project folders, live output, cancellable runs, workflow queues,
sample data, and generated-file management.

## Run the packaged application

1. Download and extract `QTL-Cartographer-v1.1.0-Windows-x64.zip`.
2. Run `QTL-Cartographer.exe`.
3. Choose an analysis tool, select a working directory, enable the desired
   options, and click **Run analysis**.

The **Workflow queue** can execute multiple tools in sequence. **Load sample
workflow** creates and configures the original nine-stage sample analysis.
Every original command-line executable is also included under `tools`.

## GUI capabilities

- All 16 analysis and data-preparation programs
- Option forms generated directly from each program's native help output
- Input/output file pickers and command preview
- Automatic and non-verbose execution modes
- Live stdout/stderr console with cancellation and log export
- Reusable multi-program workflow queue
- Project file browser and packaged sample data
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
`build\QTL-Cartographer-v1.1.0-Windows-x64.zip`.

## License and attribution

This repository is licensed under the GNU General Public License v3. The QTL
Cartographer analysis engine was created by C. J. Basten, B. S. Weir, and
Z.-B. Zeng. Its upstream source is available at
[cbasten/qtlcart](https://github.com/cbasten/qtlcart).
