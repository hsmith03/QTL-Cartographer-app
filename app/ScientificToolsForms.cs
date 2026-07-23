using System;
using System.Collections.Generic;
using System.Drawing;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace QTLCartographer.Gui
{
    internal sealed class BenchmarkResult
    {
        public int Matched { get; set; }
        public double MaximumPositionDifference { get; set; }
        public double MaximumScoreDifference { get; set; }
        public bool Passed { get; set; }
        public string Summary { get; set; }
    }

    internal static class ScientificBenchmarkSuite
    {
        public static BenchmarkResult Compare(AnalysisResults actual, string referenceCsv, double positionTolerance, double scoreTolerance)
        {
            List<Tuple<int, double, double>> expected = File.ReadLines(referenceCsv).Skip(1)
                .Where(line => !string.IsNullOrWhiteSpace(line))
                .Select(line => line.Split(','))
                .Where(p => p.Length >= 3)
                .Select(p => Tuple.Create(int.Parse(p[0]), double.Parse(p[1], CultureInfo.InvariantCulture), double.Parse(p[2], CultureInfo.InvariantCulture)))
                .ToList();
            BenchmarkResult result = new BenchmarkResult();
            foreach (Tuple<int, double, double> item in expected)
            {
                QtlPeak match = actual.Peaks.Where(p => p.Chromosome == item.Item1)
                    .OrderBy(p => Math.Abs(p.PositionCm - item.Item2)).FirstOrDefault();
                if (match == null) continue;
                result.Matched++;
                result.MaximumPositionDifference = Math.Max(result.MaximumPositionDifference, Math.Abs(match.PositionCm - item.Item2));
                result.MaximumScoreDifference = Math.Max(result.MaximumScoreDifference, Math.Abs(match.LikelihoodRatio - item.Item3));
            }
            result.Passed = result.Matched == expected.Count &&
                result.MaximumPositionDifference <= positionTolerance && result.MaximumScoreDifference <= scoreTolerance;
            result.Summary = result.Matched + "/" + expected.Count + " reference peaks matched; maximum position difference " +
                result.MaximumPositionDifference.ToString("0.####") + " cM; maximum LR difference " +
                result.MaximumScoreDifference.ToString("0.####") + ".";
            return result;
        }
    }

    internal sealed class BenchmarkForm : Form
    {
        private readonly string directory;
        private readonly string stem;
        private TextBox output;

        public BenchmarkForm(string directory, string stem)
        {
            this.directory = directory;
            this.stem = stem;
            Text = "Scientific benchmark suite";
            Size = new Size(820, 560);
            StartPosition = FormStartPosition.CenterParent;
            FlowLayoutPanel tools = new FlowLayoutPanel { Dock = DockStyle.Top, Height = 44 };
            Button bundled = new Button { Text = "Run bundled golden benchmark", AutoSize = true };
            Button independent = new Button { Text = "Compare R/qtl2 CSV...", AutoSize = true };
            bundled.Click += delegate { RunBundled(); };
            independent.Click += delegate { RunIndependent(); };
            tools.Controls.Add(bundled); tools.Controls.Add(independent);
            output = new TextBox { Dock = DockStyle.Fill, Multiline = true, ReadOnly = true, ScrollBars = ScrollBars.Vertical, Font = new Font("Consolas", 10F) };
            Controls.Add(output); Controls.Add(tools);
        }

        private void RunBundled()
        {
            string reference = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "example", "sample-golden-peaks.csv");
            ShowResult("Bundled published QTL Cartographer example", ScientificBenchmarkSuite.Compare(ResultsParser.LoadProject(directory, stem), reference, .02, .002));
        }

        private void RunIndependent()
        {
            using (OpenFileDialog dialog = new OpenFileDialog { Filter = "R/qtl2 reference CSV (*.csv)|*.csv", InitialDirectory = directory })
                if (dialog.ShowDialog(this) == DialogResult.OK)
                    ShowResult("Independent implementation comparison", ScientificBenchmarkSuite.Compare(ResultsParser.LoadProject(directory, stem), dialog.FileName, 5.0, 3.0));
        }

        private void ShowResult(string label, BenchmarkResult result)
        {
            output.AppendText(label + "\r\n" + result.Summary + "\r\nStatus: " + (result.Passed ? "PASS" : "REVIEW") +
                "\r\n\r\nIndependent CSV columns: chromosome, position_cM, LR. Convert LOD to LR with LR = 2 ln(10) LOD when needed.\r\n\r\n");
        }
    }

    internal sealed class ModernImportForm : Form
    {
        private readonly string directory;
        private readonly string stem;
        private readonly CrossData data;
        private ComboBox format;
        private TextBox primary;
        private TextBox secondary;
        private TextBox result;

        public ModernImportForm(string directory, string stem, CrossData data)
        {
            this.directory = directory; this.stem = stem; this.data = data;
            Text = "Modern format interoperability";
            Size = new Size(820, 560); StartPosition = FormStartPosition.CenterParent; AutoScaleMode = AutoScaleMode.Dpi;
            TableLayoutPanel panel = new TableLayoutPanel { Dock = DockStyle.Fill, Padding = new Padding(18), ColumnCount = 3, RowCount = 6 };
            panel.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 150)); panel.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 100)); panel.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 95));
            format = new ComboBox { Dock = DockStyle.Fill, DropDownStyle = ComboBoxStyle.DropDownList };
            format.Items.AddRange(new object[] { "CSV", "R/qtl", "R/qtl2", "PLINK", "VCF" }); format.SelectedIndex = 0;
            panel.Controls.Add(new Label { Text = "Format", Dock = DockStyle.Fill }, 0, 0); panel.Controls.Add(format, 1, 0); panel.SetColumnSpan(format, 2);
            primary = AddPath(panel, 1, "Primary file");
            secondary = AddPath(panel, 2, "Secondary file");
            FlowLayoutPanel actions = new FlowLayoutPanel { AutoSize = true };
            Button import = new Button { Text = "Import", AutoSize = true };
            Button export = new Button { Text = "Export current project CSV", AutoSize = true };
            Button templates = new Button { Text = "Create CSV templates", AutoSize = true };
            import.Click += delegate { Import(); }; export.Click += delegate { ModernFormatConverter.ExportCrossCsv(data, directory, stem); result.Text = "Current native cross exported to interoperable genotype, phenotype, and map CSV files."; };
            templates.Click += delegate { ModernFormatConverter.WriteTemplates(directory); result.Text = "CSV templates created in " + directory; };
            actions.Controls.Add(import); actions.Controls.Add(export); actions.Controls.Add(templates);
            panel.Controls.Add(actions, 1, 3); panel.SetColumnSpan(actions, 2);
            result = new TextBox { Dock = DockStyle.Fill, Multiline = true, ReadOnly = true, ScrollBars = ScrollBars.Vertical };
            panel.Controls.Add(result, 0, 4); panel.SetColumnSpan(result, 3);
            Controls.Add(panel);
        }

        private void Import()
        {
            try
            {
                ModernImportResult imported = ModernFormatConverter.Import(format.Text, primary.Text, secondary.Text, directory, stem);
                result.Text = "Imported " + imported.Individuals + " individuals, " + imported.Markers + " markers, and " + imported.Traits + " traits.\r\n" +
                    string.Join("\r\n", imported.Warnings.ToArray());
            }
            catch (Exception ex) { result.Text = "Import failed: " + ex.Message; }
        }

        private TextBox AddPath(TableLayoutPanel panel, int row, string label)
        {
            TextBox box = new TextBox { Dock = DockStyle.Fill };
            Button browse = new Button { Text = "Browse...", Dock = DockStyle.Fill };
            browse.Click += delegate { using (OpenFileDialog dialog = new OpenFileDialog()) if (dialog.ShowDialog(this) == DialogResult.OK) box.Text = dialog.FileName; };
            panel.Controls.Add(new Label { Text = label, Dock = DockStyle.Fill }, 0, row); panel.Controls.Add(box, 1, row); panel.Controls.Add(browse, 2, row);
            return box;
        }
    }

    internal sealed class TutorialForm : Form
    {
        private readonly Action loadExample;
        private readonly string[] pages =
        {
            "1. Create a project\n\nUse File > New project. Choose Example data to copy the validated sample map and genotype/phenotype data into your project folder.",
            "2. Validate before analysis\n\nOpen Analysis > Pre-analysis diagnostics. Review missingness, allele frequencies, segregation distortion, and phenotype distributions before mapping.",
            "3. Run the workflow\n\nThe workflow queue preserves completed stages, reports elapsed time, and can restart from a failed stage. Every native command remains visible for reproducibility.",
            "4. Interpret peaks\n\nOpen the Results dashboard. Select a peak to see genotype-group phenotype distributions, sample counts, effects, uncertainty bars, and 1.5-LOD, 2-LOD, or bootstrap intervals.",
            "5. Establish significance\n\nUse a deterministic background permutation job. Experiment-wide control is appropriate for genome-wide claims; chromosome-wide control is less conservative and must be reported clearly.",
            "6. Share reproducibly\n\nExport a methods paragraph, interoperable CSV files, and a reproducibility bundle containing inputs, settings, seeds, versions, hashes, logs, and a rerunnable script."
        };
        private int index;
        private Label content;
        private Label progress;

        public TutorialForm(Action loadExample)
        {
            this.loadExample = loadExample;
            Text = "Interactive example tutorial"; Size = new Size(760, 500); StartPosition = FormStartPosition.CenterParent;
            content = new Label { Dock = DockStyle.Fill, Padding = new Padding(35), Font = new Font("Segoe UI", 13F), AutoSize = false };
            progress = new Label { Dock = DockStyle.Top, Height = 32, TextAlign = ContentAlignment.MiddleCenter };
            FlowLayoutPanel actions = new FlowLayoutPanel { Dock = DockStyle.Bottom, Height = 48, FlowDirection = FlowDirection.RightToLeft };
            Button next = new Button { Text = "Next", AutoSize = true }; Button back = new Button { Text = "Back", AutoSize = true };
            Button example = new Button { Text = "Load example analysis", AutoSize = true };
            next.Click += delegate { index = Math.Min(pages.Length - 1, index + 1); UpdatePage(); };
            back.Click += delegate { index = Math.Max(0, index - 1); UpdatePage(); };
            example.Click += delegate { Close(); loadExample(); };
            actions.Controls.Add(next); actions.Controls.Add(back); actions.Controls.Add(example);
            Controls.Add(content); Controls.Add(progress); Controls.Add(actions); UpdatePage();
        }

        private void UpdatePage() { content.Text = pages[index]; progress.Text = "Tutorial " + (index + 1) + " of " + pages.Length; }
    }

    internal sealed class AnalysisDesignForm : Form
    {
        private readonly string directory;
        private readonly string stem;
        private ComboBox model;
        private TextBox covariates;
        private TextBox environment;
        private CheckBox interaction;
        private CheckBox structure;
        private Label guidance;

        public string RecommendedTool
        {
            get
            {
                if (model.SelectedIndex == 1) return "JZmapqtl";
                if (model.SelectedIndex == 2) return "MultiRegress";
                return "Zmapqtl";
            }
        }

        public AnalysisDesignForm(string directory, string stem)
        {
            this.directory = directory; this.stem = stem;
            Text = "Covariates and analysis design"; Size = new Size(760, 500); StartPosition = FormStartPosition.CenterParent;
            TableLayoutPanel panel = new TableLayoutPanel { Dock = DockStyle.Fill, Padding = new Padding(20), ColumnCount = 2, RowCount = 7 };
            panel.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 190)); panel.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 100));
            model = new ComboBox { Dock = DockStyle.Fill, DropDownStyle = ComboBoxStyle.DropDownList };
            model.Items.AddRange(new object[] { "Composite interval mapping (model 6)", "Joint/multi-trait mapping", "Multiple regression / structure adjustment" }); model.SelectedIndex = 0;
            covariates = new TextBox { Dock = DockStyle.Fill };
            environment = new TextBox { Dock = DockStyle.Fill };
            interaction = new CheckBox { Text = "Request genotype-by-environment terms", AutoSize = true };
            structure = new CheckBox { Text = "Include supplied population-structure covariates", AutoSize = true };
            Add(panel, 0, "Underlying model", model); Add(panel, 1, "Covariate columns", covariates); Add(panel, 2, "Environmental factor", environment);
            panel.Controls.Add(interaction, 1, 3); panel.Controls.Add(structure, 1, 4);
            guidance = new Label { Dock = DockStyle.Fill, AutoSize = true, MaximumSize = new Size(500, 0), ForeColor = Color.FromArgb(75, 85, 95) };
            panel.Controls.Add(guidance, 1, 5);
            Button save = new Button { Text = "Save design and add supported analysis", AutoSize = true };
            save.Click += delegate { Save(); };
            panel.Controls.Add(save, 1, 6);
            model.SelectedIndexChanged += delegate { UpdateGuidance(); };
            interaction.CheckedChanged += delegate { UpdateGuidance(); };
            Controls.Add(panel); UpdateGuidance();
        }

        private void UpdateGuidance()
        {
            bool advanced = model.SelectedIndex > 0;
            interaction.Enabled = advanced;
            structure.Enabled = model.SelectedIndex == 2;
            guidance.Text = model.SelectedIndex == 0
                ? "Model 6 can account for linked background markers. Arbitrary environmental interactions are not added because the native model does not estimate them."
                : model.SelectedIndex == 1
                    ? "JZmapqtl supports joint/multiple-trait interpretation. Encode environmental responses as validated traits or other-trait columns."
                    : "MultiRegress is the appropriate native route for supplied regression covariates. Interaction and population-structure terms require numeric columns in the prepared data.";
        }

        private void Save()
        {
            string path = Path.Combine(directory, (string.IsNullOrEmpty(stem) ? "qtlcart" : stem) + ".analysis-design.json");
            string json = "{\"version\":\"" + ProductInfo.Version + "\",\"model\":\"" + Escape(model.Text) +
                "\",\"covariates\":\"" + Escape(covariates.Text) + "\",\"environment\":\"" + Escape(environment.Text) +
                "\",\"genotypeEnvironment\":" + interaction.Checked.ToString().ToLowerInvariant() +
                ",\"populationStructure\":" + structure.Checked.ToString().ToLowerInvariant() + "}";
            string temp = path + ".tmp"; File.WriteAllText(temp, json, new UTF8Encoding(true));
            if (File.Exists(path)) File.Replace(temp, path, path + ".bak", true); else File.Move(temp, path);
            DialogResult = DialogResult.OK; Close();
        }

        private static void Add(TableLayoutPanel panel, int row, string label, Control control)
        {
            panel.Controls.Add(new Label { Text = label, Dock = DockStyle.Fill, TextAlign = ContentAlignment.MiddleLeft }, 0, row);
            panel.Controls.Add(control, 1, row);
        }

        private static string Escape(string value) { return (value ?? "").Replace("\\", "\\\\").Replace("\"", "\\\""); }
    }
}
