using System;
using System.Collections.Generic;
using System.Drawing;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;

namespace QTLCartographer.Gui
{
    internal sealed class ResultsDashboardForm : Form
    {
        private readonly string directory;
        private readonly string stem;
        private AnalysisResults results;
        private Chart chart;
        private DataGridView table;
        private NumericUpDown thresholdBox;
        private NumericUpDown alphaBox;
        private NumericUpDown permutationsBox;
        private ComboBox metricBox;
        private Label summary;

        public ResultsDashboardForm(string directory, string stem)
        {
            this.directory = directory;
            this.stem = stem;
            Text = "Results dashboard — " + (string.IsNullOrEmpty(stem) ? "QTL Cartographer" : stem);
            Size = new Size(1180, 800);
            MinimumSize = new Size(900, 620);
            StartPosition = FormStartPosition.CenterParent;
            AutoScaleMode = AutoScaleMode.Dpi;
            Build();
            Reload();
        }

        private void Build()
        {
            ToolStrip tools = new ToolStrip { GripStyle = ToolStripGripStyle.Hidden, Dock = DockStyle.Top };
            tools.Items.Add(Button("Refresh", delegate { Reload(); }));
            tools.Items.Add(Button("Export chart PNG", delegate { ExportChart("png"); }));
            tools.Items.Add(Button("Export chart SVG", delegate { ExportChart("svg"); }));
            tools.Items.Add(Button("Export peaks CSV", delegate { ExportCsv(); }));
            tools.Items.Add(Button("HTML report", delegate { ExportReport(); }));

            FlowLayoutPanel settings = new FlowLayoutPanel { Dock = DockStyle.Top, Height = 42, Padding = new Padding(6) };
            settings.Controls.Add(new Label { Text = "Metric", AutoSize = true, Margin = new Padding(3, 8, 3, 3) });
            metricBox = new ComboBox { DropDownStyle = ComboBoxStyle.DropDownList, Width = 70 };
            metricBox.Items.AddRange(new object[] { "LR", "LOD" });
            metricBox.SelectedIndex = 0;
            metricBox.SelectedIndexChanged += delegate { Draw(); };
            settings.Controls.Add(metricBox);
            settings.Controls.Add(new Label { Text = "Threshold", AutoSize = true, Margin = new Padding(14, 8, 3, 3) });
            thresholdBox = new NumericUpDown { DecimalPlaces = 2, Maximum = 100000, Width = 85 };
            thresholdBox.ValueChanged += delegate { Draw(); };
            settings.Controls.Add(thresholdBox);
            settings.Controls.Add(new Label { Text = "α", AutoSize = true, Margin = new Padding(14, 8, 3, 3) });
            alphaBox = new NumericUpDown { DecimalPlaces = 3, Minimum = .001M, Maximum = .500M, Increment = .005M, Value = .050M, Width = 75 };
            alphaBox.ValueChanged += delegate { ApplyEmpiricalThreshold(); };
            settings.Controls.Add(alphaBox);
            settings.Controls.Add(new Label { Text = "Permutations", AutoSize = true, Margin = new Padding(14, 8, 3, 3) });
            permutationsBox = new NumericUpDown { Minimum = 10, Maximum = 100000, Increment = 100, Value = 1000, Width = 80 };
            settings.Controls.Add(permutationsBox);
            Button permutation = new Button { Text = "Run permutation test", AutoSize = true };
            permutation.Click += delegate { RunPermutationTest(); };
            settings.Controls.Add(permutation);
            summary = new Label { AutoSize = true, Margin = new Padding(14, 8, 3, 3), ForeColor = Color.FromArgb(70, 90, 105) };
            settings.Controls.Add(summary);

            SplitContainer split = new SplitContainer { Dock = DockStyle.Fill, Orientation = Orientation.Horizontal, SplitterDistance = 430 };
            chart = new Chart { Dock = DockStyle.Fill, BackColor = Color.White };
            chart.ChartAreas.Add(new ChartArea("Results"));
            chart.Legends.Add(new Legend("Chromosomes"));
            table = new DataGridView
            {
                Dock = DockStyle.Fill, ReadOnly = true, AllowUserToAddRows = false,
                AutoGenerateColumns = true, AutoSizeColumnsMode = DataGridViewAutoSizeColumnsMode.Fill,
                SelectionMode = DataGridViewSelectionMode.FullRowSelect
            };
            split.Panel1.Controls.Add(chart);
            split.Panel2.Controls.Add(table);
            Controls.Add(split);
            Controls.Add(settings);
            Controls.Add(tools);
        }

        private void Reload()
        {
            results = ResultsParser.LoadProject(directory, stem);
            ApplyEmpiricalThreshold();
            BindTable();
            Draw();
        }

        private void ApplyEmpiricalThreshold()
        {
            if (results == null) return;
            double value = results.EmpiricalThreshold((double)alphaBox.Value);
            if (!double.IsNaN(value))
                thresholdBox.Value = (decimal)Math.Min((double)thresholdBox.Maximum, value);
            summary.Text = results.Points.Count + " positions • " + results.Peaks.Count + " peaks • " +
                results.PermutationMaxima.Count + " permutations";
            BindTable();
        }

        private void BindTable()
        {
            if (results == null) return;
            double threshold = (double)thresholdBox.Value;
            foreach (QtlPeak peak in results.Peaks)
                peak.Significant = peak.LikelihoodRatio >= threshold && threshold > 0;
            table.DataSource = results.Peaks.Select(p => new
            {
                Chromosome = p.Chromosome,
                Position_cM = Math.Round(p.PositionCm, 3),
                Flanking_markers = p.LeftMarker + " — " + p.RightMarker,
                LR = Math.Round(p.LikelihoodRatio, 4),
                LOD = Math.Round(p.Lod, 4),
                Additive = Math.Round(p.Additive, 4),
                Dominance = Math.Round(p.Dominance, 4),
                Significant = p.Significant ? "Yes" : "No"
            }).ToList();
        }

        private void Draw()
        {
            if (results == null) return;
            chart.Series.Clear();
            bool lod = metricBox.SelectedItem != null && metricBox.SelectedItem.ToString() == "LOD";
            foreach (IGrouping<int, ResultPoint> chromosome in results.Points.GroupBy(p => p.Chromosome))
            {
                Series series = new Series("Chr " + chromosome.Key)
                {
                    ChartType = SeriesChartType.Line,
                    BorderWidth = 2,
                    XValueType = ChartValueType.Double,
                    YValueType = ChartValueType.Double
                };
                foreach (ResultPoint point in chromosome.OrderBy(p => p.PositionCm))
                    series.Points.AddXY(point.PositionCm, lod ? point.Lod : point.LikelihoodRatio);
                chart.Series.Add(series);
            }
            ChartArea area = chart.ChartAreas[0];
            area.AxisX.Title = "Position (cM)";
            area.AxisY.Title = lod ? "LOD score" : "Likelihood-ratio statistic";
            area.AxisX.MajorGrid.LineColor = Color.Gainsboro;
            area.AxisY.MajorGrid.LineColor = Color.Gainsboro;
            area.AxisY.StripLines.Clear();
            double threshold = (double)thresholdBox.Value;
            if (threshold > 0)
            {
                if (lod) threshold /= (2.0 * Math.Log(10.0));
                area.AxisY.StripLines.Add(new StripLine
                {
                    IntervalOffset = threshold, BorderColor = Color.Firebrick, BorderWidth = 2,
                    BorderDashStyle = ChartDashStyle.Dash,
                    Text = "Threshold " + threshold.ToString("0.00", CultureInfo.InvariantCulture),
                    TextAlignment = StringAlignment.Near
                });
            }
            BindTable();
        }

        private void RunPermutationTest()
        {
            string executable = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "tools", "Zmapqtl.exe");
            if (!File.Exists(executable))
            {
                MessageBox.Show(this, "Zmapqtl.exe is missing.", Text, MessageBoxButtons.OK, MessageBoxIcon.Error);
                return;
            }
            int count = (int)permutationsBox.Value;
            DialogResult confirmation = MessageBox.Show(this,
                "Run " + count + " phenotype permutations? This can take a long time. The dashboard will calculate the " +
                (100 * (1 - (double)alphaBox.Value)).ToString("0.0") + "th percentile of global maxima.",
                "Permutation test", MessageBoxButtons.OKCancel, MessageBoxIcon.Information);
            if (confirmation != DialogResult.OK) return;
            ProcessDialog dialog = new ProcessDialog(executable,
                "-X " + Quote(stem) + " -r " + count + " -A -V", directory, "Permutation test");
            if (dialog.ShowDialog(this) == DialogResult.OK)
                Reload();
        }

        private void ExportChart(string format)
        {
            using (SaveFileDialog dialog = new SaveFileDialog())
            {
                dialog.Filter = format == "png" ? "PNG image (*.png)|*.png" : "SVG image (*.svg)|*.svg";
                dialog.FileName = (string.IsNullOrEmpty(stem) ? "qtlcart" : stem) + "-results." + format;
                if (dialog.ShowDialog(this) != DialogResult.OK) return;
                if (format == "png") chart.SaveImage(dialog.FileName, ChartImageFormat.Png);
                else ResultsReport.WriteSvg(dialog.FileName, results, (double)thresholdBox.Value, metricBox.Text == "LOD");
            }
        }

        private void ExportCsv()
        {
            using (SaveFileDialog dialog = new SaveFileDialog { Filter = "CSV table (*.csv)|*.csv", FileName = stem + "-peaks.csv" })
                if (dialog.ShowDialog(this) == DialogResult.OK)
                    ResultsReport.WriteCsv(dialog.FileName, results.Peaks);
        }

        private void ExportReport()
        {
            using (SaveFileDialog dialog = new SaveFileDialog { Filter = "Self-contained HTML (*.html)|*.html", FileName = stem + "-report.html" })
                if (dialog.ShowDialog(this) == DialogResult.OK)
                {
                    ResultsReport.WriteHtml(dialog.FileName, results, (double)thresholdBox.Value);
                    MessageBox.Show(this, "Report created:\r\n" + dialog.FileName, Text, MessageBoxButtons.OK, MessageBoxIcon.Information);
                }
        }

        private static ToolStripButton Button(string text, EventHandler click)
        {
            ToolStripButton button = new ToolStripButton(text);
            button.Click += click;
            return button;
        }

        private static string Quote(string value) { return "\"" + (value ?? "").Replace("\"", "\\\"") + "\""; }
    }
}
