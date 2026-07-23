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
        private static readonly Color[] AccessibleColors =
        {
            Color.FromArgb(0, 114, 178), Color.FromArgb(213, 94, 0),
            Color.FromArgb(0, 158, 115), Color.FromArgb(204, 121, 167),
            Color.FromArgb(230, 159, 0), Color.FromArgb(86, 180, 233),
            Color.FromArgb(240, 228, 66), Color.Black
        };

        private readonly string directory;
        private readonly string stem;
        private AnalysisResults results;
        private CrossData crossData;
        private readonly List<AnalysisResults> comparisons = new List<AnalysisResults>();
        private Chart profileChart;
        private Chart effectChart;
        private Chart phenotypeChart;
        private DataGridView table;
        private NumericUpDown thresholdBox;
        private NumericUpDown alphaBox;
        private NumericUpDown permutationsBox;
        private NumericUpDown seedBox;
        private ComboBox metricBox;
        private ComboBox chromosomeBox;
        private ComboBox thresholdScopeBox;
        private CheckBox support15Box;
        private CheckBox support20Box;
        private CheckBox bootstrapBox;
        private CheckBox markersBox;
        private Label summary;
        private Label effectSummary;
        private Label thresholdExplanation;
        private ProgressBar jobProgress;
        private Label jobStatus;
        private Button pauseJobButton;
        private BackgroundPermutationJob permutationJob;

        public ResultsDashboardForm(string directory, string stem)
        {
            this.directory = directory;
            this.stem = stem;
            Text = "Results dashboard - " + (string.IsNullOrEmpty(stem) ? "QTL Cartographer" : stem);
            Size = new Size(1320, 900);
            MinimumSize = new Size(1080, 720);
            StartPosition = FormStartPosition.CenterParent;
            AutoScaleMode = AutoScaleMode.Dpi;
            Build();
            Reload();
            RestorePermutationJob();
        }

        private void Build()
        {
            ToolStrip tools = new ToolStrip { GripStyle = ToolStripGripStyle.Hidden, Dock = DockStyle.Top };
            tools.Items.Add(Button("Refresh", delegate { Reload(); }));
            tools.Items.Add(Button("Compare run...", delegate { AddComparison(); }));
            tools.Items.Add(new ToolStripSeparator());
            tools.Items.Add(Button("Export PNG", delegate { ExportChart("png"); }));
            tools.Items.Add(Button("Export SVG", delegate { ExportChart("svg"); }));
            tools.Items.Add(Button("Export peaks CSV", delegate { ExportCsv(); }));
            tools.Items.Add(Button("HTML report", delegate { ExportReport(); }));
            tools.Items.Add(Button("Methods paragraph", delegate { ExportMethods(); }));

            FlowLayoutPanel settings = new FlowLayoutPanel
            {
                Dock = DockStyle.Top, Height = 78, Padding = new Padding(6),
                AutoScroll = true, WrapContents = true, AccessibleName = "Results display settings"
            };
            settings.Controls.Add(Label("Metric"));
            metricBox = Choice(70, "LR", "LOD");
            metricBox.SelectedIndexChanged += delegate { Draw(); };
            settings.Controls.Add(metricBox);
            settings.Controls.Add(Label("Chromosome"));
            chromosomeBox = Choice(90, "All");
            chromosomeBox.SelectedIndexChanged += delegate { Draw(); };
            settings.Controls.Add(chromosomeBox);
            settings.Controls.Add(Label("Threshold"));
            thresholdScopeBox = Choice(150, "Experiment-wide", "Chromosome-wide", "Manual");
            thresholdScopeBox.SelectedIndexChanged += delegate { ApplyEmpiricalThreshold(); };
            settings.Controls.Add(thresholdScopeBox);
            thresholdBox = new NumericUpDown { DecimalPlaces = 2, Maximum = 100000, Width = 78, AccessibleName = "Likelihood-ratio threshold" };
            thresholdBox.ValueChanged += delegate { Draw(); };
            settings.Controls.Add(thresholdBox);
            settings.Controls.Add(Label("alpha"));
            alphaBox = new NumericUpDown { DecimalPlaces = 3, Minimum = .001M, Maximum = .500M, Increment = .005M, Value = .050M, Width = 70, AccessibleName = "Significance alpha" };
            alphaBox.ValueChanged += delegate { ApplyEmpiricalThreshold(); };
            settings.Controls.Add(alphaBox);
            support15Box = Check("1.5-LOD interval", true);
            support20Box = Check("2-LOD interval", false);
            bootstrapBox = Check("Bootstrap CI", true);
            markersBox = Check("Marker labels", false);
            foreach (CheckBox box in new[] { support15Box, support20Box, bootstrapBox, markersBox })
            {
                box.CheckedChanged += delegate { Draw(); };
                settings.Controls.Add(box);
            }
            summary = new Label { AutoSize = true, Margin = new Padding(12, 8, 3, 3), ForeColor = Color.FromArgb(70, 90, 105) };
            settings.Controls.Add(summary);

            thresholdExplanation = new Label
            {
                Dock = DockStyle.Top, Height = 30, Padding = new Padding(10, 6, 3, 3),
                ForeColor = Color.FromArgb(75, 85, 95), AutoEllipsis = true,
                AccessibleName = "Multiple testing threshold explanation"
            };

            TabControl views = new TabControl { Dock = DockStyle.Fill, Multiline = true };
            TabPage profilePage = new TabPage("Chromosome profile");
            profileChart = NewChart("Profile");
            profileChart.ChartAreas[0].CursorX.IsUserEnabled = true;
            profileChart.ChartAreas[0].CursorX.IsUserSelectionEnabled = true;
            profileChart.ChartAreas[0].AxisX.ScaleView.Zoomable = true;
            profileChart.MouseWheel += delegate(object sender, MouseEventArgs e) { ZoomProfile(e.Delta); };
            profilePage.Controls.Add(profileChart);

            TabPage effectPage = new TabPage("Effect and uncertainty");
            SplitContainer effectSplit = new SplitContainer { Dock = DockStyle.Fill, Orientation = Orientation.Horizontal, SplitterDistance = 370 };
            effectChart = NewChart("Effects");
            phenotypeChart = NewChart("Phenotype");
            effectSummary = new Label { Dock = DockStyle.Top, Height = 54, Padding = new Padding(8), AutoEllipsis = true, AccessibleName = "Selected peak effect summary" };
            effectSplit.Panel1.Controls.Add(effectChart);
            effectSplit.Panel1.Controls.Add(effectSummary);
            effectSplit.Panel2.Controls.Add(phenotypeChart);
            effectPage.Controls.Add(effectSplit);

            TabPage tablePage = new TabPage("Detected peaks");
            table = new DataGridView
            {
                Dock = DockStyle.Fill, ReadOnly = true, AllowUserToAddRows = false,
                AutoGenerateColumns = true, AutoSizeColumnsMode = DataGridViewAutoSizeColumnsMode.Fill,
                SelectionMode = DataGridViewSelectionMode.FullRowSelect, MultiSelect = false,
                AccessibleName = "Detected QTL peak table"
            };
            table.SelectionChanged += delegate { ShowSelectedPeak(); };
            tablePage.Controls.Add(table);

            TabPage permutationPage = new TabPage("Significance jobs");
            permutationPage.Controls.Add(BuildPermutationPanel());

            views.TabPages.Add(profilePage);
            views.TabPages.Add(effectPage);
            views.TabPages.Add(tablePage);
            views.TabPages.Add(permutationPage);
            Controls.Add(views);
            Controls.Add(thresholdExplanation);
            Controls.Add(settings);
            Controls.Add(tools);
        }

        private Control BuildPermutationPanel()
        {
            TableLayoutPanel panel = new TableLayoutPanel { Dock = DockStyle.Top, Padding = new Padding(20), AutoSize = true, ColumnCount = 2 };
            panel.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 180));
            panel.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 100));
            permutationsBox = new NumericUpDown { Minimum = 10, Maximum = 100000, Increment = 100, Value = 1000, Width = 130 };
            seedBox = new NumericUpDown { Minimum = 1, Maximum = 2147483647, Value = 1357911, Width = 130 };
            panel.Controls.Add(Label("Permutation count"), 0, 0); panel.Controls.Add(permutationsBox, 1, 0);
            panel.Controls.Add(Label("Deterministic seed"), 0, 1); panel.Controls.Add(seedBox, 1, 1);
            jobProgress = new ProgressBar { Dock = DockStyle.Fill, Height = 22 };
            panel.Controls.Add(Label("Progress"), 0, 2); panel.Controls.Add(jobProgress, 1, 2);
            jobStatus = new Label { Dock = DockStyle.Fill, Height = 45, AutoEllipsis = true, Text = "No background permutation job is active." };
            panel.Controls.Add(Label("Status"), 0, 3); panel.Controls.Add(jobStatus, 1, 3);
            FlowLayoutPanel actions = new FlowLayoutPanel { AutoSize = true };
            Button start = new Button { Text = "Start / recover", AutoSize = true };
            pauseJobButton = new Button { Text = "Pause", AutoSize = true, Enabled = false };
            Button cancel = new Button { Text = "Cancel", AutoSize = true, Enabled = false };
            start.Click += delegate { StartPermutationJob(); cancel.Enabled = true; };
            pauseJobButton.Click += delegate { TogglePause(); };
            cancel.Click += delegate { if (permutationJob != null) permutationJob.Cancel(); };
            actions.Controls.Add(start); actions.Controls.Add(pauseJobButton); actions.Controls.Add(cancel);
            panel.Controls.Add(actions, 1, 4);
            Label note = new Label
            {
                AutoSize = true, MaximumSize = new Size(760, 0),
                Text = "Permutation jobs run outside the main UI. Existing partial native output is retained, and Start / recover continues from it. " +
                       "Experiment-wide thresholds use the global maximum per permutation; chromosome-wide thresholds use the selected chromosome maximum."
            };
            panel.Controls.Add(note, 1, 5);
            return panel;
        }

        private void Reload()
        {
            results = ResultsParser.LoadProject(directory, stem);
            string map = FindProjectFile(".map");
            string cross = FindProjectFile(".cro");
            crossData = CrossDataParser.Load(cross, map);
            int selected = chromosomeBox.SelectedIndex;
            chromosomeBox.Items.Clear();
            chromosomeBox.Items.Add("All");
            foreach (int chromosome in results.Points.Select(p => p.Chromosome).Distinct().OrderBy(x => x))
                chromosomeBox.Items.Add("Chr " + chromosome);
            chromosomeBox.SelectedIndex = Math.Min(Math.Max(0, selected), chromosomeBox.Items.Count - 1);
            ApplyEmpiricalThreshold();
            BindTable();
            Draw();
        }

        private void ApplyEmpiricalThreshold()
        {
            if (results == null) return;
            string scope = thresholdScopeBox.Text;
            double value = double.NaN;
            if (scope == "Experiment-wide")
                value = results.EmpiricalThreshold((double)alphaBox.Value);
            else if (scope == "Chromosome-wide")
            {
                int chromosome = SelectedChromosome();
                value = results.ChromosomeThreshold(chromosome, (double)alphaBox.Value);
            }
            thresholdBox.Enabled = scope == "Manual" || double.IsNaN(value);
            if (!double.IsNaN(value)) thresholdBox.Value = (decimal)Math.Min((double)thresholdBox.Maximum, value);
            thresholdExplanation.Text = scope == "Experiment-wide"
                ? "Experiment-wide control compares each permutation's maximum over the whole genome and controls the family-wise false-positive rate."
                : scope == "Chromosome-wide"
                    ? "Chromosome-wide control uses maxima for the selected chromosome; it is less conservative and does not control the genome-wide error rate."
                    : "Manual thresholds are user supplied and should be justified in the analysis report.";
            summary.Text = results.Points.Count + " positions | " + results.Peaks.Count + " peaks | " +
                results.PermutationMaxima.Count + " permutations | " + comparisons.Count + " comparison runs";
            BindTable();
            Draw();
        }

        private void BindTable()
        {
            if (results == null || table == null) return;
            double threshold = (double)thresholdBox.Value;
            foreach (QtlPeak peak in results.Peaks) peak.Significant = peak.LikelihoodRatio >= threshold && threshold > 0;
            table.DataSource = results.Peaks.Select(p => new
            {
                Chromosome = p.Chromosome,
                Position_cM = Math.Round(p.PositionCm, 3),
                Flanking_markers = p.LeftMarker + " - " + p.RightMarker,
                Support_1_5_LOD = Interval(p.Support15Left, p.Support15Right),
                Support_2_LOD = Interval(p.Support20Left, p.Support20Right),
                Bootstrap_95_CI = Interval(p.BootstrapLeft, p.BootstrapRight),
                LR = Math.Round(p.LikelihoodRatio, 4),
                LOD = Math.Round(p.Lod, 4),
                Additive = Math.Round(p.Additive, 4),
                Dominance = Math.Round(p.Dominance, 4),
                Significant = p.Significant ? "Yes" : "No"
            }).ToList();
        }

        private void Draw()
        {
            if (results == null || profileChart == null) return;
            profileChart.Series.Clear();
            ChartArea area = profileChart.ChartAreas[0];
            area.AxisX.StripLines.Clear();
            area.AxisY.StripLines.Clear();
            area.AxisX.CustomLabels.Clear();
            bool lod = metricBox.Text == "LOD";
            int selectedChromosome = SelectedChromosome();
            IEnumerable<IGrouping<int, ResultPoint>> groups = results.Points
                .Where(p => selectedChromosome == 0 || p.Chromosome == selectedChromosome)
                .GroupBy(p => p.Chromosome);
            int colorIndex = 0;
            foreach (IGrouping<int, ResultPoint> chromosome in groups)
                AddProfileSeries(profileChart, chromosome, "Chr " + chromosome.Key, AccessibleColors[colorIndex++ % AccessibleColors.Length], lod, ChartDashStyle.Solid);
            foreach (AnalysisResults comparison in comparisons)
            {
                foreach (IGrouping<int, ResultPoint> chromosome in comparison.Points
                    .Where(p => selectedChromosome == 0 || p.Chromosome == selectedChromosome).GroupBy(p => p.Chromosome))
                    AddProfileSeries(profileChart, chromosome, Path.GetFileNameWithoutExtension(comparison.SourceFile) + " Chr " + chromosome.Key,
                        AccessibleColors[colorIndex++ % AccessibleColors.Length], lod, ChartDashStyle.Dash);
            }
            area.AxisX.Title = "Position (cM)";
            area.AxisY.Title = lod ? "LOD score" : "Likelihood-ratio statistic";
            area.AxisX.MajorGrid.LineColor = Color.Gainsboro;
            area.AxisY.MajorGrid.LineColor = Color.Gainsboro;
            double threshold = (double)thresholdBox.Value;
            if (lod) threshold /= (2.0 * Math.Log(10.0));
            if (threshold > 0)
                area.AxisY.StripLines.Add(new StripLine
                {
                    IntervalOffset = threshold, BorderColor = Color.FromArgb(213, 94, 0), BorderWidth = 2,
                    BorderDashStyle = ChartDashStyle.Dash, Text = "Threshold " + threshold.ToString("0.00", CultureInfo.InvariantCulture)
                });
            DrawConfidenceBands(area, selectedChromosome);
            BindTable();
        }

        private static void AddProfileSeries(Chart target, IEnumerable<ResultPoint> points, string name, Color color, bool lod, ChartDashStyle dash)
        {
            Series series = new Series(name)
            {
                ChartType = SeriesChartType.Line, BorderWidth = 2, Color = color,
                BorderDashStyle = dash, XValueType = ChartValueType.Double, YValueType = ChartValueType.Double
            };
            foreach (ResultPoint point in points.OrderBy(p => p.PositionCm))
            {
                DataPoint data = new DataPoint(point.PositionCm, lod ? point.Lod : point.LikelihoodRatio);
                data.ToolTip = name + "\n" + point.PositionCm.ToString("0.00") + " cM\n" + (lod ? point.Lod : point.LikelihoodRatio).ToString("0.000");
                series.Points.Add(data);
            }
            target.Series.Add(series);
        }

        private void DrawConfidenceBands(ChartArea area, int chromosome)
        {
            IEnumerable<QtlPeak> peaks = results.Peaks.Where(p => chromosome == 0 || p.Chromosome == chromosome);
            foreach (QtlPeak peak in peaks)
            {
                if (support20Box.Checked) AddBand(area, peak.Support20Left, peak.Support20Right, Color.FromArgb(22, 86, 180, 233), "2-LOD");
                if (support15Box.Checked) AddBand(area, peak.Support15Left, peak.Support15Right, Color.FromArgb(34, 0, 158, 115), "1.5-LOD");
                if (bootstrapBox.Checked) AddBand(area, peak.BootstrapLeft, peak.BootstrapRight, Color.FromArgb(28, 204, 121, 167), "Bootstrap 95% CI");
                if (markersBox.Checked)
                {
                    area.AxisX.StripLines.Add(new StripLine
                    {
                        IntervalOffset = peak.PositionCm, BorderWidth = 1, BorderColor = Color.Gray,
                        Text = peak.LeftMarker, TextOrientation = TextOrientation.Rotated270
                    });
                }
            }
        }

        private static void AddBand(ChartArea area, double left, double right, Color color, string label)
        {
            if (right <= left) return;
            area.AxisX.StripLines.Add(new StripLine
            {
                IntervalOffset = left, StripWidth = right - left, BackColor = color,
                Text = label, TextAlignment = StringAlignment.Near
            });
        }

        private void ShowSelectedPeak()
        {
            if (table == null || table.CurrentRow == null || results == null) return;
            int row = table.CurrentRow.Index;
            if (row < 0 || row >= results.Peaks.Count) return;
            QtlPeak peak = results.Peaks[row];
            List<GenotypeSummary> groups = CrossDataParser.Summarize(crossData, peak.Chromosome, peak.Marker, 0);
            DrawEffects(peak, groups);
        }

        private void DrawEffects(QtlPeak peak, List<GenotypeSummary> groups)
        {
            effectChart.Series.Clear();
            phenotypeChart.Series.Clear();
            effectChart.ChartAreas[0].AxisX.Title = "Marker genotype";
            effectChart.ChartAreas[0].AxisY.Title = crossData.TraitNames.Count > 0 ? crossData.TraitNames[0] : "Phenotype";
            Series means = new Series("Mean phenotype") { ChartType = SeriesChartType.Column, Color = AccessibleColors[0] };
            Series errors = new Series("Standard error") { ChartType = SeriesChartType.ErrorBar, Color = Color.Black };
            Series observations = new Series("Individuals") { ChartType = SeriesChartType.Point, MarkerSize = 5, Color = Color.FromArgb(110, AccessibleColors[1]) };
            int index = 1;
            foreach (GenotypeSummary group in groups)
            {
                means.Points.AddXY(index, group.Mean);
                means.Points[means.Points.Count - 1].AxisLabel = GenotypeLabel(group.Genotype) + "\n(n=" + group.Count + ")";
                errors.Points.AddXY(index, group.Mean - group.StandardError, group.Mean + group.StandardError);
                int jitter = 0;
                foreach (double value in group.Values)
                {
                    observations.Points.AddXY(index + ((jitter++ % 7) - 3) * .025, value);
                }
                index++;
            }
            effectChart.Series.Add(means); effectChart.Series.Add(errors); effectChart.Series.Add(observations);
            DrawHistogram(phenotypeChart, crossData.Individuals.Where(i => i.Traits.Length > 0 && i.Traits[0].HasValue).Select(i => i.Traits[0].Value).ToList());
            effectSummary.Text = "Chr " + peak.Chromosome + ", " + peak.PositionCm.ToString("0.00") + " cM | " +
                peak.LeftMarker + " - " + peak.RightMarker + " | 1.5-LOD: " + Interval(peak.Support15Left, peak.Support15Right) +
                " | Bootstrap 95% CI: " + Interval(peak.BootstrapLeft, peak.BootstrapRight) +
                " | additive " + peak.Additive.ToString("0.###") + ", dominance " + peak.Dominance.ToString("0.###");
        }

        private static void DrawHistogram(Chart target, List<double> values)
        {
            target.ChartAreas[0].AxisX.Title = "Phenotype";
            target.ChartAreas[0].AxisY.Title = "Count";
            if (values.Count == 0) return;
            int bins = Math.Max(5, (int)Math.Sqrt(values.Count));
            double min = values.Min(), max = values.Max(), width = Math.Max(.0001, (max - min) / bins);
            int[] counts = new int[bins];
            foreach (double value in values) counts[Math.Min(bins - 1, (int)((value - min) / width))]++;
            Series histogram = new Series("Phenotype distribution") { ChartType = SeriesChartType.Column, Color = AccessibleColors[2] };
            for (int i = 0; i < bins; i++) histogram.Points.AddXY(min + (i + .5) * width, counts[i]);
            target.Series.Add(histogram);
        }

        private void AddComparison()
        {
            using (OpenFileDialog dialog = new OpenFileDialog { Filter = "Zmapqtl results (*.z)|*.z|All files (*.*)|*.*", InitialDirectory = directory, Multiselect = true })
            {
                if (dialog.ShowDialog(this) != DialogResult.OK) return;
                foreach (string file in dialog.FileNames)
                {
                    AnalysisResults comparison = new AnalysisResults();
                    ResultsParser.ParseZ(file, comparison);
                    comparisons.Add(comparison);
                }
                Draw();
            }
        }

        private void StartPermutationJob()
        {
            if (permutationJob != null && permutationJob.IsRunning) return;
            string executable = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "tools", "Zmapqtl.exe");
            permutationJob = new BackgroundPermutationJob(executable, directory, stem, (int)permutationsBox.Value, (int)seedBox.Value);
            permutationJob.ProgressChanged += PermutationProgress;
            permutationJob.Completed += delegate(object sender, EventArgs e) { BeginInvoke((MethodInvoker)delegate { pauseJobButton.Enabled = false; Reload(); }); };
            permutationJob.StartOrRecover();
            pauseJobButton.Enabled = true;
        }

        private void RestorePermutationJob()
        {
            PermutationJobState state = BackgroundPermutationJob.LoadState(directory, stem);
            if (state == null) return;
            permutationsBox.Value = Math.Min(permutationsBox.Maximum, Math.Max(permutationsBox.Minimum, state.Requested));
            seedBox.Value = Math.Min(seedBox.Maximum, Math.Max(seedBox.Minimum, state.Seed));
            jobStatus.Text = state.Completed + " of " + state.Requested + " permutations retained; select Start / recover to continue.";
            jobProgress.Maximum = Math.Max(1, state.Requested);
            jobProgress.Value = Math.Min(jobProgress.Maximum, state.Completed);
        }

        private void PermutationProgress(object sender, PermutationProgressEventArgs e)
        {
            if (InvokeRequired) { BeginInvoke((MethodInvoker)delegate { PermutationProgress(sender, e); }); return; }
            jobProgress.Maximum = Math.Max(1, e.Total);
            jobProgress.Value = Math.Min(jobProgress.Maximum, e.Completed);
            jobStatus.Text = e.Message + " | estimated remaining " + e.EstimatedRemaining;
            pauseJobButton.Text = e.Paused ? "Resume" : "Pause";
        }

        private void TogglePause()
        {
            if (permutationJob == null) return;
            if (permutationJob.IsPaused) permutationJob.Resume(); else permutationJob.Pause();
        }

        private void ZoomProfile(int delta)
        {
            AxisScaleView view = profileChart.ChartAreas[0].AxisX.ScaleView;
            if (delta < 0) { view.ZoomReset(); return; }
            double min = profileChart.ChartAreas[0].AxisX.Minimum;
            double max = profileChart.ChartAreas[0].AxisX.Maximum;
            if (!double.IsNaN(min) && !double.IsNaN(max)) view.Zoom(min + (max - min) * .2, max - (max - min) * .2);
        }

        private int SelectedChromosome()
        {
            if (chromosomeBox == null || chromosomeBox.SelectedIndex <= 0) return 0;
            int value;
            return int.TryParse(chromosomeBox.Text.Replace("Chr ", ""), out value) ? value : 0;
        }

        private void ExportChart(string format)
        {
            using (SaveFileDialog dialog = new SaveFileDialog
            {
                Filter = format == "png" ? "PNG image (*.png)|*.png" : "SVG image (*.svg)|*.svg",
                FileName = stem + "-results." + format
            })
                if (dialog.ShowDialog(this) == DialogResult.OK)
                {
                    if (format == "png") profileChart.SaveImage(dialog.FileName, ChartImageFormat.Png);
                    else ResultsReport.WriteSvg(dialog.FileName, results, (double)thresholdBox.Value, metricBox.Text == "LOD");
                }
        }

        private void ExportCsv()
        {
            using (SaveFileDialog dialog = new SaveFileDialog { Filter = "CSV table (*.csv)|*.csv", FileName = stem + "-peaks.csv" })
                if (dialog.ShowDialog(this) == DialogResult.OK) ResultsReport.WriteCsv(dialog.FileName, results.Peaks);
        }

        private void ExportReport()
        {
            using (SaveFileDialog dialog = new SaveFileDialog { Filter = "Self-contained HTML (*.html)|*.html", FileName = stem + "-report.html" })
                if (dialog.ShowDialog(this) == DialogResult.OK) ResultsReport.WriteHtml(dialog.FileName, results, (double)thresholdBox.Value);
        }

        private void ExportMethods()
        {
            string paragraph = ScientificMethods.Generate(directory, stem, results, crossData, thresholdScopeBox.Text,
                (double)thresholdBox.Value, (double)alphaBox.Value);
            Clipboard.SetText(paragraph);
            MessageBox.Show(this, paragraph + "\r\n\r\nThe paragraph has been copied to the clipboard.", "Methods paragraph", MessageBoxButtons.OK, MessageBoxIcon.Information);
        }

        private string FindProjectFile(string extension)
        {
            string exact = Path.Combine(directory, (string.IsNullOrEmpty(stem) ? "qtlcart" : stem) + extension);
            if (File.Exists(exact)) return exact;
            return Directory.Exists(directory) ? Directory.GetFiles(directory, "*" + extension).OrderByDescending(File.GetLastWriteTimeUtc).FirstOrDefault() ?? "" : "";
        }

        private static Chart NewChart(string name)
        {
            Chart chart = new Chart { Dock = DockStyle.Fill, BackColor = Color.White, AccessibleName = name + " chart", Palette = ChartColorPalette.None };
            chart.ChartAreas.Add(new ChartArea(name));
            chart.Legends.Add(new Legend(name + " legend"));
            return chart;
        }

        private static ToolStripButton Button(string text, EventHandler click) { ToolStripButton button = new ToolStripButton(text); button.Click += click; return button; }
        private static Label Label(string text) { return new Label { Text = text, AutoSize = true, Margin = new Padding(8, 8, 3, 3) }; }
        private static CheckBox Check(string text, bool value) { return new CheckBox { Text = text, Checked = value, AutoSize = true, Margin = new Padding(8, 6, 3, 3) }; }
        private static ComboBox Choice(int width, params string[] items) { ComboBox box = new ComboBox { Width = width, DropDownStyle = ComboBoxStyle.DropDownList }; box.Items.AddRange(items); box.SelectedIndex = 0; return box; }
        private static string Interval(double left, double right) { return left.ToString("0.00") + "-" + right.ToString("0.00") + " cM"; }
        private static string GenotypeLabel(int genotype) { return genotype == 0 ? "AA" : genotype == 1 ? "AB/H" : genotype == 2 ? "BB" : genotype.ToString(); }
    }
}
