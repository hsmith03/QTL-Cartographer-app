using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;

namespace QTLCartographer.Gui
{
    internal sealed class DiagnosticsForm : Form
    {
        private readonly CrossData data;

        public DiagnosticsForm(CrossData data)
        {
            this.data = data;
            Text = "Pre-analysis diagnostics";
            Size = new Size(1200, 820);
            MinimumSize = new Size(900, 650);
            StartPosition = FormStartPosition.CenterParent;
            AutoScaleMode = AutoScaleMode.Dpi;
            Build();
        }

        private void Build()
        {
            TabControl tabs = new TabControl { Dock = DockStyle.Fill, Multiline = true };
            tabs.TabPages.Add(Page("Missingness heatmap", BuildMissingness()));
            tabs.TabPages.Add(Page("Allele frequencies", BuildAlleleFrequencies()));
            tabs.TabPages.Add(Page("Segregation distortion", BuildDistortion()));
            tabs.TabPages.Add(Page("Phenotype histograms", BuildPhenotypes()));
            Controls.Add(tabs);
        }

        private Control BuildMissingness()
        {
            DataGridView grid = new DataGridView
            {
                Dock = DockStyle.Fill, ReadOnly = true, AllowUserToAddRows = false,
                RowHeadersWidth = 75, AutoSizeColumnsMode = DataGridViewAutoSizeColumnsMode.None,
                AccessibleName = "Marker missingness heatmap"
            };
            int markers = Math.Min(100, data.Markers.Count);
            foreach (MarkerMetadata marker in data.Markers.Take(markers)) grid.Columns.Add(marker.Name, marker.Name);
            foreach (CrossIndividual individual in data.Individuals.Take(300))
            {
                int row = grid.Rows.Add();
                grid.Rows[row].HeaderCell.Value = individual.Id.ToString();
                for (int m = 0; m < markers; m++)
                {
                    bool missing = m >= individual.Genotypes.Length || individual.Genotypes[m] < 0;
                    grid.Rows[row].Cells[m].Value = missing ? "X" : "";
                    grid.Rows[row].Cells[m].Style.BackColor = missing ? Color.FromArgb(213, 94, 0) : Color.FromArgb(230, 245, 250);
                    grid.Columns[m].Width = 38;
                }
            }
            return grid;
        }

        private Control BuildAlleleFrequencies()
        {
            Chart chart = Chart("Allele frequency", "Marker", "Frequency");
            Series series = new Series("Alternate allele frequency") { ChartType = SeriesChartType.Column, Color = Color.FromArgb(0, 114, 178) };
            for (int marker = 0; marker < data.Markers.Count; marker++)
            {
                List<int> values = data.Individuals.Where(i => marker < i.Genotypes.Length && i.Genotypes[marker] >= 0).Select(i => i.Genotypes[marker]).ToList();
                double frequency = values.Count == 0 ? 0 : values.Sum() / (2.0 * values.Count);
                DataPoint point = new DataPoint(marker + 1, frequency) { AxisLabel = data.Markers[marker].Name };
                series.Points.Add(point);
            }
            chart.Series.Add(series);
            return chart;
        }

        private Control BuildDistortion()
        {
            DataGridView grid = new DataGridView { Dock = DockStyle.Fill, ReadOnly = true, AutoGenerateColumns = true, AutoSizeColumnsMode = DataGridViewAutoSizeColumnsMode.Fill };
            var rows = new List<object>();
            for (int marker = 0; marker < data.Markers.Count; marker++)
            {
                int[] counts = new int[3];
                foreach (CrossIndividual individual in data.Individuals)
                    if (marker < individual.Genotypes.Length && individual.Genotypes[marker] >= 0 && individual.Genotypes[marker] <= 2)
                        counts[individual.Genotypes[marker]]++;
                int n = counts.Sum();
                double[] expected = data.CrossType.IndexOf("F2", StringComparison.OrdinalIgnoreCase) >= 0
                    ? new[] { n * .25, n * .5, n * .25 } : new[] { n * .5, 0.0, n * .5 };
                double chi = 0;
                for (int i = 0; i < 3; i++) if (expected[i] > 0) chi += (counts[i] - expected[i]) * (counts[i] - expected[i]) / expected[i];
                double approximateP = Math.Exp(-chi / 2.0);
                rows.Add(new { Marker = data.Markers[marker].Name, AA = counts[0], AB = counts[1], BB = counts[2], Chi_square = Math.Round(chi, 4), Approx_P = Math.Round(approximateP, 6), Flag = approximateP < .01 ? "Review" : "" });
            }
            grid.DataSource = rows;
            return grid;
        }

        private Control BuildPhenotypes()
        {
            Chart chart = Chart("Phenotypes", "Phenotype", "Count");
            for (int trait = 0; trait < data.TraitNames.Count; trait++)
            {
                List<double> values = data.Individuals.Where(i => trait < i.Traits.Length && i.Traits[trait].HasValue).Select(i => i.Traits[trait].Value).ToList();
                if (values.Count == 0) continue;
                int bins = Math.Max(5, (int)Math.Sqrt(values.Count));
                double min = values.Min(), max = values.Max(), width = Math.Max(.0001, (max - min) / bins);
                int[] counts = new int[bins];
                foreach (double value in values) counts[Math.Min(bins - 1, (int)((value - min) / width))]++;
                Series series = new Series(data.TraitNames[trait]) { ChartType = SeriesChartType.Line, BorderWidth = 2 };
                for (int b = 0; b < bins; b++) series.Points.AddXY(min + (b + .5) * width, counts[b]);
                chart.Series.Add(series);
            }
            return chart;
        }

        private static Chart Chart(string name, string x, string y)
        {
            Chart chart = new Chart { Dock = DockStyle.Fill, AccessibleName = name + " chart" };
            ChartArea area = new ChartArea(name); area.AxisX.Title = x; area.AxisY.Title = y;
            chart.ChartAreas.Add(area); chart.Legends.Add(new Legend());
            return chart;
        }

        private static TabPage Page(string title, Control content) { TabPage page = new TabPage(title); page.Controls.Add(content); return page; }
    }
}
