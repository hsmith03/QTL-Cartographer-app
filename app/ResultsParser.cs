using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace QTLCartographer.Gui
{
    internal sealed class ResultPoint
    {
        public int Chromosome { get; set; }
        public int Marker { get; set; }
        public double PositionCm { get; set; }
        public double LikelihoodRatio { get; set; }
        public double Lod { get { return LikelihoodRatio / (2.0 * Math.Log(10.0)); } }
        public double Additive { get; set; }
        public double Dominance { get; set; }
    }

    internal sealed class QtlPeak
    {
        public int Chromosome { get; set; }
        public int Marker { get; set; }
        public double PositionCm { get; set; }
        public double LikelihoodRatio { get; set; }
        public double Lod { get { return LikelihoodRatio / (2.0 * Math.Log(10.0)); } }
        public double Additive { get; set; }
        public double Dominance { get; set; }
        public string LeftMarker { get; set; }
        public string RightMarker { get; set; }
        public bool Significant { get; set; }
        public double Support15Left { get; set; }
        public double Support15Right { get; set; }
        public double Support20Left { get; set; }
        public double Support20Right { get; set; }
        public double BootstrapLeft { get; set; }
        public double BootstrapRight { get; set; }
    }

    internal sealed class AnalysisResults
    {
        public string SourceFile { get; set; }
        public string Trait { get; set; }
        public List<ResultPoint> Points { get; private set; }
        public List<QtlPeak> Peaks { get; private set; }
        public List<double> PermutationMaxima { get; private set; }
        public Dictionary<int, List<double>> ChromosomePermutationMaxima { get; private set; }

        public AnalysisResults()
        {
            Trait = "";
            Points = new List<ResultPoint>();
            Peaks = new List<QtlPeak>();
            PermutationMaxima = new List<double>();
            ChromosomePermutationMaxima = new Dictionary<int, List<double>>();
        }

        public double EmpiricalThreshold(double alpha)
        {
            if (PermutationMaxima.Count == 0)
                return double.NaN;
            List<double> sorted = PermutationMaxima.OrderBy(x => x).ToList();
            int index = (int)Math.Ceiling((1.0 - alpha) * sorted.Count) - 1;
            return sorted[Math.Max(0, Math.Min(sorted.Count - 1, index))];
        }

        public double ChromosomeThreshold(int chromosome, double alpha)
        {
            List<double> values;
            if (chromosome <= 0 || !ChromosomePermutationMaxima.TryGetValue(chromosome, out values) || values.Count == 0)
                return EmpiricalThreshold(alpha);
            List<double> sorted = values.OrderBy(x => x).ToList();
            int index = (int)Math.Ceiling((1.0 - alpha) * sorted.Count) - 1;
            return sorted[Math.Max(0, Math.Min(sorted.Count - 1, index))];
        }
    }

    internal static class ResultsParser
    {
        private static readonly CultureInfo Invariant = CultureInfo.InvariantCulture;

        public static AnalysisResults LoadProject(string directory, string stem)
        {
            AnalysisResults results = new AnalysisResults();
            string z = Find(directory, stem, ".z");
            if (!string.IsNullOrEmpty(z))
                ParseZ(z, results);
            string eqt = Find(directory, stem, ".eqt");
            if (!string.IsNullOrEmpty(eqt))
                ParseEqt(eqt, results);
            string map = Find(directory, stem, ".map");
            string activeStem = string.IsNullOrEmpty(stem) ? "qtlcart" : stem;
            foreach (string file in Directory.Exists(directory)
                ? Directory.GetFiles(directory).Where(path =>
                    Path.GetFileName(path).StartsWith(activeStem + ".z", StringComparison.OrdinalIgnoreCase) &&
                    path.EndsWith("e", StringComparison.OrdinalIgnoreCase)) : new string[0])
                ParsePermutationMaxima(file, results);
            if (results.Peaks.Count == 0)
                DetectLocalPeaks(results);
            if (!string.IsNullOrEmpty(map))
                ApplyMarkerNames(map, results);
            CalculateSupportIntervals(results);
            ApplyBootstrapIntervals(directory, stem, results);
            return results;
        }

        public static void ParseZ(string file, AnalysisResults results)
        {
            results.SourceFile = file;
            foreach (string raw in File.ReadLines(file))
            {
                string line = raw.Trim();
                if (line.StartsWith("-trait", StringComparison.OrdinalIgnoreCase))
                {
                    Match match = Regex.Match(line, @"\[([^\]]+)\]");
                    if (match.Success) results.Trait = match.Groups[1].Value;
                }
                if (line.Length == 0 || line[0] == '#' || line[0] == '-')
                    continue;
                string[] p = Regex.Split(line, @"\s+");
                if (p.Length < 11)
                    continue;
                int chrom;
                int marker;
                double position;
                double lr;
                double additive;
                double dominance;
                if (!int.TryParse(p[0], out chrom) || !int.TryParse(p[1], out marker) ||
                    !double.TryParse(p[2], NumberStyles.Float, Invariant, out position) ||
                    !double.TryParse(p[10], NumberStyles.Float, Invariant, out lr) ||
                    !double.TryParse(p[6], NumberStyles.Float, Invariant, out additive) ||
                    !double.TryParse(p[8], NumberStyles.Float, Invariant, out dominance))
                    continue;
                results.Points.Add(new ResultPoint
                {
                    Chromosome = chrom,
                    Marker = marker,
                    PositionCm = position * 100.0,
                    LikelihoodRatio = lr,
                    Additive = additive,
                    Dominance = dominance
                });
            }
        }

        public static void ParseEqt(string file, AnalysisResults results)
        {
            bool inTable = false;
            foreach (string raw in File.ReadLines(file))
            {
                string line = raw.Trim();
                if (line.Contains("..Chrom..Markr.")) { inTable = true; continue; }
                if (!inTable || line.Length == 0 || line[0] == '#')
                    continue;
                if (line[0] == '-') break;
                string[] p = Regex.Split(line, @"\s+");
                int chrom;
                int marker;
                double position;
                double lr;
                double additive;
                double dominance;
                if (p.Length < 7 || !int.TryParse(p[1], out chrom) || !int.TryParse(p[2], out marker) ||
                    !double.TryParse(p[3], NumberStyles.Float, Invariant, out position) ||
                    !double.TryParse(p[4], NumberStyles.Float, Invariant, out lr) ||
                    !double.TryParse(p[5], NumberStyles.Float, Invariant, out additive) ||
                    !double.TryParse(p[6], NumberStyles.Float, Invariant, out dominance))
                    continue;
                results.Peaks.Add(new QtlPeak
                {
                    Chromosome = chrom,
                    Marker = marker,
                    PositionCm = position,
                    LikelihoodRatio = lr,
                    Additive = additive,
                    Dominance = dominance,
                    LeftMarker = "Marker " + marker,
                    RightMarker = "Marker " + (marker + 1)
                });
            }
        }

        public static void ParsePermutationMaxima(string file, AnalysisResults results)
        {
            foreach (string raw in File.ReadLines(file))
            {
                string line = raw.Trim();
                if (line.Length == 0 || line[0] == '#' || line[0] == '-')
                    continue;
                string[] p = Regex.Split(line, @"\s+");
                double value;
                if (p.Length >= 2 && int.TryParse(p[0], out _) &&
                    double.TryParse(p[1], NumberStyles.Float, Invariant, out value))
                {
                    results.PermutationMaxima.Add(value);
                    for (int chromosome = 1; chromosome + 1 < p.Length; chromosome++)
                    {
                        double chromosomeMaximum;
                        if (!double.TryParse(p[chromosome + 1], NumberStyles.Float, Invariant, out chromosomeMaximum)) continue;
                        if (!results.ChromosomePermutationMaxima.ContainsKey(chromosome))
                            results.ChromosomePermutationMaxima[chromosome] = new List<double>();
                        results.ChromosomePermutationMaxima[chromosome].Add(chromosomeMaximum);
                    }
                }
            }
        }

        private static void ApplyMarkerNames(string file, AnalysisResults results)
        {
            Dictionary<string, string> names = new Dictionary<string, string>();
            bool active = false;
            foreach (string raw in File.ReadLines(file))
            {
                string line = raw.Trim();
                if (line.StartsWith("-b", StringComparison.OrdinalIgnoreCase) && line.Contains("MarkerNames")) { active = true; continue; }
                if (active && line.StartsWith("-e", StringComparison.OrdinalIgnoreCase)) break;
                if (!active) continue;
                string[] p = Regex.Split(line, @"\s+");
                int chrom;
                int marker;
                if (p.Length >= 3 && int.TryParse(p[0], out chrom) && int.TryParse(p[1], out marker))
                    names[chrom + ":" + marker] = p[2];
            }
            foreach (QtlPeak peak in results.Peaks)
            {
                string value;
                if (names.TryGetValue(peak.Chromosome + ":" + peak.Marker, out value)) peak.LeftMarker = value;
                if (names.TryGetValue(peak.Chromosome + ":" + (peak.Marker + 1), out value)) peak.RightMarker = value;
                else peak.RightMarker = "(chromosome end)";
            }
        }

        private static void DetectLocalPeaks(AnalysisResults results)
        {
            foreach (IGrouping<int, ResultPoint> chromosome in results.Points.GroupBy(p => p.Chromosome))
            {
                List<ResultPoint> points = chromosome.OrderBy(p => p.PositionCm).ToList();
                for (int i = 1; i < points.Count - 1; i++)
                {
                    if (points[i].LikelihoodRatio >= points[i - 1].LikelihoodRatio &&
                        points[i].LikelihoodRatio > points[i + 1].LikelihoodRatio)
                    {
                        ResultPoint p = points[i];
                        results.Peaks.Add(new QtlPeak
                        {
                            Chromosome = p.Chromosome,
                            Marker = p.Marker,
                            PositionCm = p.PositionCm,
                            LikelihoodRatio = p.LikelihoodRatio,
                            Additive = p.Additive,
                            Dominance = p.Dominance,
                            LeftMarker = "Marker " + p.Marker,
                            RightMarker = "Marker " + (p.Marker + 1)
                        });
                    }
                }
            }
        }

        private static void CalculateSupportIntervals(AnalysisResults results)
        {
            double lrPerLod = 2.0 * Math.Log(10.0);
            foreach (QtlPeak peak in results.Peaks)
            {
                List<ResultPoint> points = results.Points.Where(p => p.Chromosome == peak.Chromosome)
                    .OrderBy(p => p.PositionCm).ToList();
                peak.Support15Left = FindBoundary(points, peak.PositionCm, peak.LikelihoodRatio - 1.5 * lrPerLod, true);
                peak.Support15Right = FindBoundary(points, peak.PositionCm, peak.LikelihoodRatio - 1.5 * lrPerLod, false);
                peak.Support20Left = FindBoundary(points, peak.PositionCm, peak.LikelihoodRatio - 2.0 * lrPerLod, true);
                peak.Support20Right = FindBoundary(points, peak.PositionCm, peak.LikelihoodRatio - 2.0 * lrPerLod, false);
                peak.BootstrapLeft = peak.Support15Left;
                peak.BootstrapRight = peak.Support15Right;
            }
        }

        private static double FindBoundary(List<ResultPoint> points, double peak, double cutoff, bool left)
        {
            IEnumerable<ResultPoint> side = left
                ? points.Where(p => p.PositionCm <= peak).OrderByDescending(p => p.PositionCm)
                : points.Where(p => p.PositionCm >= peak).OrderBy(p => p.PositionCm);
            ResultPoint boundary = side.FirstOrDefault(p => p.LikelihoodRatio <= cutoff);
            if (boundary != null) return boundary.PositionCm;
            ResultPoint end = side.LastOrDefault();
            return end == null ? peak : end.PositionCm;
        }

        private static void ApplyBootstrapIntervals(string directory, string stem, AnalysisResults results)
        {
            if (!Directory.Exists(directory)) return;
            string prefix = string.IsNullOrEmpty(stem) ? "qtlcart" : stem;
            foreach (string file in Directory.GetFiles(directory).Where(path =>
                Path.GetFileName(path).StartsWith(prefix + ".z", StringComparison.OrdinalIgnoreCase) &&
                path.EndsWith("b", StringComparison.OrdinalIgnoreCase)))
            {
                Dictionary<int, List<double>> positions = new Dictionary<int, List<double>>();
                foreach (string raw in File.ReadLines(file))
                {
                    string[] p = Regex.Split(raw.Trim(), @"\s+");
                    int chromosome;
                    double position;
                    if (p.Length >= 3 && int.TryParse(p[0], out chromosome) &&
                        double.TryParse(p[2], NumberStyles.Float, Invariant, out position))
                    {
                        if (position < 10) position *= 100.0;
                        if (!positions.ContainsKey(chromosome)) positions[chromosome] = new List<double>();
                        positions[chromosome].Add(position);
                    }
                }
                foreach (QtlPeak peak in results.Peaks)
                {
                    List<double> values;
                    if (!positions.TryGetValue(peak.Chromosome, out values) || values.Count < 4) continue;
                    values.Sort();
                    peak.BootstrapLeft = Percentile(values, .025);
                    peak.BootstrapRight = Percentile(values, .975);
                }
            }
        }

        private static double Percentile(List<double> sorted, double probability)
        {
            double index = probability * (sorted.Count - 1);
            int lower = (int)Math.Floor(index);
            int upper = Math.Min(sorted.Count - 1, lower + 1);
            return sorted[lower] + (sorted[upper] - sorted[lower]) * (index - lower);
        }

        private static string Find(string directory, string stem, string extension)
        {
            if (!Directory.Exists(directory))
                return "";
            string exact = Path.Combine(directory, (string.IsNullOrEmpty(stem) ? "qtlcart" : stem) + extension);
            return File.Exists(exact) ? exact : "";
        }
    }
}
