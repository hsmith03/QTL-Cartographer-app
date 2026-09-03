using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace QTLCartographer.Gui
{
    internal sealed class CrossIndividual
    {
        public int Id { get; set; }
        public int[] Genotypes { get; set; }
        public double?[] Traits { get; set; }
        public Dictionary<string, double?> Covariates { get; private set; }

        public CrossIndividual() { Covariates = new Dictionary<string, double?>(StringComparer.OrdinalIgnoreCase); }
    }

    internal sealed class MarkerMetadata
    {
        public int Chromosome { get; set; }
        public int Marker { get; set; }
        public int GlobalIndex { get; set; }
        public string Name { get; set; }
        public double PositionCm { get; set; }
    }

    internal sealed class CrossData
    {
        public string CrossType { get; set; }
        public List<string> TraitNames { get; private set; }
        public List<string> CovariateNames { get; private set; }
        public List<MarkerMetadata> Markers { get; private set; }
        public List<CrossIndividual> Individuals { get; private set; }

        public CrossData()
        {
            CrossType = "";
            TraitNames = new List<string>();
            CovariateNames = new List<string>();
            Markers = new List<MarkerMetadata>();
            Individuals = new List<CrossIndividual>();
        }

        public MarkerMetadata FindMarker(int chromosome, int marker)
        {
            return Markers.FirstOrDefault(m => m.Chromosome == chromosome && m.Marker == marker);
        }
    }

    internal sealed class GenotypeSummary
    {
        public int Genotype { get; set; }
        public int Count { get; set; }
        public double Mean { get; set; }
        public double StandardError { get; set; }
        public List<double> Values { get; set; }
    }

    internal static class CrossDataParser
    {
        private static readonly CultureInfo Invariant = CultureInfo.InvariantCulture;

        public static CrossData Load(string crossFile, string mapFile)
        {
            CrossData data = new CrossData();
            if (File.Exists(mapFile)) ParseMarkers(mapFile, data);
            if (!File.Exists(crossFile)) return data;
            string[] lines = File.ReadAllLines(crossFile);
            int markerCount = data.Markers.Count;
            int traitCount = 0;
            bool inData = false;
            for (int i = 0; i < lines.Length; i++)
            {
                string line = lines[i].Trim();
                Match meta;
                if ((meta = Regex.Match(line, @"^-p\s+(\d+)")).Success && markerCount == 0)
                    markerCount = Math.Max(0, int.Parse(meta.Groups[1].Value) - 1);
                else if ((meta = Regex.Match(line, @"^-traits\s+(\d+)")).Success)
                    traitCount = int.Parse(meta.Groups[1].Value);
                else if ((meta = Regex.Match(line, @"^-cross\s+(\S+)")).Success)
                    data.CrossType = meta.Groups[1].Value;
                else if (line.StartsWith("-Names of the traits", StringComparison.OrdinalIgnoreCase))
                {
                    for (int t = 0; t < traitCount && i + 1 < lines.Length; t++)
                    {
                        string[] name = Regex.Split(lines[++i].Trim(), @"\s+");
                        data.TraitNames.Add(name.Length > 1 ? name[1] : "Trait " + (t + 1));
                    }
                }
                else if (line == "-s") inData = true;
                else if (inData)
                {
                    Match start = Regex.Match(line, @"^(\d+)\s+[12]\s*$");
                    if (!start.Success) continue;
                    CrossIndividual individual = new CrossIndividual { Id = int.Parse(start.Groups[1].Value) };
                    List<string> values = new List<string>();
                    int j = i + 1;
                    while (j < lines.Length && !Regex.IsMatch(lines[j].Trim(), @"^\d+\s+[12]\s*$"))
                    {
                        values.AddRange(Regex.Split(lines[j].Trim(), @"\s+").Where(v => v.Length > 0));
                        j++;
                    }
                    individual.Genotypes = new int[markerCount];
                    for (int m = 0; m < markerCount; m++)
                    {
                        int genotype;
                        individual.Genotypes[m] = m < values.Count && int.TryParse(values[m], out genotype) ? genotype : -1;
                    }
                    individual.Traits = new double?[traitCount];
                    for (int t = 0; t < traitCount; t++)
                    {
                        double phenotype;
                        string token = markerCount + t < values.Count ? values[markerCount + t] : ".";
                        if (double.TryParse(token, NumberStyles.Float, Invariant, out phenotype)) individual.Traits[t] = phenotype;
                    }
                    data.Individuals.Add(individual);
                    i = j - 1;
                }
            }
            while (data.TraitNames.Count < traitCount) data.TraitNames.Add("Trait " + (data.TraitNames.Count + 1));
            return data;
        }

        public static List<GenotypeSummary> Summarize(CrossData data, int chromosome, int marker, int traitIndex)
        {
            MarkerMetadata metadata = data.FindMarker(chromosome, marker);
            if (metadata == null || traitIndex < 0) return new List<GenotypeSummary>();
            return data.Individuals
                .Where(i => metadata.GlobalIndex < i.Genotypes.Length && i.Genotypes[metadata.GlobalIndex] >= 0 &&
                            traitIndex < i.Traits.Length && i.Traits[traitIndex].HasValue)
                .GroupBy(i => i.Genotypes[metadata.GlobalIndex])
                .OrderBy(g => g.Key)
                .Select(g =>
                {
                    List<double> values = g.Select(i => i.Traits[traitIndex].Value).ToList();
                    double mean = values.Average();
                    double variance = values.Count > 1 ? values.Sum(v => (v - mean) * (v - mean)) / (values.Count - 1) : 0;
                    return new GenotypeSummary
                    {
                        Genotype = g.Key, Count = values.Count, Mean = mean,
                        StandardError = values.Count > 0 ? Math.Sqrt(variance / values.Count) : 0,
                        Values = values
                    };
                }).ToList();
        }

        private static void ParseMarkers(string mapFile, CrossData data)
        {
            string[] lines = File.ReadAllLines(mapFile);
            int chromosomeCount = 0;
            bool distancesAreMorgans = false;
            foreach (string raw in lines)
            {
                Match chromosomeMatch = Regex.Match(raw.Trim(), @"^-c\s+(\d+)");
                if (chromosomeMatch.Success) chromosomeCount = int.Parse(chromosomeMatch.Groups[1].Value);
                Match unitMatch = Regex.Match(raw.Trim(), @"^-u\s+(\S+)");
                if (unitMatch.Success) distancesAreMorgans = unitMatch.Groups[1].Value.StartsWith("m", StringComparison.OrdinalIgnoreCase);
            }

            Dictionary<int, Dictionary<int, double>> intervals = new Dictionary<int, Dictionary<int, double>>();
            for (int chromosome = 1; chromosome <= chromosomeCount; chromosome++)
                intervals[chromosome] = new Dictionary<int, double>();
            foreach (string raw in lines)
            {
                Match intervalMatch = Regex.Match(raw, @"^-l\s+(\d+)\s+\|");
                if (!intervalMatch.Success) continue;
                int interval = int.Parse(intervalMatch.Groups[1].Value);
                string fields = raw.Substring(intervalMatch.Index + intervalMatch.Length);
                for (int chromosome = 1; chromosome <= chromosomeCount; chromosome++)
                {
                    int start = chromosome == 1 ? 0 : 9 + (chromosome - 2) * 8;
                    int width = chromosome == 1 ? 9 : 8;
                    if (start >= fields.Length) continue;
                    string field = fields.Substring(start, Math.Min(width, fields.Length - start)).Trim();
                    double distance;
                    if (double.TryParse(field, NumberStyles.Float, Invariant, out distance))
                        intervals[chromosome][interval] = distancesAreMorgans ? distance * 100.0 : distance;
                }
            }

            bool active = false;
            foreach (string raw in lines)
            {
                string line = raw.Trim();
                if (line.StartsWith("-b", StringComparison.OrdinalIgnoreCase) && line.Contains("MarkerNames")) { active = true; continue; }
                if (active && line.StartsWith("-e", StringComparison.OrdinalIgnoreCase)) break;
                if (!active) continue;
                string[] p = Regex.Split(line, @"\s+");
                int chromosome;
                int marker;
                if (p.Length >= 3 && int.TryParse(p[0], out chromosome) && int.TryParse(p[1], out marker))
                {
                    double position = 0;
                    Dictionary<int, double> chromosomeIntervals;
                    if (intervals.TryGetValue(chromosome, out chromosomeIntervals))
                        for (int interval = 1; interval < marker; interval++)
                        {
                            double distance;
                            if (chromosomeIntervals.TryGetValue(interval, out distance)) position += distance;
                        }
                    data.Markers.Add(new MarkerMetadata
                    {
                        Chromosome = chromosome, Marker = marker,
                        GlobalIndex = data.Markers.Count, Name = p[2], PositionCm = position
                    });
                }
            }
        }
    }
}
