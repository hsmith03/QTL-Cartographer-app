using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace QTLCartographer.Gui
{
    internal static class ScientificMethods
    {
        public static string Generate(string directory, string stem, AnalysisResults results, CrossData data,
            string thresholdMethod, double threshold, double alpha)
        {
            string model = ReadResource(directory, "-Model", "unspecified");
            string window = ReadResource(directory, "-window", "unspecified");
            string background = ReadResource(directory, "-background", "unspecified");
            int sampleSize = data == null ? 0 : data.Individuals.Count;
            string trait = results == null || string.IsNullOrWhiteSpace(results.Trait) ? "the selected quantitative trait" : results.Trait;
            return string.Format(CultureInfo.InvariantCulture,
                "Quantitative trait loci for {0} were analyzed with QTL Cartographer for Windows {1} " +
                "(validated native engine 1.17j) using model {2}, a {3}-cM window, and {4} background parameters. " +
                "The analysis included {5} individuals from a {6} cross. Statistical significance was assessed using " +
                "{7} control at alpha={8:0.###}, corresponding to an LR threshold of {9:0.###}. " +
                "Peak uncertainty was summarized with 1.5-LOD and 2-LOD support intervals and, where available, bootstrap percentile intervals.",
                trait, ProductInfo.Version, model, window, background, sampleSize,
                data == null || string.IsNullOrEmpty(data.CrossType) ? "unspecified" : data.CrossType,
                thresholdMethod.ToLowerInvariant(), alpha, threshold);
        }

        private static string ReadResource(string directory, string key, string fallback)
        {
            string file = Path.Combine(directory, "qtlcart.rc");
            if (!File.Exists(file)) return fallback;
            string line = File.ReadLines(file).LastOrDefault(value => value.TrimStart().StartsWith(key, StringComparison.OrdinalIgnoreCase));
            if (line == null) return fallback;
            Match match = Regex.Match(line, @"^\s*\S+\s+(\S+)");
            return match.Success ? match.Groups[1].Value : fallback;
        }
    }
}
