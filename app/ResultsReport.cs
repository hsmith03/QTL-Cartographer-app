using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net;
using System.Text;

namespace QTLCartographer.Gui
{
    internal static class ResultsReport
    {
        private static readonly CultureInfo Invariant = CultureInfo.InvariantCulture;

        public static void WriteCsv(string path, IEnumerable<QtlPeak> peaks)
        {
            StringBuilder text = new StringBuilder("Chromosome,Position_cM,Left_marker,Right_marker,Support_1.5_LOD_left,Support_1.5_LOD_right,Support_2_LOD_left,Support_2_LOD_right,Bootstrap_95_left,Bootstrap_95_right,LR,LOD,Additive,Dominance,Significant\r\n");
            foreach (QtlPeak p in peaks)
                text.AppendFormat(Invariant, "{0},{1:0.####},{2},{3},{4:0.####},{5:0.####},{6:0.####},{7:0.####},{8:0.####},{9:0.####},{10:0.####},{11:0.####},{12:0.####},{13:0.####},{14}\r\n",
                    p.Chromosome, p.PositionCm, Csv(p.LeftMarker), Csv(p.RightMarker),
                    p.Support15Left, p.Support15Right, p.Support20Left, p.Support20Right,
                    p.BootstrapLeft, p.BootstrapRight, p.LikelihoodRatio, p.Lod,
                    p.Additive, p.Dominance, p.Significant ? "Yes" : "No");
            File.WriteAllText(path, text.ToString(), new UTF8Encoding(true));
        }

        public static void WriteSvg(string path, AnalysisResults results, double threshold, bool lod)
        {
            const double left = 75, top = 35, plotWidth = 1080, plotHeight = 585;
            double maxX = Math.Max(1, results.Points.Count == 0 ? 1 : results.Points.Max(p => p.PositionCm));
            double maxY = Math.Max(1, results.Points.Count == 0 ? 1 : results.Points.Max(p => lod ? p.Lod : p.LikelihoodRatio));
            StringBuilder svg = new StringBuilder();
            svg.AppendLine("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1200\" height=\"700\" viewBox=\"0 0 1200 700\">");
            svg.AppendLine("<rect width=\"100%\" height=\"100%\" fill=\"white\"/><g font-family=\"Arial\" font-size=\"16\">");
            svg.AppendFormat(Invariant, "<line x1=\"{0}\" y1=\"{1}\" x2=\"{0}\" y2=\"{2}\" stroke=\"#333\"/><line x1=\"{0}\" y1=\"{2}\" x2=\"{3}\" y2=\"{2}\" stroke=\"#333\"/>",
                left, top, top + plotHeight, left + plotWidth).AppendLine();
            foreach (QtlPeak peak in results.Peaks)
            {
                double bandX = left + peak.Support15Left / maxX * plotWidth;
                double bandWidth = Math.Max(1, (peak.Support15Right - peak.Support15Left) / maxX * plotWidth);
                svg.AppendFormat(Invariant, "<rect x=\"{0:0.0}\" y=\"{1}\" width=\"{2:0.0}\" height=\"{3}\" fill=\"#009e73\" fill-opacity=\"0.12\"/>",
                    bandX, top, bandWidth, plotHeight).AppendLine();
            }
            string[] colors = { "#0072b2", "#d55e00", "#009e73", "#cc79a7", "#e69f00", "#56b4e9", "#000000" };
            int color = 0;
            foreach (var group in results.Points.GroupBy(p => p.Chromosome))
            {
                string points = string.Join(" ", group.OrderBy(p => p.PositionCm).Select(p =>
                    (left + p.PositionCm / maxX * plotWidth).ToString("0.0", Invariant) + "," +
                    (top + plotHeight - (lod ? p.Lod : p.LikelihoodRatio) / maxY * plotHeight).ToString("0.0", Invariant)));
                svg.AppendFormat("<polyline fill=\"none\" stroke=\"{0}\" stroke-width=\"2\" points=\"{1}\"/>", colors[color++ % colors.Length], points).AppendLine();
            }
            double shownThreshold = lod ? threshold / (2 * Math.Log(10)) : threshold;
            if (shownThreshold > 0)
            {
                double y = top + plotHeight - shownThreshold / maxY * plotHeight;
                svg.AppendFormat(Invariant, "<line x1=\"{0}\" y1=\"{1:0.0}\" x2=\"{2}\" y2=\"{1:0.0}\" stroke=\"#d55e00\" stroke-width=\"2\" stroke-dasharray=\"8 5\"/>", left, y, left + plotWidth).AppendLine();
            }
            svg.AppendFormat("<text x=\"600\" y=\"680\" text-anchor=\"middle\">Position (cM)</text><text x=\"18\" y=\"350\" transform=\"rotate(-90 18 350)\" text-anchor=\"middle\">{0}</text>", lod ? "LOD score" : "Likelihood ratio");
            svg.AppendLine("</g></svg>");
            File.WriteAllText(path, svg.ToString(), new UTF8Encoding(false));
        }

        public static void WriteHtml(string path, AnalysisResults results, double threshold)
        {
            string svgPath = Path.GetTempFileName();
            WriteSvg(svgPath, results, threshold, false);
            string svg = File.ReadAllText(svgPath);
            File.Delete(svgPath);
            StringBuilder rows = new StringBuilder();
            foreach (QtlPeak p in results.Peaks)
                rows.AppendFormat(Invariant,
                    "<tr><td>{0}</td><td>{1:0.###}</td><td>{2}</td><td>{3:0.###}-{4:0.###}</td><td>{5:0.###}-{6:0.###}</td><td>{7:0.###}</td><td>{8:0.###}</td><td>{9:0.###}</td><td>{10}</td></tr>",
                    p.Chromosome, p.PositionCm, WebUtility.HtmlEncode(p.LeftMarker + " - " + p.RightMarker),
                    p.Support15Left, p.Support15Right, p.BootstrapLeft, p.BootstrapRight,
                    p.LikelihoodRatio, p.Lod, p.Additive, p.Significant ? "Yes" : "No");
            string html = "<!doctype html><html><head><meta charset=\"utf-8\"><title>QTL Cartographer report</title>" +
                "<style>body{font:15px Arial;max-width:1200px;margin:30px auto;color:#1a2b44}table{border-collapse:collapse;width:100%}th,td{border:1px solid #ccd6df;padding:7px;text-align:right}th{background:#edf3f7}svg{max-width:100%;height:auto}</style></head><body>" +
                "<h1>QTL Cartographer analysis report</h1><p>Windows GUI " + ProductInfo.Version + " | Trait: " +
                WebUtility.HtmlEncode(results.Trait) + " | Empirical/manual LR threshold: " + threshold.ToString("0.###", Invariant) +
                " | Permutations: " + results.PermutationMaxima.Count + "</p>" + svg +
                "<h2>Detected peaks</h2><table><thead><tr><th>Chromosome</th><th>Position (cM)</th><th>Flanking markers</th><th>1.5-LOD interval</th><th>Bootstrap 95% CI</th><th>LR</th><th>LOD</th><th>Additive</th><th>Significant</th></tr></thead><tbody>" +
                rows + "</tbody></table><p>Generated " + DateTime.Now.ToString("u") + ". Interpret peaks relative to a validated significance threshold.</p></body></html>";
            File.WriteAllText(path, html, new UTF8Encoding(true));
        }

        private static string Csv(string value) { return "\"" + (value ?? "").Replace("\"", "\"\"") + "\""; }
    }
}
