using System;
using System.Linq;
using System.Windows.Forms;

namespace QTLCartographer.Gui
{
    internal static class Program
    {
        [STAThread]
        private static void Main(string[] args)
        {
            if (Environment.OSVersion.Version.Major >= 6)
                SetProcessDPIAware();
            if (args.Length == 1 && args[0] == "--smoke-test")
            {
                foreach (ToolDefinition tool in ToolCatalog.Create())
                {
                    string executable = System.IO.Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "tools", tool.Name + ".exe");
                    ToolHelp help = HelpParser.Load(executable);
                    if (help.Options.Count == 0)
                        throw new InvalidOperationException("No options were discovered for " + tool.Name + ".");
                }
                ProjectStateSmokeTest();
                return;
            }
            if (args.Length == 2 && args[0] == "--feature-smoke-test")
            {
                FeatureSmokeTest(args[1]);
                return;
            }
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new MainForm());
        }

        [System.Runtime.InteropServices.DllImport("user32.dll")]
        private static extern bool SetProcessDPIAware();

        private static void ProjectStateSmokeTest()
        {
            string directory = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "qtlcart-gui-state-" + Guid.NewGuid().ToString("N"));
            System.IO.Directory.CreateDirectory(directory);
            try
            {
                System.IO.File.WriteAllText(System.IO.Path.Combine(directory, "qtlcart.rc"),
                    "-stem qtlcart # legacy default\r\n-map five_chromosomes.map # linkage map\r\n-chrom 5\r\n");
                ProjectState state = ProjectState.Load(directory, "");
                if (state.Stem != "five_chromosomes")
                    throw new InvalidOperationException("Project filename stem was not synchronized from the resource file.");
                if (!state.ResolveProjectFile(directory, "-map", "").EndsWith("five_chromosomes.map"))
                    throw new InvalidOperationException("Project map was not loaded from the resource file.");
            }
            finally
            {
                System.IO.Directory.Delete(directory, true);
            }
        }

        private static void FeatureSmokeTest(string directory)
        {
            System.IO.File.WriteAllText(System.IO.Path.Combine(directory, "qtlcart.z3e"),
                "# permutation global maxima\r\n1 10.0\r\n2 12.0\r\n3 14.0\r\n4 16.0\r\n");
            AnalysisResults results = ResultsParser.LoadProject(directory, "qtlcart");
            if (results.Points.Count == 0 || results.Peaks.Count == 0)
                throw new InvalidOperationException("Results dashboard parser did not load Zmapqtl/Eqtl output.");
            if (results.Points.Select(p => p.Chromosome).Distinct().Count() < 2)
                throw new InvalidOperationException("Results parser collapsed chromosome series.");
            if (results.EmpiricalThreshold(.05) != 16.0)
                throw new InvalidOperationException("Empirical permutation threshold calculation failed.");

            string csv = System.IO.Path.Combine(directory, "feature-test-peaks.csv");
            string svg = System.IO.Path.Combine(directory, "feature-test-chart.svg");
            string html = System.IO.Path.Combine(directory, "feature-test-report.html");
            ResultsReport.WriteCsv(csv, results.Peaks);
            ResultsReport.WriteSvg(svg, results, 13.0, false);
            ResultsReport.WriteHtml(html, results, 13.0);
            foreach (string file in new[] { csv, svg, html })
                if (!System.IO.File.Exists(file) || new System.IO.FileInfo(file).Length == 0)
                    throw new InvalidOperationException("Results export failed: " + file);

            ImportInspection inspection = DataImportValidator.Inspect(
                System.IO.Path.Combine(directory, "sample.mps"),
                System.IO.Path.Combine(directory, "sample.raw"), "F2");
            if (!inspection.IsValid || inspection.Chromosomes != 2)
                throw new InvalidOperationException("Guided import validation failed.");

            string project = System.IO.Path.Combine(directory, "feature-test.qtlproject");
            ProjectDocument document = new ProjectDocument { WorkingDirectory = directory, Stem = "qtlcart" };
            document.Queue.Add(new ProjectStep { Tool = "Zmapqtl", Arguments = "-A", Status = "Succeeded" });
            ProjectStore.Save(project, document);
            ProjectDocument loaded = ProjectStore.Load(project);
            if (loaded.Stem != "qtlcart" || loaded.Queue.Count != 1)
                throw new InvalidOperationException("Project save/open round trip failed.");
        }
    }
}
