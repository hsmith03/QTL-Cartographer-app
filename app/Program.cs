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
                try { FeatureSmokeTest(args[1]); }
                catch (Exception ex)
                {
                    System.IO.File.WriteAllText(System.IO.Path.Combine(args[1], "feature-smoke-error.txt"), ex.ToString());
                    Environment.ExitCode = 1;
                }
                return;
            }
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            try
            {
                Application.Run(new MainForm());
            }
            catch (Exception ex)
            {
                string crashDirectory = System.IO.Path.Combine(
                    Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData),
                    "QTL Cartographer");
                System.IO.Directory.CreateDirectory(crashDirectory);
                string crashFile = System.IO.Path.Combine(crashDirectory, "qtlcart-gui-crash.txt");
                System.IO.File.WriteAllText(crashFile, ex.ToString());
                MessageBox.Show(
                    "QTL Cartographer could not start.\r\n\r\nTechnical details were written to:\r\n" + crashFile,
                    "QTL Cartographer startup error",
                    MessageBoxButtons.OK,
                    MessageBoxIcon.Error);
            }
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
            if (results.Peaks.Any(p => p.Support15Right < p.Support15Left || p.Support20Right < p.Support20Left))
                throw new InvalidOperationException("QTL support interval calculation failed.");

            CrossData cross = CrossDataParser.Load(
                System.IO.Path.Combine(directory, "qtlcart.cro"),
                System.IO.Path.Combine(directory, "qtlcart.map"));
            if (cross.Individuals.Count != 333 || cross.Markers.Count != 12 || cross.TraitNames.Count != 1)
                throw new InvalidOperationException("Genotype/phenotype effect parser failed.");
            if (CrossDataParser.Summarize(cross, 1, 1, 0).Count == 0)
                throw new InvalidOperationException("Marker genotype effect summary failed.");

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
            loaded.InputHashes = ProjectStore.HashInputs(directory);
            ProjectStore.Save(project, loaded);
            if (!System.IO.File.Exists(project + ".bak"))
                throw new InvalidOperationException("Atomic project backup was not created.");
            string hashProbe = System.IO.Path.Combine(directory, "hash-probe.inp");
            System.IO.File.WriteAllText(hashProbe, "original");
            loaded.InputHashes = ProjectStore.HashInputs(directory);
            System.IO.File.AppendAllText(hashProbe, "-changed");
            if (!ProjectStore.DetectChangedInputs(loaded).Contains("hash-probe.inp"))
                throw new InvalidOperationException("Changed input hash was not detected.");

            string methods = ScientificMethods.Generate(directory, "qtlcart", results, cross, "Experiment-wide", 13.0, .05);
            if (!methods.Contains(ProductInfo.Version) || !methods.Contains("333"))
                throw new InvalidOperationException("Methods paragraph generation failed.");

            string interoperability = System.IO.Path.Combine(directory, "interoperability");
            ModernFormatConverter.WriteTemplates(interoperability);
            ModernFormatConverter.ExportCrossCsv(cross, interoperability, "sample");
            if (!System.IO.File.Exists(System.IO.Path.Combine(interoperability, "sample-genotypes.csv")))
                throw new InvalidOperationException("Modern CSV export failed.");
            TestModernImports(interoperability);

            string reference = System.IO.Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "example", "sample-golden-peaks.csv");
            BenchmarkResult benchmark = ScientificBenchmarkSuite.Compare(results, reference, .02, .002);
            if (!benchmark.Passed)
                throw new InvalidOperationException("Golden scientific benchmark failed: " + benchmark.Summary);

            string bundle = System.IO.Path.Combine(directory, "feature-reproducibility.zip");
            ReproducibilityBundle.Create(directory, "qtlcart", project,
                new[] { new CommandRequest { Tool = new ToolDefinition("Zmapqtl", "", ""), Arguments = "-X qtlcart -A", WorkingDirectory = directory } }, bundle);
            if (!System.IO.File.Exists(bundle) || new System.IO.FileInfo(bundle).Length == 0)
                throw new InvalidOperationException("Reproducibility bundle export failed.");

            FuzzResultParser(directory);
            TestCrossTypeValidation(directory);
        }

        private static void TestModernImports(string directory)
        {
            string vcf = System.IO.Path.Combine(directory, "tiny.vcf");
            System.IO.File.WriteAllText(vcf,
                "##fileformat=VCFv4.2\r\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\r\n1\t10\tm1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\r\n");
            ModernImportResult vcfResult = ModernFormatConverter.Import("VCF", vcf, "", directory, "vcf-test");
            if (vcfResult.Individuals != 2 || vcfResult.Markers != 1)
                throw new InvalidOperationException("VCF import failed.");

            string rqtl2 = System.IO.Path.Combine(directory, "rqtl2.csv");
            System.IO.File.WriteAllText(rqtl2, "id,m1,m2\r\ns1,0,1\r\ns2,2,NA\r\n");
            ModernImportResult csvResult = ModernFormatConverter.Import("R/qtl2", rqtl2,
                System.IO.Path.Combine(directory, "qtl-phenotypes-template.csv"), directory, "rqtl2-test");
            if (csvResult.Individuals != 2 || csvResult.Markers != 2)
                throw new InvalidOperationException("R/qtl2 import failed.");

            string map = System.IO.Path.Combine(directory, "tiny.map");
            string ped = System.IO.Path.Combine(directory, "tiny.ped");
            System.IO.File.WriteAllText(map, "1 m1 0 10\r\n1 m2 5 20\r\n");
            System.IO.File.WriteAllText(ped, "f s1 0 0 1 5.2 A A A G\r\nf s2 0 0 2 6.1 G G G G\r\n");
            ModernImportResult plink = ModernFormatConverter.Import("PLINK", ped, map, directory, "plink-test");
            if (plink.Individuals != 2 || plink.Markers != 2 || plink.Traits != 1)
                throw new InvalidOperationException("PLINK import failed.");

            string malformed = System.IO.Path.Combine(directory, "malformed.csv");
            System.IO.File.WriteAllText(malformed, "wrong,header\r\n");
            bool rejected = false;
            try { ModernFormatConverter.Import("CSV", malformed, "", directory, "bad"); }
            catch (System.IO.InvalidDataException) { rejected = true; }
            if (!rejected) throw new InvalidOperationException("Malformed modern CSV was not rejected.");
            ImportInspection missing = DataImportValidator.Inspect(
                System.IO.Path.Combine(directory, "does-not-exist.map"),
                System.IO.Path.Combine(directory, "does-not-exist.cross"), "");
            if (missing.IsValid) throw new InvalidOperationException("Missing import files were not rejected.");
        }

        private static void FuzzResultParser(string directory)
        {
            Random random = new Random(314159);
            for (int iteration = 0; iteration < 50; iteration++)
            {
                string path = System.IO.Path.Combine(directory, "fuzz-" + iteration + ".z");
                System.Text.StringBuilder text = new System.Text.StringBuilder();
                for (int line = 0; line < 100; line++)
                {
                    int length = random.Next(0, 80);
                    for (int i = 0; i < length; i++) text.Append((char)random.Next(32, 127));
                    text.AppendLine();
                }
                System.IO.File.WriteAllText(path, text.ToString());
                AnalysisResults parsed = new AnalysisResults();
                ResultsParser.ParseZ(path, parsed);
            }
        }

        private static void TestCrossTypeValidation(string directory)
        {
            string map = System.IO.Path.Combine(directory, "cross-types-map.inp");
            System.IO.File.WriteAllText(map, "-chromosomes 2\r\n-Chromosome 1\r\n-marker m1\r\n-Chromosome 2\r\n-marker m2\r\n");
            foreach (string crossType in new[] { "B1", "B2", "F2", "RI0", "RI1", "SF2", "RF2" })
            {
                string cross = System.IO.Path.Combine(directory, "cross-" + crossType + ".inp");
                System.IO.File.WriteAllText(cross, "-Cross " + crossType + "\r\n-traits 2\r\n-SampleSize 10000\r\n-start markers\r\n");
                ImportInspection inspection = DataImportValidator.Inspect(map, cross, crossType);
                if (!inspection.IsValid || inspection.Individuals != 10000 || inspection.Traits != 2)
                    throw new InvalidOperationException("Cross-type validation failed for " + crossType + ".");
            }
        }
    }
}
