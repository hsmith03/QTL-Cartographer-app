using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace QTLCartographer.Gui
{
    internal sealed class ModernImportResult
    {
        public string GenotypeCsv { get; set; }
        public string PhenotypeCsv { get; set; }
        public string MapCsv { get; set; }
        public int Individuals { get; set; }
        public int Markers { get; set; }
        public int Traits { get; set; }
        public List<string> Warnings { get; private set; }

        public ModernImportResult() { Warnings = new List<string>(); }
    }

    internal static class ModernFormatConverter
    {
        private static readonly CultureInfo Invariant = CultureInfo.InvariantCulture;

        public static ModernImportResult Import(string format, string primaryFile, string secondaryFile, string outputDirectory, string stem)
        {
            Directory.CreateDirectory(outputDirectory);
            switch ((format ?? "").ToUpperInvariant())
            {
                case "PLINK": return ImportPlink(primaryFile, secondaryFile, outputDirectory, stem);
                case "VCF": return ImportVcf(primaryFile, outputDirectory, stem);
                case "R/QTL2": return ImportRqtl2(primaryFile, secondaryFile, outputDirectory, stem);
                case "R/QTL": return ImportWideCsv(primaryFile, outputDirectory, stem, "R/qtl");
                default: return ImportWideCsv(primaryFile, outputDirectory, stem, "CSV");
            }
        }

        public static void ExportCrossCsv(CrossData data, string directory, string stem)
        {
            string genotype = Path.Combine(directory, stem + "-genotypes.csv");
            string phenotype = Path.Combine(directory, stem + "-phenotypes.csv");
            string map = Path.Combine(directory, stem + "-map.csv");
            StringBuilder g = new StringBuilder("id");
            foreach (MarkerMetadata marker in data.Markers) g.Append(",").Append(Csv(marker.Name));
            g.AppendLine();
            foreach (CrossIndividual individual in data.Individuals)
            {
                g.Append(individual.Id);
                foreach (int value in individual.Genotypes) g.Append(",").Append(value < 0 ? "NA" : value.ToString(Invariant));
                g.AppendLine();
            }
            File.WriteAllText(genotype, g.ToString(), new UTF8Encoding(true));

            StringBuilder p = new StringBuilder("id");
            foreach (string trait in data.TraitNames) p.Append(",").Append(Csv(trait));
            p.AppendLine();
            foreach (CrossIndividual individual in data.Individuals)
            {
                p.Append(individual.Id);
                foreach (double? value in individual.Traits) p.Append(",").Append(value.HasValue ? value.Value.ToString("R", Invariant) : "NA");
                p.AppendLine();
            }
            File.WriteAllText(phenotype, p.ToString(), new UTF8Encoding(true));

            StringBuilder m = new StringBuilder("marker,chr,pos_cM\r\n");
            foreach (MarkerMetadata marker in data.Markers)
                m.Append(Csv(marker.Name)).Append(",").Append(marker.Chromosome).Append(",").AppendLine("");
            File.WriteAllText(map, m.ToString(), new UTF8Encoding(true));
        }

        public static void WriteTemplates(string directory)
        {
            Directory.CreateDirectory(directory);
            File.WriteAllText(Path.Combine(directory, "qtl-genotypes-template.csv"), "id,marker_1,marker_2\r\nsample_1,0,1\r\nsample_2,2,NA\r\n");
            File.WriteAllText(Path.Combine(directory, "qtl-phenotypes-template.csv"), "id,trait_1,covariate_environment\r\nsample_1,5.2,A\r\nsample_2,6.1,B\r\n");
            File.WriteAllText(Path.Combine(directory, "qtl-map-template.csv"), "marker,chr,pos_cM\r\nmarker_1,1,0\r\nmarker_2,1,10\r\n");
        }

        private static ModernImportResult ImportWideCsv(string file, string directory, string stem, string source)
        {
            string[] lines = File.ReadAllLines(file);
            if (lines.Length < 2) throw new InvalidDataException(source + " CSV must contain a header and at least one individual.");
            string[] header = SplitCsv(lines[0]);
            if (header.Length < 2 || !header[0].Equals("id", StringComparison.OrdinalIgnoreCase))
                throw new InvalidDataException(source + " CSV must be wide format with 'id' as the first column.");
            string output = Path.Combine(directory, stem + "-genotypes.csv");
            File.Copy(file, output, true);
            ModernImportResult result = new ModernImportResult
            {
                GenotypeCsv = output, Individuals = lines.Length - 1, Markers = header.Length - 1
            };
            result.Warnings.Add("A separate phenotype CSV is required before native QTL analysis.");
            return result;
        }

        private static ModernImportResult ImportRqtl2(string genotypeFile, string phenotypeFile, string directory, string stem)
        {
            ModernImportResult result = ImportWideCsv(genotypeFile, directory, stem, "R/qtl2 genotype");
            if (!string.IsNullOrEmpty(phenotypeFile) && File.Exists(phenotypeFile))
            {
                string destination = Path.Combine(directory, stem + "-phenotypes.csv");
                File.Copy(phenotypeFile, destination, true);
                result.PhenotypeCsv = destination;
                string[] header = SplitCsv(File.ReadLines(phenotypeFile).First());
                result.Traits = Math.Max(0, header.Length - 1);
                result.Warnings.Clear();
            }
            return result;
        }

        private static ModernImportResult ImportPlink(string pedFile, string mapFile, string directory, string stem)
        {
            if (!File.Exists(pedFile) || !File.Exists(mapFile)) throw new FileNotFoundException("PLINK import requires both PED and MAP files.");
            string[] markers = File.ReadLines(mapFile).Where(line => !string.IsNullOrWhiteSpace(line))
                .Select(line => line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries))
                .Select(p => p.Length > 1 ? p[1] : "marker").ToArray();
            string genotypeFile = Path.Combine(directory, stem + "-genotypes.csv");
            string phenotypeFile = Path.Combine(directory, stem + "-phenotypes.csv");
            StringBuilder genotypes = new StringBuilder("id," + string.Join(",", markers.Select(Csv)) + "\r\n");
            StringBuilder phenotypes = new StringBuilder("id,phenotype\r\n");
            int individuals = 0;
            string[] referenceAlleles = new string[markers.Length];
            foreach (string line in File.ReadLines(pedFile))
            {
                string[] p = line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
                if (p.Length < 6 + markers.Length * 2) continue;
                string id = p[1];
                genotypes.Append(Csv(id));
                for (int m = 0; m < markers.Length; m++)
                {
                    string a = p[6 + m * 2], b = p[7 + m * 2];
                    if (referenceAlleles[m] == null && a != "0") referenceAlleles[m] = a;
                    string value = a == "0" || b == "0" ? "NA" :
                        a != b ? "1" : a == referenceAlleles[m] ? "0" : "2";
                    genotypes.Append(",").Append(value);
                }
                genotypes.AppendLine();
                phenotypes.Append(Csv(id)).Append(",").AppendLine(p[5] == "-9" ? "NA" : p[5]);
                individuals++;
            }
            File.WriteAllText(genotypeFile, genotypes.ToString(), new UTF8Encoding(true));
            File.WriteAllText(phenotypeFile, phenotypes.ToString(), new UTF8Encoding(true));
            string mapDestination = Path.Combine(directory, stem + "-map.csv");
            StringBuilder map = new StringBuilder("marker,chr,pos_cM\r\n");
            foreach (string line in File.ReadLines(mapFile))
            {
                string[] p = line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
                if (p.Length >= 4) map.Append(Csv(p[1])).Append(",").Append(p[0]).Append(",").AppendLine(p[2]);
            }
            File.WriteAllText(mapDestination, map.ToString(), new UTF8Encoding(true));
            return new ModernImportResult { GenotypeCsv = genotypeFile, PhenotypeCsv = phenotypeFile, MapCsv = mapDestination, Individuals = individuals, Markers = markers.Length, Traits = 1 };
        }

        private static ModernImportResult ImportVcf(string file, string directory, string stem)
        {
            string[] samples = new string[0];
            List<string> markers = new List<string>();
            List<string[]> genotypes = new List<string[]>();
            foreach (string line in File.ReadLines(file))
            {
                if (line.StartsWith("#CHROM"))
                {
                    string[] p = line.Split('\t');
                    samples = p.Skip(9).ToArray();
                }
                else if (!line.StartsWith("#"))
                {
                    string[] p = line.Split('\t');
                    if (p.Length < 10) continue;
                    markers.Add(p[2] == "." ? p[0] + ":" + p[1] : p[2]);
                    string[] row = new string[samples.Length];
                    for (int i = 0; i < samples.Length; i++)
                    {
                        string gt = p[9 + i].Split(':')[0].Replace('|', '/');
                        row[i] = gt == "0/0" ? "0" : gt == "0/1" || gt == "1/0" ? "1" : gt == "1/1" ? "2" : "NA";
                    }
                    genotypes.Add(row);
                }
            }
            if (samples.Length == 0) throw new InvalidDataException("VCF contains no sample columns.");
            StringBuilder csv = new StringBuilder("id," + string.Join(",", markers.Select(Csv)) + "\r\n");
            for (int sample = 0; sample < samples.Length; sample++)
            {
                csv.Append(Csv(samples[sample]));
                foreach (string[] marker in genotypes) csv.Append(",").Append(marker[sample]);
                csv.AppendLine();
            }
            string output = Path.Combine(directory, stem + "-genotypes.csv");
            File.WriteAllText(output, csv.ToString(), new UTF8Encoding(true));
            ModernImportResult result = new ModernImportResult { GenotypeCsv = output, Individuals = samples.Length, Markers = markers.Count };
            result.Warnings.Add("VCF genotypes were imported. Add phenotype and genetic-map CSV files before QTL analysis.");
            return result;
        }

        private static string[] SplitCsv(string line) { return line.Split(',').Select(v => v.Trim().Trim('"')).ToArray(); }
        private static string Csv(string value) { return "\"" + (value ?? "").Replace("\"", "\"\"") + "\""; }
    }
}
