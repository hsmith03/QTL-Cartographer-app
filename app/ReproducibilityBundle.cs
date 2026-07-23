using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Json;
using System.Security.Cryptography;
using System.Text;

namespace QTLCartographer.Gui
{
    [DataContract]
    internal sealed class ReproducibilityManifest
    {
        [DataMember] public string ApplicationVersion { get; set; }
        [DataMember] public string NativeEngineVersion { get; set; }
        [DataMember] public string CreatedUtc { get; set; }
        [DataMember] public string OperatingSystem { get; set; }
        [DataMember] public string Stem { get; set; }
        [DataMember] public Dictionary<string, string> Sha256 { get; set; }
        [DataMember] public List<string> Commands { get; set; }
    }

    internal static class ReproducibilityBundle
    {
        public static void Create(string directory, string stem, string projectFile, IEnumerable<CommandRequest> requests, string destination)
        {
            ReproducibilityManifest manifest = new ReproducibilityManifest
            {
                ApplicationVersion = ProductInfo.Version,
                NativeEngineVersion = "QTL Cartographer 1.17j",
                CreatedUtc = DateTime.UtcNow.ToString("o"),
                OperatingSystem = Environment.OSVersion.ToString(),
                Stem = stem,
                Sha256 = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase),
                Commands = requests.Select(request => request.DisplayCommand).ToList()
            };
            string temp = destination + ".tmp";
            if (File.Exists(temp)) File.Delete(temp);
            using (FileStream stream = File.Create(temp))
            using (ZipArchive archive = new ZipArchive(stream, ZipArchiveMode.Create))
            using (SHA256 sha = SHA256.Create())
            {
                foreach (string file in Directory.GetFiles(directory))
                {
                    if (string.Equals(Path.GetFullPath(file), Path.GetFullPath(destination), StringComparison.OrdinalIgnoreCase)) continue;
                    if (string.Equals(Path.GetFullPath(file), Path.GetFullPath(temp), StringComparison.OrdinalIgnoreCase)) continue;
                    string name = Path.GetFileName(file);
                    using (FileStream input = File.Open(file, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))
                        manifest.Sha256[name] = BitConverter.ToString(sha.ComputeHash(input)).Replace("-", "").ToLowerInvariant();
                    archive.CreateEntryFromFile(file, "project/" + name, CompressionLevel.Optimal);
                }
                if (!string.IsNullOrEmpty(projectFile) && File.Exists(projectFile) &&
                    !directory.Equals(Path.GetDirectoryName(projectFile), StringComparison.OrdinalIgnoreCase))
                    archive.CreateEntryFromFile(projectFile, "project/" + Path.GetFileName(projectFile), CompressionLevel.Optimal);
                WriteEntry(archive, "rerun-analysis.cmd", BuildScript(manifest.Commands));
                WriteEntry(archive, "README.txt",
                    "QTL Cartographer reproducibility bundle\r\n\r\nExtract this archive on Windows. Review rerun-analysis.cmd, ensure the packaged QTL Cartographer tools are on PATH, and run it from the project directory.\r\nChecksums are stored in manifest.json.\r\n");
                ZipArchiveEntry manifestEntry = archive.CreateEntry("manifest.json");
                using (Stream output = manifestEntry.Open())
                    new DataContractJsonSerializer(typeof(ReproducibilityManifest)).WriteObject(output, manifest);
            }
            if (File.Exists(destination)) File.Replace(temp, destination, destination + ".bak", true);
            else File.Move(temp, destination);
        }

        private static string BuildScript(IEnumerable<string> commands)
        {
            StringBuilder script = new StringBuilder("@echo off\r\nsetlocal\r\ncd /d \"%~dp0project\"\r\n");
            foreach (string command in commands) script.AppendLine(command);
            script.AppendLine("if errorlevel 1 (echo Analysis failed. & exit /b 1)");
            script.AppendLine("echo Analysis completed.");
            return script.ToString();
        }

        private static void WriteEntry(ZipArchive archive, string name, string content)
        {
            ZipArchiveEntry entry = archive.CreateEntry(name);
            using (StreamWriter writer = new StreamWriter(entry.Open(), new UTF8Encoding(false))) writer.Write(content);
        }
    }
}
