using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Json;
using System.Security.Cryptography;

namespace QTLCartographer.Gui
{
    [DataContract]
    internal sealed class ProjectDocument
    {
        [DataMember] public string Version { get; set; }
        [DataMember] public string WorkingDirectory { get; set; }
        [DataMember] public string Stem { get; set; }
        [DataMember] public string ResourceFile { get; set; }
        [DataMember] public string SelectedTool { get; set; }
        [DataMember] public string SelectedArguments { get; set; }
        [DataMember] public List<ProjectStep> Queue { get; set; }
        [DataMember] public List<string> Results { get; set; }
        [DataMember] public Dictionary<string, string> InputHashes { get; set; }

        public ProjectDocument()
        {
            Version = ProductInfo.Version;
            SelectedTool = "";
            SelectedArguments = "";
            Queue = new List<ProjectStep>();
            Results = new List<string>();
            InputHashes = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
        }
    }

    [DataContract]
    internal sealed class ProjectStep
    {
        [DataMember] public string Tool { get; set; }
        [DataMember] public string Arguments { get; set; }
        [DataMember] public string Status { get; set; }
        [DataMember] public double ElapsedSeconds { get; set; }
    }

    internal static class ProjectStore
    {
        public static void Save(string path, ProjectDocument document)
        {
            DataContractJsonSerializer serializer = new DataContractJsonSerializer(typeof(ProjectDocument));
            string temp = path + ".tmp";
            using (FileStream stream = File.Create(temp))
                serializer.WriteObject(stream, document);
            if (File.Exists(path))
                File.Replace(temp, path, path + ".bak", true);
            else
                File.Move(temp, path);
        }

        public static ProjectDocument Load(string path)
        {
            DataContractJsonSerializer serializer = new DataContractJsonSerializer(typeof(ProjectDocument));
            using (FileStream stream = File.OpenRead(path))
                return (ProjectDocument)serializer.ReadObject(stream);
        }

        public static Dictionary<string, string> HashInputs(string directory)
        {
            Dictionary<string, string> hashes = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            if (!Directory.Exists(directory)) return hashes;
            string[] extensions = { ".inp", ".raw", ".mps", ".map", ".cro", ".csv", ".ped", ".vcf" };
            using (SHA256 sha = SHA256.Create())
            {
                foreach (string file in Directory.GetFiles(directory).Where(path => extensions.Contains(Path.GetExtension(path).ToLowerInvariant())))
                {
                    using (FileStream stream = File.OpenRead(file))
                        hashes[Path.GetFileName(file)] = BitConverter.ToString(sha.ComputeHash(stream)).Replace("-", "").ToLowerInvariant();
                }
            }
            return hashes;
        }

        public static List<string> DetectChangedInputs(ProjectDocument document)
        {
            Dictionary<string, string> current = HashInputs(document.WorkingDirectory);
            List<string> changed = new List<string>();
            foreach (KeyValuePair<string, string> entry in document.InputHashes ?? new Dictionary<string, string>())
            {
                string value;
                if (!current.TryGetValue(entry.Key, out value) || !string.Equals(value, entry.Value, StringComparison.OrdinalIgnoreCase))
                    changed.Add(entry.Key);
            }
            return changed;
        }
    }
}
