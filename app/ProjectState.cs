using System;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;

namespace QTLCartographer.Gui
{
    internal sealed class ProjectState
    {
        private static readonly Regex ResourceLine = new Regex(
            @"^\s*(?<key>-\S+)\s+(?<value>\S+)", RegexOptions.Compiled);
        private static readonly Regex OutputLine = new Regex(
            @"This output file \((?<file>[^)]+)\)", RegexOptions.Compiled | RegexOptions.IgnoreCase);

        public string Stem { get; private set; }
        public Dictionary<string, string> Values { get; private set; }

        private ProjectState()
        {
            Stem = "";
            Values = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
        }

        public static ProjectState Load(string workingDirectory, string requestedResource)
        {
            ProjectState state = new ProjectState();
            string resource = ResolvePath(workingDirectory,
                string.IsNullOrWhiteSpace(requestedResource) ? "qtlcart.rc" : requestedResource);
            if (File.Exists(resource))
            {
                foreach (string line in File.ReadAllLines(resource))
                {
                    Match match = ResourceLine.Match(line);
                    if (match.Success)
                        state.Values[match.Groups["key"].Value] = match.Groups["value"].Value;
                }
            }

            string value;
            if (state.Values.TryGetValue("-stem", out value))
                state.Stem = value;

            if (string.IsNullOrEmpty(state.Stem) ||
                string.Equals(state.Stem, "qtlcart", StringComparison.OrdinalIgnoreCase))
            {
                foreach (string key in new[] { "-map", "-ifile", "-error" })
                {
                    if (state.Values.TryGetValue(key, out value))
                    {
                        string inferred = StemFromFile(value);
                        if (!string.IsNullOrEmpty(inferred) &&
                            !string.Equals(inferred, "qtlcart", StringComparison.OrdinalIgnoreCase))
                        {
                            state.Stem = inferred;
                            break;
                        }
                    }
                }
            }

            if (string.IsNullOrEmpty(state.Stem) ||
                string.Equals(state.Stem, "qtlcart", StringComparison.OrdinalIgnoreCase))
            {
                string logStem = StemFromLatestLog(workingDirectory);
                if (!string.IsNullOrEmpty(logStem))
                    state.Stem = logStem;
            }
            return state;
        }

        public string ResolveProjectFile(string workingDirectory, string key, string fallback)
        {
            string value;
            if (!Values.TryGetValue(key, out value) || string.IsNullOrWhiteSpace(value))
                value = fallback;
            return ResolvePath(workingDirectory, value);
        }

        public static string StemFromFile(string value)
        {
            if (string.IsNullOrWhiteSpace(value))
                return "";
            return Path.GetFileNameWithoutExtension(value.Trim().Trim('"'));
        }

        private static string StemFromLatestLog(string directory)
        {
            if (!Directory.Exists(directory))
                return "";
            FileInfo latest = null;
            foreach (string path in Directory.GetFiles(directory, "*.log"))
            {
                FileInfo candidate = new FileInfo(path);
                if (latest == null || candidate.LastWriteTimeUtc > latest.LastWriteTimeUtc)
                    latest = candidate;
            }
            if (latest == null)
                return "";
            string lastStem = "";
            foreach (string line in File.ReadAllLines(latest.FullName))
            {
                Match match = OutputLine.Match(line);
                if (match.Success)
                {
                    string stem = StemFromFile(match.Groups["file"].Value);
                    if (!string.IsNullOrEmpty(stem))
                        lastStem = stem;
                }
            }
            return lastStem;
        }

        private static string ResolvePath(string directory, string value)
        {
            string path = (value ?? "").Trim().Trim('"');
            if (path.Length == 0)
                return "";
            return Path.IsPathRooted(path) ? path : Path.Combine(directory ?? "", path);
        }
    }
}
