using System;
using System.Diagnostics;
using System.IO;
using System.Text.RegularExpressions;

namespace QTLCartographer.Gui
{
    internal static class HelpParser
    {
        private static readonly Regex OptionLine = new Regex(
            @"^\s*\[\s*(?<inside>-\S+(?:\s+.*?)?)\s*\]\s*(?<description>.*?)\s*$",
            RegexOptions.Compiled);

        public static ToolHelp Load(string executable)
        {
            ToolHelp result = new ToolHelp();
            ProcessStartInfo start = new ProcessStartInfo(executable, "-h");
            start.WorkingDirectory = Path.GetDirectoryName(executable);
            start.UseShellExecute = false;
            start.CreateNoWindow = true;
            start.RedirectStandardOutput = true;
            start.RedirectStandardError = true;

            using (Process process = Process.Start(start))
            {
                string output = process.StandardOutput.ReadToEnd();
                string error = process.StandardError.ReadToEnd();
                if (!process.WaitForExit(10000))
                {
                    process.Kill();
                    throw new InvalidOperationException("Timed out while reading program options.");
                }

                result.FullText = output + error;
            }

            using (StringReader reader = new StringReader(result.FullText))
            {
                string line;
                while ((line = reader.ReadLine()) != null)
                {
                    if (line.StartsWith("PURPOSE:", StringComparison.OrdinalIgnoreCase))
                        result.Purpose = line.Substring(8).Trim();

                    Match match = OptionLine.Match(line);
                    if (!match.Success)
                        continue;

                    string inside = match.Groups["inside"].Value.Trim();
                    int separator = IndexOfWhitespace(inside);
                    string flag = separator < 0 ? inside : inside.Substring(0, separator);
                    string value = separator < 0 ? "" : inside.Substring(separator).Trim();
                    result.Options.Add(new OptionDefinition
                    {
                        Flag = flag,
                        DefaultValue = value,
                        Description = match.Groups["description"].Value.Trim()
                    });
                }
            }

            return result;
        }

        private static int IndexOfWhitespace(string value)
        {
            for (int i = 0; i < value.Length; i++)
                if (char.IsWhiteSpace(value[i]))
                    return i;
            return -1;
        }
    }
}
