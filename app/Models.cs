using System;
using System.Collections.Generic;

namespace QTLCartographer.Gui
{
    internal sealed class ToolDefinition
    {
        public string Name { get; private set; }
        public string Category { get; private set; }
        public string Summary { get; private set; }

        public ToolDefinition(string name, string category, string summary)
        {
            Name = name;
            Category = category;
            Summary = summary;
        }

        public override string ToString()
        {
            return Name;
        }
    }

    internal sealed class OptionDefinition
    {
        public string Flag { get; set; }
        public string DefaultValue { get; set; }
        public string Description { get; set; }
    }

    internal sealed class ToolHelp
    {
        public string Purpose { get; set; }
        public string FullText { get; set; }
        public List<OptionDefinition> Options { get; private set; }

        public ToolHelp()
        {
            Options = new List<OptionDefinition>();
        }
    }

    internal sealed class OptionEditor
    {
        public OptionDefinition Definition { get; set; }
        public System.Windows.Forms.CheckBox Enabled { get; set; }
        public System.Windows.Forms.Control Value { get; set; }
    }

    internal sealed class CommandRequest
    {
        public ToolDefinition Tool { get; set; }
        public string Arguments { get; set; }
        public string WorkingDirectory { get; set; }
        public string ResourceFile { get; set; }
        public string RequestedStem { get; set; }
        public Dictionary<string, string> Options { get; set; }
        public string Status { get; set; }
        public DateTime StartedAt { get; set; }
        public TimeSpan Elapsed { get; set; }
        public string FailureDetails { get; set; }
        public bool OverwriteConfirmed { get; set; }

        public string DisplayCommand
        {
            get { return Tool.Name + ".exe" + (Arguments.Length == 0 ? "" : " " + Arguments); }
        }
    }
}
