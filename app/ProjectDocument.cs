using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Json;

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

        public ProjectDocument()
        {
            Version = ProductInfo.Version;
            SelectedTool = "";
            SelectedArguments = "";
            Queue = new List<ProjectStep>();
            Results = new List<string>();
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
            using (FileStream stream = File.Create(path))
                serializer.WriteObject(stream, document);
        }

        public static ProjectDocument Load(string path)
        {
            DataContractJsonSerializer serializer = new DataContractJsonSerializer(typeof(ProjectDocument));
            using (FileStream stream = File.OpenRead(path))
                return (ProjectDocument)serializer.ReadObject(stream);
        }
    }
}
