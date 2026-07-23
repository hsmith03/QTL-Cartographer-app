using System;
using System.Windows.Forms;

namespace QTLCartographer.Gui
{
    internal static class Program
    {
        [STAThread]
        private static void Main(string[] args)
        {
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
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new MainForm());
        }

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
    }
}
