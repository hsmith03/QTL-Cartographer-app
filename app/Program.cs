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
                return;
            }
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new MainForm());
        }
    }
}
