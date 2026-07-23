using System;
using System.Diagnostics;
using System.Drawing;
using System.Text;
using System.Windows.Forms;

namespace QTLCartographer.Gui
{
    internal sealed class ProcessDialog : Form
    {
        private readonly string executable;
        private readonly string arguments;
        private readonly string directory;
        private RichTextBox output;
        private ProgressBar progress;
        private Process process;

        public ProcessDialog(string executable, string arguments, string directory, string title)
        {
            this.executable = executable;
            this.arguments = arguments;
            this.directory = directory;
            Text = title;
            Size = new Size(850, 560);
            StartPosition = FormStartPosition.CenterParent;
            AutoScaleMode = AutoScaleMode.Dpi;
            output = new RichTextBox { Dock = DockStyle.Fill, ReadOnly = true, Font = new Font("Consolas", 9F), BackColor = Color.FromArgb(20, 29, 38), ForeColor = Color.White };
            progress = new ProgressBar { Dock = DockStyle.Top, Height = 8, Style = ProgressBarStyle.Marquee };
            Button cancel = new Button { Text = "Cancel", Dock = DockStyle.Bottom, Height = 34 };
            cancel.Click += delegate { if (process != null && !process.HasExited) process.Kill(); };
            Controls.Add(output);
            Controls.Add(progress);
            Controls.Add(cancel);
            Shown += delegate { Start(); };
        }

        private void Start()
        {
            StringBuilder errors = new StringBuilder();
            process = new Process();
            process.StartInfo = new ProcessStartInfo(executable, arguments)
            {
                WorkingDirectory = directory, UseShellExecute = false, CreateNoWindow = true,
                RedirectStandardOutput = true, RedirectStandardError = true
            };
            process.EnableRaisingEvents = true;
            process.OutputDataReceived += delegate(object sender, DataReceivedEventArgs e) { if (e.Data != null) Append(e.Data); };
            process.ErrorDataReceived += delegate(object sender, DataReceivedEventArgs e) { if (e.Data != null) { errors.AppendLine(e.Data); Append(e.Data); } };
            process.Exited += delegate
            {
                process.WaitForExit();
                int code = process.ExitCode;
                BeginInvoke((MethodInvoker)delegate
                {
                    progress.Style = ProgressBarStyle.Blocks;
                    if (code == 0) { DialogResult = DialogResult.OK; Close(); }
                    else
                    {
                        Append("Command failed (exit " + code + "):\r\n" + executable + " " + arguments +
                            "\r\n\r\nSuggested remedy: verify the project map/cross files and review the messages above.");
                    }
                });
            };
            process.Start();
            process.BeginOutputReadLine();
            process.BeginErrorReadLine();
        }

        private void Append(string line)
        {
            if (output.InvokeRequired) { output.BeginInvoke((MethodInvoker)delegate { Append(line); }); return; }
            output.AppendText(line + Environment.NewLine);
            output.ScrollToCaret();
        }
    }
}
