using System;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Json;
using System.Threading;

namespace QTLCartographer.Gui
{
    [DataContract]
    internal sealed class PermutationJobState
    {
        [DataMember] public int Requested { get; set; }
        [DataMember] public int Completed { get; set; }
        [DataMember] public int Seed { get; set; }
        [DataMember] public string Status { get; set; }
        [DataMember] public DateTime StartedUtc { get; set; }
    }

    internal sealed class PermutationProgressEventArgs : EventArgs
    {
        public int Completed { get; set; }
        public int Total { get; set; }
        public bool Paused { get; set; }
        public string Message { get; set; }
        public string EstimatedRemaining { get; set; }
    }

    internal sealed class BackgroundPermutationJob
    {
        private readonly string executable;
        private readonly string directory;
        private readonly string stem;
        private readonly int requested;
        private readonly int seed;
        private Process process;
        private Timer timer;
        private DateTime started;
        private int startCount;

        public event EventHandler<PermutationProgressEventArgs> ProgressChanged;
        public event EventHandler Completed;
        public bool IsRunning { get { return process != null && !process.HasExited; } }
        public bool IsPaused { get; private set; }

        public BackgroundPermutationJob(string executable, string directory, string stem, int requested, int seed)
        {
            this.executable = executable;
            this.directory = directory;
            this.stem = string.IsNullOrWhiteSpace(stem) ? "qtlcart" : stem;
            this.requested = requested;
            this.seed = seed;
        }

        public void StartOrRecover()
        {
            if (IsRunning) return;
            startCount = CountCompleted(directory, stem);
            started = DateTime.UtcNow;
            process = new Process();
            process.StartInfo = new ProcessStartInfo(executable,
                "-X \"" + stem.Replace("\"", "\\\"") + "\" -s " + seed.ToString(CultureInfo.InvariantCulture) +
                " -r " + requested.ToString(CultureInfo.InvariantCulture) + " -A -V")
            {
                WorkingDirectory = directory, UseShellExecute = false, CreateNoWindow = true,
                RedirectStandardOutput = true, RedirectStandardError = true
            };
            process.EnableRaisingEvents = true;
            process.OutputDataReceived += delegate { };
            process.ErrorDataReceived += delegate { };
            process.Exited += delegate
            {
                if (timer != null) timer.Dispose();
                Update();
                SaveState(directory, stem, new PermutationJobState
                {
                    Requested = requested, Completed = CountCompleted(directory, stem), Seed = seed,
                    StartedUtc = started, Status = process.ExitCode == 0 ? "Completed" : "Interrupted"
                });
                EventHandler handler = Completed;
                if (handler != null) handler(this, EventArgs.Empty);
            };
            process.Start();
            process.BeginOutputReadLine();
            process.BeginErrorReadLine();
            timer = new Timer(delegate { Update(); }, null, 0, 1000);
        }

        public void Pause()
        {
            if (!IsRunning || IsPaused) return;
            NtSuspendProcess(process.Handle);
            IsPaused = true;
            Update();
        }

        public void Resume()
        {
            if (!IsRunning || !IsPaused) return;
            NtResumeProcess(process.Handle);
            IsPaused = false;
            Update();
        }

        public void Cancel()
        {
            if (!IsRunning) return;
            if (IsPaused) Resume();
            process.Kill();
        }

        private void Update()
        {
            int complete = CountCompleted(directory, stem);
            TimeSpan elapsed = DateTime.UtcNow - started;
            double perItem = complete > startCount ? elapsed.TotalSeconds / (complete - startCount) : 0;
            TimeSpan remaining = perItem > 0 ? TimeSpan.FromSeconds(Math.Max(0, requested - complete) * perItem) : TimeSpan.Zero;
            SaveState(directory, stem, new PermutationJobState
            {
                Requested = requested, Completed = complete, Seed = seed,
                StartedUtc = started, Status = IsPaused ? "Paused" : IsRunning ? "Running" : "Stopped"
            });
            EventHandler<PermutationProgressEventArgs> handler = ProgressChanged;
            if (handler != null)
                handler(this, new PermutationProgressEventArgs
                {
                    Completed = complete, Total = requested, Paused = IsPaused,
                    Message = IsPaused ? "Paused; partial results are safe" : complete + " of " + requested + " complete",
                    EstimatedRemaining = perItem <= 0 ? "calculating..." : remaining.ToString(@"hh\:mm\:ss")
                });
        }

        public static PermutationJobState LoadState(string directory, string stem)
        {
            string path = StatePath(directory, stem);
            if (!File.Exists(path)) return null;
            try
            {
                DataContractJsonSerializer serializer = new DataContractJsonSerializer(typeof(PermutationJobState));
                using (FileStream stream = File.OpenRead(path)) return (PermutationJobState)serializer.ReadObject(stream);
            }
            catch { return null; }
        }

        private static void SaveState(string directory, string stem, PermutationJobState state)
        {
            string path = StatePath(directory, stem);
            string temp = path + ".tmp";
            DataContractJsonSerializer serializer = new DataContractJsonSerializer(typeof(PermutationJobState));
            using (FileStream stream = File.Create(temp)) serializer.WriteObject(stream, state);
            if (File.Exists(path)) File.Replace(temp, path, path + ".bak", true);
            else File.Move(temp, path);
        }

        private static int CountCompleted(string directory, string stem)
        {
            string file = Directory.Exists(directory)
                ? Directory.GetFiles(directory, stem + ".z*e").OrderByDescending(File.GetLastWriteTimeUtc).FirstOrDefault()
                : null;
            if (string.IsNullOrEmpty(file)) return 0;
            return File.ReadLines(file).Count(line =>
            {
                string[] p = line.Trim().Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
                int repetition;
                return p.Length >= 2 && int.TryParse(p[0], out repetition);
            });
        }

        private static string StatePath(string directory, string stem)
        {
            return Path.Combine(directory, (string.IsNullOrWhiteSpace(stem) ? "qtlcart" : stem) + ".permutation-job.json");
        }

        [DllImport("ntdll.dll")]
        private static extern int NtSuspendProcess(IntPtr processHandle);

        [DllImport("ntdll.dll")]
        private static extern int NtResumeProcess(IntPtr processHandle);
    }
}
