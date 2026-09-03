using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Windows.Forms;

namespace QTLCartographer.Gui
{
    internal sealed class MainForm : Form
    {
        private static readonly Color Navy = Color.FromArgb(26, 43, 68);
        private static readonly Color Blue = Color.FromArgb(38, 116, 165);
        private static readonly Color Pale = Color.FromArgb(242, 246, 249);
        private static readonly Color Border = Color.FromArgb(211, 221, 229);

        private readonly List<ToolDefinition> tools = ToolCatalog.Create();
        private readonly Dictionary<string, ToolHelp> helpCache = new Dictionary<string, ToolHelp>();
        private readonly List<OptionEditor> editors = new List<OptionEditor>();
        private readonly List<CommandRequest> queue = new List<CommandRequest>();

        private TreeView toolTree;
        private TextBox searchBox;
        private Label titleLabel;
        private Label summaryLabel;
        private TextBox workingDirectoryBox;
        private TextBox stemBox;
        private TextBox resourceBox;
        private TextBox rawArgumentsBox;
        private CheckBox automaticBox;
        private CheckBox quietBox;
        private TableLayoutPanel optionTable;
        private RichTextBox outputBox;
        private TextBox commandPreview;
        private Button runButton;
        private Button cancelButton;
        private Button helpButton;
        private ListView queueView;
        private ListView filesView;
        private RichTextBox filePreview;
        private ComboBox fileFilter;
        private TabControl tabs;
        private ProgressBar workflowProgress;
        private ToolStripStatusLabel elapsedLabel;
        private Timer elapsedTimer;
        private string currentProjectFile;
        private readonly StringBuilder currentError = new StringBuilder();
        private ToolStripStatusLabel statusLabel;
        private ToolDefinition selectedTool;
        private Process currentProcess;
        private CommandRequest currentRequest;
        private bool runningQueue;
        private bool synchronizingProject;
        private int queueIndex;
        private string lastHelpText = "";

        private string ToolsDirectory
        {
            get { return Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "tools"); }
        }

        public MainForm()
        {
            Text = "QTL Cartographer " + ProductInfo.Version;
            Icon = SystemIcons.Application;
            MinimumSize = new Size(1180, 760);
            Size = new Size(1400, 900);
            StartPosition = FormStartPosition.CenterScreen;
            BackColor = Pale;
            Font = new Font("Segoe UI", 9F);
            AutoScaleMode = AutoScaleMode.Dpi;
            KeyPreview = true;

            BuildInterface();
            BuildApplicationMenu();
            PopulateToolTree("");
            if (toolTree.Nodes.Count > 0 && toolTree.Nodes[0].Nodes.Count > 0)
                toolTree.SelectedNode = toolTree.Nodes[0].Nodes[0];

            string projects = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments), "QTL Cartographer Projects");
            Directory.CreateDirectory(projects);
            workingDirectoryBox.Text = projects;

            FormClosing += delegate
            {
                if (currentProcess != null && !currentProcess.HasExited)
                    currentProcess.Kill();
            };
            elapsedTimer = new Timer { Interval = 500 };
            elapsedTimer.Tick += delegate { UpdateElapsedTime(); };
        }

        private void BuildInterface()
        {
            Panel banner = new Panel { Dock = DockStyle.Top, Height = 92, BackColor = Navy };
            Label appName = new Label
            {
                Text = "QTL Cartographer",
                ForeColor = Color.White,
                Font = new Font("Segoe UI Semibold", 21F),
                AutoSize = true,
                Location = new Point(24, 14)
            };
            Label subtitle = new Label
            {
                Text = "Quantitative trait locus analysis workspace",
                ForeColor = Color.FromArgb(190, 211, 227),
                AutoSize = false,
                Location = new Point(27, 57),
                Size = new Size(640, 24),
                TextAlign = ContentAlignment.MiddleLeft
            };
            banner.Controls.Add(appName);
            banner.Controls.Add(subtitle);
            Controls.Add(banner);

            StatusStrip status = new StatusStrip { SizingGrip = false };
            statusLabel = new ToolStripStatusLabel("Ready") { Spring = true, TextAlign = ContentAlignment.MiddleLeft };
            elapsedLabel = new ToolStripStatusLabel("Elapsed 00:00");
            workflowProgress = new ProgressBar { Width = 150, Height = 16 };
            ToolStripControlHost progressHost = new ToolStripControlHost(workflowProgress);
            status.Items.Add(statusLabel);
            status.Items.Add(elapsedLabel);
            status.Items.Add(progressHost);
            Controls.Add(status);

            SplitContainer split = new SplitContainer
            {
                Dock = DockStyle.None,
                FixedPanel = FixedPanel.Panel1,
                BackColor = Border
            };
            Controls.Add(split);
            Shown += delegate
            {
                PositionWorkspace(split, banner, status);
                split.Panel1MinSize = 280;
                split.Panel2MinSize = 760;
                split.SplitterDistance = 300;
            };
            Resize += delegate { PositionWorkspace(split, banner, status); };

            TableLayoutPanel navigation = new TableLayoutPanel
            {
                Dock = DockStyle.Fill,
                BackColor = Color.White,
                Padding = new Padding(14),
                ColumnCount = 1,
                RowCount = 2
            };
            navigation.ColumnStyles.Clear();
            navigation.RowStyles.Clear();
            navigation.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 100));
            navigation.RowStyles.Add(new RowStyle(SizeType.Absolute, 68));
            navigation.RowStyles.Add(new RowStyle(SizeType.Percent, 100));
            Panel navigationHeader = new Panel { Dock = DockStyle.Top, Height = 68, BackColor = Color.White };
            Label toolsLabel = new Label
            {
                Text = "ANALYSIS TOOLS",
                ForeColor = Color.FromArgb(90, 105, 118),
                Font = new Font("Segoe UI Semibold", 8F),
                AutoSize = true,
                Location = new Point(0, 2)
            };
            searchBox = new TextBox { Anchor = AnchorStyles.Top | AnchorStyles.Left | AnchorStyles.Right, Location = new Point(0, 27), Width = 220 };
            searchBox.TextChanged += delegate { PopulateToolTree(searchBox.Text); };
            navigationHeader.Resize += delegate { searchBox.Width = navigationHeader.ClientSize.Width; };
            navigationHeader.Controls.Add(toolsLabel);
            navigationHeader.Controls.Add(searchBox);
            toolTree = new TreeView
            {
                Dock = DockStyle.Fill,
                BorderStyle = BorderStyle.None,
                FullRowSelect = true,
                HideSelection = false,
                ItemHeight = 25,
                ShowLines = false,
                ShowPlusMinus = false,
                ShowRootLines = false,
                BackColor = Color.White
            };
            toolTree.AfterSelect += ToolTreeAfterSelect;
            navigation.Controls.Add(navigationHeader, 0, 0);
            navigation.Controls.Add(toolTree, 0, 1);
            split.Panel1.Controls.Add(navigation);

            TableLayoutPanel content = new TableLayoutPanel
            {
                Dock = DockStyle.Fill,
                BackColor = Pale,
                Padding = new Padding(18),
                ColumnCount = 1,
                RowCount = 2
            };
            content.ColumnStyles.Clear();
            content.RowStyles.Clear();
            content.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 100));
            content.RowStyles.Add(new RowStyle(SizeType.Absolute, 74));
            content.RowStyles.Add(new RowStyle(SizeType.Percent, 100));
            Panel heading = new Panel { Dock = DockStyle.Fill, BackColor = Pale };
            titleLabel = new Label
            {
                Text = "Choose an analysis tool",
                ForeColor = Navy,
                Font = new Font("Segoe UI Semibold", 18F),
                AutoSize = true,
                Location = new Point(2, 0)
            };
            summaryLabel = new Label
            {
                Text = "Select a program on the left to configure its complete command-line interface.",
                ForeColor = Color.FromArgb(80, 92, 104),
                AutoEllipsis = true,
                Location = new Point(5, 40),
                Size = new Size(730, 25)
            };
            helpButton = MakeButton("Program help", false);
            helpButton.Anchor = AnchorStyles.Top | AnchorStyles.Right;
            helpButton.Location = new Point(760, 9);
            helpButton.Size = new Size(115, 34);
            helpButton.Enabled = false;
            helpButton.Click += delegate { ShowProgramHelp(); };
            heading.Resize += delegate { helpButton.Left = heading.ClientSize.Width - helpButton.Width; };
            heading.Controls.Add(titleLabel);
            heading.Controls.Add(summaryLabel);
            heading.Controls.Add(helpButton);

            tabs = new TabControl
            {
                Dock = DockStyle.Fill,
                Multiline = true,
                SizeMode = TabSizeMode.Normal,
                Padding = new Point(14, 5)
            };
            tabs.TabPages.Add(BuildConfigurePage());
            tabs.TabPages.Add(BuildOutputPage());
            tabs.TabPages.Add(BuildQueuePage());
            tabs.TabPages.Add(BuildFilesPage());
            tabs.TabPages.Add(BuildAboutPage());
            content.Controls.Add(heading, 0, 0);
            content.Controls.Add(tabs, 0, 1);
            split.Panel2.Controls.Add(content);
        }

        private void PositionWorkspace(SplitContainer split, Panel banner, StatusStrip status)
        {
            int top = Math.Max(banner.Bottom, MainMenuStrip == null ? 0 : MainMenuStrip.Bottom);
            int bottom = status.Top;
            split.SetBounds(0, top, ClientSize.Width, Math.Max(100, bottom - top));
        }

        private void BuildApplicationMenu()
        {
            MenuStrip menu = new MenuStrip { Dock = DockStyle.Top };
            ToolStripMenuItem file = new ToolStripMenuItem("&File");
            ((ToolStripMenuItem)file.DropDownItems.Add("&New project…", null, delegate { NewProject(); })).ShortcutKeys = Keys.Control | Keys.N;
            ((ToolStripMenuItem)file.DropDownItems.Add("&Open project…", null, delegate { OpenProject(); })).ShortcutKeys = Keys.Control | Keys.O;
            ((ToolStripMenuItem)file.DropDownItems.Add("&Save project", null, delegate { SaveProject(false); })).ShortcutKeys = Keys.Control | Keys.S;
            file.DropDownItems.Add("Save project &as…", null, delegate { SaveProject(true); });
            file.DropDownItems.Add("Export &reproducibility bundle…", null, delegate { ExportReproducibilityBundle(); });
            file.DropDownItems.Add(new ToolStripSeparator());
            file.DropDownItems.Add("E&xit", null, delegate { Close(); });
            ToolStripMenuItem analysis = new ToolStripMenuItem("&Analysis");
            ((ToolStripMenuItem)analysis.DropDownItems.Add("&Results dashboard…", null, delegate { ShowResultsDashboard(); })).ShortcutKeys = Keys.Control | Keys.D;
            analysis.DropDownItems.Add("&Load example analysis", null, delegate { LoadSampleWorkflow(); });
            analysis.DropDownItems.Add("&Pre-analysis diagnostics…", null, delegate { ShowDiagnostics(); });
            analysis.DropDownItems.Add("Covariates and analysis &design…", null, delegate { ShowAnalysisDesign(); });
            analysis.DropDownItems.Add("&Scientific benchmarks…", null, delegate { ShowBenchmarks(); });
            ToolStripMenuItem data = new ToolStripMenuItem("&Data");
            data.DropDownItems.Add("&Modern format import/export…", null, delegate { ShowModernFormats(); });
            ToolStripMenuItem help = new ToolStripMenuItem("&Help");
            ((ToolStripMenuItem)help.DropDownItems.Add("Selected program help", null, delegate { ShowProgramHelp(); })).ShortcutKeys = Keys.F1;
            help.DropDownItems.Add("&Interactive example tutorial…", null, delegate { ShowTutorial(); });
            menu.Items.Add(file);
            menu.Items.Add(analysis);
            menu.Items.Add(data);
            menu.Items.Add(help);
            MainMenuStrip = menu;
            Controls.Add(menu);
            menu.BringToFront();
        }

        private TabPage BuildConfigurePage()
        {
            TabPage page = new TabPage("Configure") { BackColor = Color.White, Padding = new Padding(14) };
            Panel actions = new Panel { Dock = DockStyle.Bottom, Height = 92, BackColor = Color.White };

            commandPreview = new TextBox
            {
                Dock = DockStyle.Top,
                ReadOnly = true,
                BackColor = Color.FromArgb(247, 249, 251),
                ForeColor = Navy,
                Height = 28
            };
            runButton = MakeButton("Run analysis", true);
            runButton.Enabled = false;
            runButton.Location = new Point(0, 43);
            runButton.Size = new Size(128, 36);
            runButton.Click += delegate { RunCurrent(); };
            Button addQueue = MakeButton("Add to queue", false);
            addQueue.Location = new Point(138, 43);
            addQueue.Size = new Size(120, 36);
            addQueue.Click += delegate { AddCurrentToQueue(); };
            actions.Controls.Add(commandPreview);
            actions.Controls.Add(runButton);
            actions.Controls.Add(addQueue);

            Panel scroll = new Panel { Dock = DockStyle.Fill, AutoScroll = true, BackColor = Color.White };
            TableLayoutPanel settings = new TableLayoutPanel
            {
                Dock = DockStyle.Top,
                AutoSize = true,
                ColumnCount = 3,
                RowCount = 6,
                Padding = new Padding(2),
                BackColor = Color.White
            };
            settings.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 155));
            settings.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 100));
            settings.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 92));

            workingDirectoryBox = AddPathRow(settings, 0, "Working directory", true);
            stemBox = AddTextRow(settings, 1, "Filename stem", "", "Optional -X value");
            resourceBox = AddPathRow(settings, 2, "Resource file", false);
            automaticBox = new CheckBox { Text = "Automatic mode (-A)", Checked = true, AutoSize = true, Margin = new Padding(3, 8, 3, 8) };
            quietBox = new CheckBox { Text = "Non-verbose (-V)", AutoSize = true, Margin = new Padding(18, 8, 3, 8) };
            FlowLayoutPanel commonFlags = new FlowLayoutPanel { Dock = DockStyle.Fill, AutoSize = true };
            commonFlags.Controls.Add(automaticBox);
            commonFlags.Controls.Add(quietBox);
            settings.Controls.Add(MakeFieldLabel("Execution"), 0, 3);
            settings.Controls.Add(commonFlags, 1, 3);
            settings.SetColumnSpan(commonFlags, 2);

            Label optionsHeader = new Label
            {
                Text = "PROGRAM OPTIONS",
                ForeColor = Color.FromArgb(90, 105, 118),
                Font = new Font("Segoe UI Semibold", 8F),
                AutoSize = true,
                Margin = new Padding(3, 18, 3, 8)
            };
            settings.Controls.Add(optionsHeader, 0, 4);
            settings.SetColumnSpan(optionsHeader, 3);

            optionTable = new TableLayoutPanel
            {
                Dock = DockStyle.Top,
                AutoSize = true,
                ColumnCount = 4,
                BackColor = Color.White,
                Margin = new Padding(0)
            };
            optionTable.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 84));
            optionTable.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 190));
            optionTable.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 36));
            optionTable.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 100));
            settings.Controls.Add(optionTable, 0, 5);
            settings.SetColumnSpan(optionTable, 3);

            Label rawLabel = MakeFieldLabel("Additional arguments");
            rawLabel.Margin = new Padding(3, 18, 3, 3);
            rawArgumentsBox = new TextBox { Dock = DockStyle.Fill, Margin = new Padding(3, 14, 3, 6) };
            settings.RowCount = 7;
            settings.Controls.Add(rawLabel, 0, 6);
            settings.Controls.Add(rawArgumentsBox, 1, 6);
            settings.SetColumnSpan(rawArgumentsBox, 2);

            EventHandler update = delegate { UpdateCommandPreview(); };
            workingDirectoryBox.TextChanged += delegate
            {
                UpdateCommandPreview();
                SynchronizeProjectState(null);
            };
            stemBox.TextChanged += update;
            resourceBox.TextChanged += update;
            rawArgumentsBox.TextChanged += update;
            automaticBox.CheckedChanged += update;
            quietBox.CheckedChanged += update;

            scroll.Controls.Add(settings);
            page.Controls.Add(scroll);
            page.Controls.Add(actions);
            return page;
        }

        private TabPage BuildOutputPage()
        {
            TabPage page = new TabPage("Console") { BackColor = Color.White, Padding = new Padding(10) };
            FlowLayoutPanel toolbar = new FlowLayoutPanel { Dock = DockStyle.Top, Height = 42 };
            cancelButton = MakeButton("Cancel", false);
            cancelButton.Enabled = false;
            cancelButton.Click += delegate { CancelCurrentProcess(); };
            Button clear = MakeButton("Clear", false);
            clear.Click += delegate { outputBox.Clear(); };
            Button save = MakeButton("Save log", false);
            save.Click += delegate { SaveConsoleLog(); };
            toolbar.Controls.Add(cancelButton);
            toolbar.Controls.Add(clear);
            toolbar.Controls.Add(save);
            outputBox = new RichTextBox
            {
                Dock = DockStyle.Fill,
                ReadOnly = true,
                BackColor = Color.FromArgb(20, 29, 38),
                ForeColor = Color.FromArgb(220, 232, 239),
                Font = new Font("Consolas", 9.5F),
                BorderStyle = BorderStyle.None,
                DetectUrls = false
            };
            page.Controls.Add(outputBox);
            page.Controls.Add(toolbar);
            return page;
        }

        private TabPage BuildQueuePage()
        {
            TabPage page = new TabPage("Workflow queue") { BackColor = Color.White, Padding = new Padding(10) };
            FlowLayoutPanel toolbar = new FlowLayoutPanel { Dock = DockStyle.Top, Height = 44 };
            Button run = MakeButton("Run queue", true);
            run.Click += delegate { RunQueue(); };
            Button sample = MakeButton("Load example analysis", false);
            sample.Click += delegate { LoadSampleWorkflow(); };
            sample.AutoSize = true;
            Button remove = MakeButton("Remove selected", false);
            remove.Click += delegate { RemoveSelectedQueueItems(); };
            Button clear = MakeButton("Clear", false);
            clear.Click += delegate { queue.Clear(); RefreshQueue(); };
            Button retry = MakeButton("Retry failed stage", false);
            retry.AutoSize = true;
            retry.Click += delegate { RetryFailedStage(); };
            toolbar.Controls.Add(run);
            toolbar.Controls.Add(sample);
            toolbar.Controls.Add(remove);
            toolbar.Controls.Add(retry);
            toolbar.Controls.Add(clear);

            queueView = new ListView
            {
                Dock = DockStyle.Fill,
                View = View.Details,
                FullRowSelect = true,
                GridLines = true
            };
            queueView.Columns.Add("#", 45);
            queueView.Columns.Add("Program", 125);
            queueView.Columns.Add("Status", 95);
            queueView.Columns.Add("Elapsed", 75);
            queueView.Columns.Add("Arguments", 410);
            queueView.Columns.Add("Working directory", 220);
            page.Controls.Add(queueView);
            page.Controls.Add(toolbar);
            return page;
        }

        private TabPage BuildFilesPage()
        {
            TabPage page = new TabPage("Project files") { BackColor = Color.White, Padding = new Padding(10) };
            FlowLayoutPanel toolbar = new FlowLayoutPanel { Dock = DockStyle.Top, Height = 44 };
            Button refresh = MakeButton("Refresh", false);
            refresh.Click += delegate { RefreshFiles(); };
            Button open = MakeButton("Open selected", false);
            open.Click += delegate { OpenSelectedFile(); };
            Button folder = MakeButton("Open folder", false);
            folder.Click += delegate { OpenWorkingDirectory(); };
            fileFilter = new ComboBox { DropDownStyle = ComboBoxStyle.DropDownList, Width = 110, AccessibleName = "Project file filter" };
            fileFilter.Items.AddRange(new object[] { "All files", "Inputs", "Results", "Logs", "Plots" });
            fileFilter.SelectedIndex = 0;
            fileFilter.SelectedIndexChanged += delegate { RefreshFiles(); };
            toolbar.Controls.Add(refresh);
            toolbar.Controls.Add(open);
            toolbar.Controls.Add(folder);
            toolbar.Controls.Add(new Label { Text = "Filter:", AutoSize = true, Margin = new Padding(12, 9, 3, 3) });
            toolbar.Controls.Add(fileFilter);

            SplitContainer split = new SplitContainer { Dock = DockStyle.Fill, SplitterDistance = 610 };
            filesView = new ListView
            {
                Dock = DockStyle.Fill,
                View = View.Details,
                FullRowSelect = true,
                GridLines = true
            };
            filesView.Columns.Add("Name", 240);
            filesView.Columns.Add("Type / purpose", 245);
            filesView.Columns.Add("Size", 80);
            filesView.Columns.Add("Modified", 135);
            filesView.DoubleClick += delegate { PreviewSelectedFile(); };
            filesView.SelectedIndexChanged += delegate { PreviewSelectedFile(); };
            filePreview = new RichTextBox
            {
                Dock = DockStyle.Fill, ReadOnly = true, WordWrap = false,
                Font = new Font("Consolas", 9F), BackColor = Color.White,
                AccessibleName = "Selected project file preview"
            };
            split.Panel1.Controls.Add(filesView);
            split.Panel2.Controls.Add(filePreview);
            page.Controls.Add(split);
            page.Controls.Add(toolbar);
            return page;
        }

        private TabPage BuildAboutPage()
        {
            TabPage page = new TabPage("About") { BackColor = Color.White, Padding = new Padding(28) };
            Label about = new Label
            {
                Dock = DockStyle.Top,
                Height = 240,
                Font = new Font("Segoe UI", 10F),
                Text = "QTL Cartographer for Windows " + ProductInfo.Version + "\r\n\r\n" +
                       "A native Windows desktop interface for QTL Cartographer 1.17. " +
                       "All statistical calculations are performed by the original GPL-licensed C engine.\r\n\r\n" +
                       "The GUI includes guided import, typed options, recoverable workflows, project persistence, " +
                       "integrated LR/LOD charts, empirical permutation thresholds, peak tables, and report export.\r\n\r\n" +
                       "QTL Cartographer authors: C. J. Basten, B. S. Weir, and Z.-B. Zeng."
            };
            LinkLabel source = new LinkLabel
            {
                Text = "Original QTL Cartographer source",
                AutoSize = true,
                Top = 250,
                Left = 28
            };
            source.Click += delegate { Process.Start("https://github.com/cbasten/qtlcart"); };
            page.Controls.Add(source);
            page.Controls.Add(about);
            return page;
        }

        private TextBox AddPathRow(TableLayoutPanel table, int row, string label, bool folder)
        {
            TextBox box = new TextBox { Dock = DockStyle.Fill, Margin = new Padding(3, 5, 3, 5) };
            Button browse = MakeButton("Browse...", false);
            browse.Dock = DockStyle.Fill;
            browse.Margin = new Padding(3, 3, 3, 3);
            browse.Click += delegate
            {
                if (folder)
                {
                    using (FolderBrowserDialog dialog = new FolderBrowserDialog())
                    {
                        dialog.SelectedPath = Directory.Exists(box.Text) ? box.Text : "";
                        if (dialog.ShowDialog(this) == DialogResult.OK)
                            box.Text = dialog.SelectedPath;
                    }
                }
                else
                {
                    using (OpenFileDialog dialog = new OpenFileDialog())
                    {
                        dialog.Filter = "All files (*.*)|*.*";
                        if (dialog.ShowDialog(this) == DialogResult.OK)
                            box.Text = dialog.FileName;
                    }
                }
            };
            table.Controls.Add(MakeFieldLabel(label), 0, row);
            table.Controls.Add(box, 1, row);
            table.Controls.Add(browse, 2, row);
            return box;
        }

        private TextBox AddTextRow(TableLayoutPanel table, int row, string label, string value, string hint)
        {
            TextBox box = new TextBox { Text = value, Dock = DockStyle.Fill, Margin = new Padding(3, 5, 3, 5) };
            Label note = new Label { Text = hint, Dock = DockStyle.Fill, ForeColor = Color.Gray, TextAlign = ContentAlignment.MiddleLeft };
            table.Controls.Add(MakeFieldLabel(label), 0, row);
            table.Controls.Add(box, 1, row);
            table.Controls.Add(note, 2, row);
            return box;
        }

        private Label MakeFieldLabel(string text)
        {
            return new Label { Text = text, Dock = DockStyle.Fill, TextAlign = ContentAlignment.MiddleLeft, ForeColor = Navy };
        }

        private Button MakeButton(string text, bool primary)
        {
            return new Button
            {
                Text = text,
                AutoSize = false,
                Size = new Size(110, 32),
                FlatStyle = FlatStyle.Flat,
                BackColor = primary ? Blue : Color.White,
                ForeColor = primary ? Color.White : Navy,
                Cursor = Cursors.Hand,
                Margin = new Padding(3, 3, 6, 3),
                FlatAppearance = { BorderColor = primary ? Blue : Border }
            };
        }

        private void PopulateToolTree(string filter)
        {
            string selectedName = selectedTool == null ? "" : selectedTool.Name;
            toolTree.BeginUpdate();
            toolTree.Nodes.Clear();
            Dictionary<string, TreeNode> categories = new Dictionary<string, TreeNode>();
            foreach (ToolDefinition tool in tools)
            {
                if (filter.Length > 0 && tool.Name.IndexOf(filter, StringComparison.OrdinalIgnoreCase) < 0 &&
                    tool.Summary.IndexOf(filter, StringComparison.OrdinalIgnoreCase) < 0)
                    continue;

                TreeNode category;
                if (!categories.TryGetValue(tool.Category, out category))
                {
                    category = new TreeNode(tool.Category) { ForeColor = Color.FromArgb(100, 112, 122), NodeFont = new Font(Font, FontStyle.Bold) };
                    categories.Add(tool.Category, category);
                    toolTree.Nodes.Add(category);
                }
                TreeNode node = new TreeNode(tool.Name) { Tag = tool, ForeColor = Navy };
                category.Nodes.Add(node);
                if (tool.Name == selectedName)
                    toolTree.SelectedNode = node;
            }
            toolTree.ExpandAll();
            toolTree.EndUpdate();
        }

        private void ToolTreeAfterSelect(object sender, TreeViewEventArgs e)
        {
            ToolDefinition tool = e.Node.Tag as ToolDefinition;
            if (tool == null)
                return;

            selectedTool = tool;
            titleLabel.Text = tool.Name;
            summaryLabel.Text = tool.Summary;
            helpButton.Enabled = true;
            runButton.Enabled = true;
            LoadOptions(tool);
        }

        private void LoadOptions(ToolDefinition tool)
        {
            optionTable.SuspendLayout();
            optionTable.Controls.Clear();
            optionTable.RowStyles.Clear();
            optionTable.RowCount = 0;
            editors.Clear();

            try
            {
                ToolHelp help;
                if (!helpCache.TryGetValue(tool.Name, out help))
                {
                    string executable = Path.Combine(ToolsDirectory, tool.Name + ".exe");
                    if (!File.Exists(executable))
                        throw new FileNotFoundException("The native analysis program is missing.", executable);
                    help = HelpParser.Load(executable);
                    helpCache.Add(tool.Name, help);
                }
                lastHelpText = help.FullText;
                if (!string.IsNullOrEmpty(help.Purpose))
                    summaryLabel.Text = help.Purpose;

                int row = 0;
                foreach (OptionDefinition option in help.Options)
                {
                    optionTable.RowCount++;
                    optionTable.RowStyles.Add(new RowStyle(SizeType.AutoSize));
                    CheckBox enabled = new CheckBox
                    {
                        Text = option.Flag,
                        AutoSize = true,
                        Font = new Font("Consolas", 9F, FontStyle.Bold),
                        Margin = new Padding(4, 8, 2, 5)
                    };
                    Control value = CreateOptionControl(option);
                    Button browse = MakeButton("…", false);
                    browse.Dock = DockStyle.Fill;
                    browse.Margin = new Padding(2, 4, 2, 4);
                    browse.Visible = LooksLikeFileOption(option);
                    Label description = new Label
                    {
                        Text = option.Description,
                        Dock = DockStyle.Fill,
                        AutoEllipsis = true,
                        TextAlign = ContentAlignment.MiddleLeft,
                        ForeColor = Color.FromArgb(70, 79, 88),
                        Margin = new Padding(8, 6, 3, 4)
                    };
                    OptionEditor editor = new OptionEditor { Definition = option, Enabled = enabled, Value = value };
                    editors.Add(editor);
                    enabled.CheckedChanged += delegate { value.Enabled = enabled.Checked; UpdateCommandPreview(); };
                    value.TextChanged += delegate
                    {
                        InferStemFromOutput(editor);
                        UpdateCommandPreview();
                    };
                    browse.Click += delegate { BrowseOption(editor); };
                    optionTable.Controls.Add(enabled, 0, row);
                    optionTable.Controls.Add(value, 1, row);
                    optionTable.Controls.Add(browse, 2, row);
                    optionTable.Controls.Add(description, 3, row);
                    row++;
                }
            }
            catch (Exception ex)
            {
                Label error = new Label { AutoSize = true, ForeColor = Color.Firebrick, Text = ex.Message };
                optionTable.Controls.Add(error, 0, 0);
                optionTable.SetColumnSpan(error, 4);
                runButton.Enabled = false;
            }
            finally
            {
                optionTable.ResumeLayout();
                SynchronizeProjectState(null);
                UpdateCommandPreview();
            }
        }

        private Control CreateOptionControl(OptionDefinition option)
        {
            string description = option.Description ?? "";
            if (!LooksLikeFileOption(option))
            {
                Match mapped = System.Text.RegularExpressions.Regex.Match(description, @"=>\s*\(([^)]+)\)");
                if (mapped.Success)
                {
                    ComboBox choice = new ComboBox
                    {
                        Dock = DockStyle.Fill, Enabled = false, DropDownStyle = ComboBoxStyle.DropDown,
                        Margin = new Padding(3, 5, 3, 5), AccessibleName = description
                    };
                    foreach (string item in mapped.Groups[1].Value.Split(','))
                        choice.Items.Add(item.Trim());
                    choice.Text = option.DefaultValue;
                    return choice;
                }
                decimal numeric;
                if (decimal.TryParse(option.DefaultValue, System.Globalization.NumberStyles.Float,
                    System.Globalization.CultureInfo.InvariantCulture, out numeric))
                {
                    decimal minimum = -2147483648;
                    decimal maximum = 2147483647;
                    Match range = Regex.Match(description, @"\[\s*(-?\d+(?:\.\d+)?)\s*[-,]\s*(-?\d+(?:\.\d+)?)\s*\]");
                    decimal parsedMinimum;
                    decimal parsedMaximum;
                    if (range.Success &&
                        decimal.TryParse(range.Groups[1].Value, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out parsedMinimum) &&
                        decimal.TryParse(range.Groups[2].Value, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out parsedMaximum))
                    {
                        minimum = parsedMinimum;
                        maximum = parsedMaximum;
                    }
                    NumericUpDown number = new NumericUpDown
                    {
                        Dock = DockStyle.Fill, Enabled = false, Minimum = minimum,
                        Maximum = maximum, DecimalPlaces = option.DefaultValue.Contains(".") ? 6 : 0,
                        Increment = option.DefaultValue.Contains(".") ? .1M : 1M,
                        Margin = new Padding(3, 5, 3, 5), AccessibleName = description,
                        Value = Math.Max(minimum, Math.Min(maximum, numeric))
                    };
                    return number;
                }
            }
            TextBox text = new TextBox
            {
                Text = option.DefaultValue, Dock = DockStyle.Fill, Enabled = false,
                Margin = new Padding(3, 5, 3, 5), AccessibleName = description
            };
            text.Validating += delegate(object sender, System.ComponentModel.CancelEventArgs e)
            {
                text.BackColor = string.IsNullOrWhiteSpace(text.Text) ? Color.MistyRose : Color.White;
            };
            return text;
        }

        private bool LooksLikeFileOption(OptionDefinition option)
        {
            string value = option.Description.ToLowerInvariant();
            return value.Contains("file") || value.Contains("input") || value.Contains("output");
        }

        private void BrowseOption(OptionEditor editor)
        {
            string description = editor.Definition.Description.ToLowerInvariant();
            bool save = description.Contains("output") || description.Contains("error") || description.Contains("results");
            FileDialog dialog = save ? (FileDialog)new SaveFileDialog() : new OpenFileDialog();
            using (dialog)
            {
                dialog.InitialDirectory = Directory.Exists(workingDirectoryBox.Text) ? workingDirectoryBox.Text : "";
                dialog.Filter = "All files (*.*)|*.*";
                if (dialog.ShowDialog(this) == DialogResult.OK)
                {
                    string value = dialog.FileName;
                    if (string.Equals(Path.GetDirectoryName(value), workingDirectoryBox.Text, StringComparison.OrdinalIgnoreCase))
                        value = Path.GetFileName(value);
                    editor.Value.Text = value;
                    editor.Enabled.Checked = true;
                }
            }
        }

        private CommandRequest BuildCurrentRequest()
        {
            if (selectedTool == null)
                return null;
            List<string> arguments = new List<string>();
            Dictionary<string, string> options = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            if (resourceBox.Text.Trim().Length > 0)
            {
                arguments.Add("-R");
                arguments.Add(Quote(resourceBox.Text.Trim()));
            }
            if (stemBox.Text.Trim().Length > 0)
            {
                arguments.Add("-X");
                arguments.Add(Quote(stemBox.Text.Trim()));
            }
            foreach (OptionEditor editor in editors)
            {
                if (!editor.Enabled.Checked)
                    continue;
                arguments.Add(editor.Definition.Flag);
                string optionValue = editor.Value.Text.Trim();
                options[editor.Definition.Flag] = optionValue;
                if (optionValue.Length > 0)
                    arguments.Add(Quote(optionValue));
            }
            if (automaticBox.Checked)
                arguments.Add("-A");
            if (quietBox.Checked)
                arguments.Add("-V");
            if (rawArgumentsBox.Text.Trim().Length > 0)
                arguments.Add(rawArgumentsBox.Text.Trim());

            return new CommandRequest
            {
                Tool = selectedTool,
                Arguments = string.Join(" ", arguments.ToArray()),
                WorkingDirectory = workingDirectoryBox.Text.Trim(),
                ResourceFile = resourceBox.Text.Trim(),
                RequestedStem = stemBox.Text.Trim(),
                Options = options
            };
        }

        private static string Quote(string value)
        {
            if (value.Length == 0)
                return "\"\"";
            if (value.IndexOfAny(new[] { ' ', '\t', '\"' }) < 0)
                return value;
            return "\"" + value.Replace("\"", "\\\"") + "\"";
        }

        private void UpdateCommandPreview()
        {
            if (commandPreview == null)
                return;
            CommandRequest request = BuildCurrentRequest();
            commandPreview.Text = request == null ? "" : request.DisplayCommand;
        }

        private void RunCurrent()
        {
            CommandRequest request = BuildCurrentRequest();
            if (request != null)
                StartRequest(request, false);
        }

        private void StartRequest(CommandRequest request, bool fromQueue)
        {
            if (currentProcess != null && !currentProcess.HasExited)
            {
                MessageBox.Show(this, "An analysis is already running.", Text, MessageBoxButtons.OK, MessageBoxIcon.Information);
                return;
            }
            if (request.WorkingDirectory.Length == 0)
            {
                MessageBox.Show(this, "Choose a working directory first.", Text, MessageBoxButtons.OK, MessageBoxIcon.Warning);
                return;
            }
            Directory.CreateDirectory(request.WorkingDirectory);
            if (!ConfirmRcrossMap(request))
                return;
            if (!ConfirmOverwrite(request))
                return;
            string executable = Path.Combine(ToolsDirectory, request.Tool.Name + ".exe");
            if (!File.Exists(executable))
            {
                MessageBox.Show(this, "Missing analysis program:\r\n" + executable, Text, MessageBoxButtons.OK, MessageBoxIcon.Error);
                return;
            }

            AppendOutput("\r\n> " + request.DisplayCommand + "\r\n", Color.FromArgb(100, 190, 235));
            currentError.Clear();
            request.Status = "Running";
            request.StartedAt = DateTime.Now;
            request.FailureDetails = "";
            RefreshQueue();
            elapsedTimer.Start();
            UpdateWorkflowProgress();
            Process process = new Process();
            process.StartInfo = new ProcessStartInfo(executable, request.Arguments)
            {
                WorkingDirectory = request.WorkingDirectory,
                UseShellExecute = false,
                CreateNoWindow = true,
                RedirectStandardOutput = true,
                RedirectStandardError = true
            };
            process.EnableRaisingEvents = true;
            process.OutputDataReceived += delegate(object sender, DataReceivedEventArgs e)
            {
                if (e.Data != null) AppendOutput(e.Data + Environment.NewLine, null);
            };
            process.ErrorDataReceived += delegate(object sender, DataReceivedEventArgs e)
            {
                if (e.Data != null)
                {
                    lock (currentError) { currentError.AppendLine(e.Data); }
                    AppendOutput(e.Data + Environment.NewLine, Color.FromArgb(255, 150, 140));
                }
            };
            process.Exited += delegate
            {
                process.WaitForExit();
                int code = process.ExitCode;
                BeginInvoke((MethodInvoker)delegate
                {
                    AppendOutput("[Exited with code " + code + "]\r\n", code == 0 ? Color.FromArgb(125, 210, 145) : Color.FromArgb(255, 150, 140));
                    currentProcess = null;
                    request.Elapsed = DateTime.Now - request.StartedAt;
                    currentRequest = null;
                    request.Status = code == 0 ? "Succeeded" : "Failed";
                    request.FailureDetails = currentError.ToString();
                    runButton.Enabled = selectedTool != null;
                    cancelButton.Enabled = false;
                    statusLabel.Text = code == 0 ? "Completed " + request.Tool.Name : request.Tool.Name + " failed (exit " + code + ")";
                    RefreshFiles();
                    RefreshQueue();
                    UpdateWorkflowProgress();
                    elapsedTimer.Stop();
                    if (code == 0)
                        SynchronizeProjectState(request);
                    if (fromQueue && runningQueue)
                    {
                        if (code == 0)
                        {
                            queueIndex++;
                            RunNextQueueItem();
                        }
                        else
                        {
                            runningQueue = false;
                            ShowDetailedFailure(request, code);
                        }
                    }
                });
            };

            try
            {
                currentProcess = process;
                currentRequest = request;
                runButton.Enabled = false;
                cancelButton.Enabled = true;
                statusLabel.Text = "Running " + request.Tool.Name + "...";
                process.Start();
                process.BeginOutputReadLine();
                process.BeginErrorReadLine();
            }
            catch (Exception ex)
            {
                currentProcess = null;
                currentRequest = null;
                request.Status = "Failed";
                request.FailureDetails = ex.Message;
                runButton.Enabled = true;
                cancelButton.Enabled = false;
                MessageBox.Show(this, ex.Message, Text, MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
        }

        private bool ConfirmOverwrite(CommandRequest request)
        {
            if (request.OverwriteConfirmed)
                return true;
            List<string> outputs = new List<string>();
            string output;
            if (request.Options != null && request.Options.TryGetValue("-o", out output) && !string.IsNullOrWhiteSpace(output))
                outputs.Add(ResolveProjectPath(request.WorkingDirectory, output));
            else if (!string.IsNullOrWhiteSpace(request.RequestedStem))
            {
                Dictionary<string, string> extensions = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase)
                {
                    { "Rmap", ".map" }, { "Rcross", ".cro" }, { "Qstats", ".qst" },
                    { "LRmapqtl", ".lr" }, { "SRmapqtl", ".sr" }, { "Zmapqtl", ".z" },
                    { "MImapqtl", ".mim" }, { "Eqtl", ".eqt" }, { "Preplot", ".plt" }
                };
                string extension;
                if (extensions.TryGetValue(request.Tool.Name, out extension))
                    outputs.Add(Path.Combine(request.WorkingDirectory, request.RequestedStem + extension));
            }
            List<string> existing = outputs.Where(File.Exists).ToList();
            if (existing.Count == 0)
                return true;
            DialogResult result = MessageBox.Show(this,
                "This analysis will overwrite:\r\n" + string.Join("\r\n", existing.ToArray()) +
                "\r\n\r\nContinue?", "Confirm overwrite", MessageBoxButtons.YesNo, MessageBoxIcon.Warning,
                MessageBoxDefaultButton.Button2);
            request.OverwriteConfirmed = result == DialogResult.Yes;
            return request.OverwriteConfirmed;
        }

        private void ShowDetailedFailure(CommandRequest request, int code)
        {
            string relevant = request.FailureDetails;
            if (relevant.Length > 1800) relevant = relevant.Substring(relevant.Length - 1800);
            MessageBox.Show(this,
                "Workflow stopped at " + request.Tool.Name + " (exit " + code + ").\r\n\r\nCommand:\r\n" +
                request.DisplayCommand + "\r\n\r\nRelevant output:\r\n" + relevant +
                "\r\nSuggested remedy: check that input and map files exist, validate the project inputs, then use Retry failed stage.",
                "Analysis failed", MessageBoxButtons.OK, MessageBoxIcon.Error);
        }

        private void InferStemFromOutput(OptionEditor editor)
        {
            if (synchronizingProject || selectedTool == null || editor.Definition.Flag != "-o")
                return;
            if (selectedTool.Name != "Rmap" && selectedTool.Name != "Rcross")
                return;
            string stem = ProjectState.StemFromFile(editor.Value.Text);
            if (!string.IsNullOrEmpty(stem) && !string.Equals(stemBox.Text, stem, StringComparison.OrdinalIgnoreCase))
                stemBox.Text = stem;
        }

        private bool ConfirmRcrossMap(CommandRequest request)
        {
            if (request.Tool.Name != "Rcross")
                return true;
            string input;
            if (request.Options == null || !request.Options.TryGetValue("-i", out input) || string.IsNullOrWhiteSpace(input))
                return true;

            string map;
            if (request.Options.TryGetValue("-m", out map) && !string.IsNullOrWhiteSpace(map))
                map = ResolveProjectPath(request.WorkingDirectory, map);
            else if (!string.IsNullOrWhiteSpace(request.RequestedStem))
                map = ResolveProjectPath(request.WorkingDirectory, request.RequestedStem + ".map");
            else
            {
                ProjectState state = ProjectState.Load(request.WorkingDirectory, request.ResourceFile);
                map = state.ResolveProjectFile(request.WorkingDirectory, "-map", "");
            }
            if (!string.IsNullOrWhiteSpace(map) && File.Exists(map))
                return true;

            DialogResult result = MessageBox.Show(this,
                "The linkage map for this cross was not found:\r\n" +
                (string.IsNullOrWhiteSpace(map) ? "(no map configured)" : map) +
                "\r\n\r\nWithout a map, Rcross places every marker on one chromosome. " +
                "Choose No, select the correct map (-m), and run again.\r\n\r\n" +
                "Continue with the one-chromosome fallback?",
                "Linkage map required", MessageBoxButtons.YesNo, MessageBoxIcon.Warning,
                MessageBoxDefaultButton.Button2);
            return result == DialogResult.Yes;
        }

        private void SynchronizeProjectState(CommandRequest completedRequest)
        {
            if (synchronizingProject || workingDirectoryBox == null)
                return;
            string directory = completedRequest == null ? workingDirectoryBox.Text.Trim() : completedRequest.WorkingDirectory;
            if (!Directory.Exists(directory))
                return;
            string resource = completedRequest == null ? resourceBox.Text.Trim() : completedRequest.ResourceFile;
            ProjectState state = ProjectState.Load(directory, resource);
            string stem = completedRequest == null ? "" : completedRequest.RequestedStem;
            string output;
            if (completedRequest != null && string.IsNullOrWhiteSpace(stem) &&
                completedRequest.Options != null && completedRequest.Options.TryGetValue("-o", out output))
                stem = ProjectState.StemFromFile(output);
            if (string.IsNullOrWhiteSpace(stem))
                stem = state.Stem;

            synchronizingProject = true;
            try
            {
                if (!string.IsNullOrWhiteSpace(stem))
                    stemBox.Text = stem;
                ApplyProjectDefaults(state);
            }
            finally
            {
                synchronizingProject = false;
            }
            UpdateCommandPreview();
        }

        private void ApplyProjectDefaults(ProjectState state)
        {
            if (selectedTool == null)
                return;
            Dictionary<string, string> mapping = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            mapping["-e"] = "-error";
            if (selectedTool.Name == "Rmap")
            {
                mapping["-i"] = "-mapin";
                mapping["-o"] = "-map";
            }
            else if (selectedTool.Name == "Rcross")
            {
                mapping["-i"] = "-iinfile";
                mapping["-o"] = "-ifile";
                mapping["-m"] = "-map";
                mapping["-q"] = "-qtl";
            }
            foreach (OptionEditor editor in editors)
            {
                string resourceKey;
                string value;
                if (mapping.TryGetValue(editor.Definition.Flag, out resourceKey) &&
                    state.Values.TryGetValue(resourceKey, out value) &&
                    !string.IsNullOrWhiteSpace(value))
                {
                    editor.Definition.DefaultValue = value;
                    if (!editor.Enabled.Checked)
                        editor.Value.Text = value;
                }
            }
        }

        private static string ResolveProjectPath(string directory, string path)
        {
            string clean = path.Trim().Trim('"');
            return Path.IsPathRooted(clean) ? clean : Path.Combine(directory, clean);
        }

        private void AppendOutput(string text, Color? color)
        {
            if (outputBox.InvokeRequired)
            {
                outputBox.BeginInvoke((MethodInvoker)delegate { AppendOutput(text, color); });
                return;
            }
            outputBox.SelectionStart = outputBox.TextLength;
            outputBox.SelectionColor = color.HasValue ? color.Value : outputBox.ForeColor;
            outputBox.AppendText(text);
            outputBox.SelectionColor = outputBox.ForeColor;
            outputBox.ScrollToCaret();
        }

        private void CancelCurrentProcess()
        {
            runningQueue = false;
            if (currentProcess != null && !currentProcess.HasExited)
            {
                currentProcess.Kill();
                statusLabel.Text = "Cancelling...";
            }
        }

        private void AddCurrentToQueue()
        {
            CommandRequest request = BuildCurrentRequest();
            if (request == null)
                return;
            request.Status = "Pending";
            queue.Add(request);
            RefreshQueue();
            statusLabel.Text = request.Tool.Name + " added to workflow queue";
        }

        private void RefreshQueue()
        {
            queueView.Items.Clear();
            for (int i = 0; i < queue.Count; i++)
            {
                CommandRequest item = queue[i];
                ListViewItem row = new ListViewItem((i + 1).ToString());
                row.SubItems.Add(item.Tool.Name);
                row.SubItems.Add(string.IsNullOrEmpty(item.Status) ? "Pending" : item.Status);
                row.SubItems.Add(item.Elapsed == TimeSpan.Zero ? "" : item.Elapsed.ToString(@"mm\:ss"));
                row.SubItems.Add(item.Arguments);
                row.SubItems.Add(item.WorkingDirectory);
                if (item.Status == "Succeeded") row.BackColor = Color.Honeydew;
                else if (item.Status == "Failed") row.BackColor = Color.MistyRose;
                else if (item.Status == "Running") row.BackColor = Color.LightCyan;
                queueView.Items.Add(row);
            }
        }

        private void RemoveSelectedQueueItems()
        {
            for (int i = queueView.SelectedIndices.Count - 1; i >= 0; i--)
                queue.RemoveAt(queueView.SelectedIndices[i]);
            RefreshQueue();
        }

        private void RunQueue()
        {
            if (queue.Count == 0)
            {
                MessageBox.Show(this, "The workflow queue is empty.", Text, MessageBoxButtons.OK, MessageBoxIcon.Information);
                return;
            }
            if (currentProcess != null && !currentProcess.HasExited)
                return;
            runningQueue = true;
            queueIndex = 0;
            foreach (CommandRequest item in queue)
            {
                item.Status = "Pending";
                item.Elapsed = TimeSpan.Zero;
                item.OverwriteConfirmed = false;
            }
            outputBox.Clear();
            workflowProgress.Minimum = 0;
            workflowProgress.Maximum = queue.Count;
            workflowProgress.Value = 0;
            RefreshQueue();
            RunNextQueueItem();
        }

        private void RunNextQueueItem()
        {
            if (queueIndex >= queue.Count)
            {
                runningQueue = false;
                statusLabel.Text = "Workflow completed";
                workflowProgress.Value = workflowProgress.Maximum;
                SaveProject(false);
                MessageBox.Show(this, "All workflow steps completed successfully.", Text, MessageBoxButtons.OK, MessageBoxIcon.Information);
                return;
            }
            StartRequest(queue[queueIndex], true);
        }

        private void RetryFailedStage()
        {
            int failed = queue.FindIndex(delegate(CommandRequest item) { return item.Status == "Failed"; });
            if (failed < 0)
            {
                MessageBox.Show(this, "There is no failed workflow stage to retry.", Text, MessageBoxButtons.OK, MessageBoxIcon.Information);
                return;
            }
            for (int i = failed; i < queue.Count; i++)
            {
                queue[i].Status = "Pending";
                queue[i].Elapsed = TimeSpan.Zero;
                queue[i].FailureDetails = "";
                queue[i].OverwriteConfirmed = false;
            }
            queueIndex = failed;
            runningQueue = true;
            RefreshQueue();
            RunNextQueueItem();
        }

        private void UpdateWorkflowProgress()
        {
            if (workflowProgress == null) return;
            int completed = queue.Count(item => item.Status == "Succeeded");
            workflowProgress.Maximum = Math.Max(1, queue.Count);
            workflowProgress.Value = Math.Min(workflowProgress.Maximum, completed);
        }

        private void UpdateElapsedTime()
        {
            if (currentProcess == null || currentProcess.HasExited) return;
            CommandRequest running = currentRequest;
            TimeSpan elapsed = running == null ? TimeSpan.Zero : DateTime.Now - running.StartedAt;
            elapsedLabel.Text = "Elapsed " + elapsed.ToString(@"hh\:mm\:ss");
            if (running != null)
            {
                running.Elapsed = elapsed;
                RefreshQueue();
            }
        }

        private ToolDefinition FindTool(string name)
        {
            return tools.Find(delegate(ToolDefinition tool) { return tool.Name == name; });
        }

        private void LoadSampleWorkflow()
        {
            string work = workingDirectoryBox.Text.Trim();
            if (work.Length == 0)
                return;
            Directory.CreateDirectory(work);
            string example = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "example");
            string map = Path.Combine(example, "sample.mps");
            string cross = Path.Combine(example, "sample.raw");
            if (!File.Exists(map) || !File.Exists(cross))
            {
                MessageBox.Show(this, "The packaged sample data could not be found.", Text, MessageBoxButtons.OK, MessageBoxIcon.Error);
                return;
            }
            File.Copy(map, Path.Combine(work, "sample.mps"), true);
            File.Copy(cross, Path.Combine(work, "sample.raw"), true);
            DialogResult explain = MessageBox.Show(this,
                "Load example analysis copies sample.mps and sample.raw into the current working directory and replaces the workflow queue. Continue?",
                "Load example analysis", MessageBoxButtons.OKCancel, MessageBoxIcon.Information);
            if (explain != DialogResult.OK) return;
            queue.Clear();
            AddSampleStep("Rmap", "-i sample.mps -A", work);
            AddSampleStep("Rcross", "-i sample.raw -A", work);
            AddSampleStep("Qstats", "-A", work);
            AddSampleStep("LRmapqtl", "-A", work);
            AddSampleStep("SRmapqtl", "-A", work);
            AddSampleStep("Zmapqtl", "-A", work);
            AddSampleStep("MImapqtl", "-A", work);
            AddSampleStep("Eqtl", "-A", work);
            AddSampleStep("Preplot", "-A", work);
            RefreshQueue();
            statusLabel.Text = "Example analysis loaded";
        }

        private void AddSampleStep(string name, string arguments, string work)
        {
            queue.Add(new CommandRequest { Tool = FindTool(name), Arguments = arguments, WorkingDirectory = work, Status = "Pending" });
        }

        private void RefreshFiles()
        {
            filesView.Items.Clear();
            string directory = workingDirectoryBox.Text.Trim();
            if (!Directory.Exists(directory))
                return;
            foreach (FileInfo file in new DirectoryInfo(directory).GetFiles())
            {
                string category;
                string purpose = DescribeFile(file.Extension, out category);
                if (fileFilter != null && fileFilter.SelectedIndex > 0 &&
                    !string.Equals(fileFilter.SelectedItem.ToString(), category, StringComparison.OrdinalIgnoreCase))
                    continue;
                ListViewItem item = new ListViewItem(file.Name) { Tag = file.FullName };
                item.SubItems.Add(purpose);
                item.SubItems.Add(FormatSize(file.Length));
                item.SubItems.Add(file.LastWriteTime.ToString("g"));
                filesView.Items.Add(item);
            }
        }

        private static string DescribeFile(string extension, out string category)
        {
            switch (extension.ToLowerInvariant())
            {
                case ".mps": case ".raw": case ".inp": case ".map": case ".cro":
                    category = "Inputs"; return extension == ".map" ? "Linkage map" : extension == ".cro" ? "Cross/genotype data" : "Imported source data";
                case ".z": category = "Results"; return "Interval/composite interval mapping results";
                case ".eqt": category = "Results"; return "Estimated QTL peaks and effects";
                case ".lr": category = "Results"; return "Linear-regression mapping results";
                case ".sr": category = "Results"; return "Stepwise-regression mapping results";
                case ".qst": case ".mim": case ".mr": case ".bys": category = "Results"; return "Statistical analysis results";
                case ".log": case ".rc": category = "Logs"; return extension == ".rc" ? "Project settings/resource file" : "Analysis log";
                case ".plt": case ".png": case ".svg": category = "Plots"; return "Plot or chart";
                default: category = "All files"; return "Project file";
            }
        }

        private static string FormatSize(long bytes)
        {
            if (bytes >= 1024 * 1024) return (bytes / (1024.0 * 1024.0)).ToString("0.0") + " MB";
            if (bytes >= 1024) return (bytes / 1024.0).ToString("0.0") + " KB";
            return bytes + " B";
        }

        private void OpenSelectedFile()
        {
            if (filesView.SelectedItems.Count == 0)
                return;
            Process.Start((string)filesView.SelectedItems[0].Tag);
        }

        private void PreviewSelectedFile()
        {
            if (filesView.SelectedItems.Count == 0 || filePreview == null)
                return;
            string path = (string)filesView.SelectedItems[0].Tag;
            string extension = Path.GetExtension(path).ToLowerInvariant();
            if (new[] { ".png", ".jpg", ".jpeg", ".gif", ".exe", ".zip" }.Contains(extension))
            {
                filePreview.Text = "Binary file. Use Open selected to view it in the associated application.";
                return;
            }
            try
            {
                using (StreamReader reader = new StreamReader(path, Encoding.Default, true))
                {
                    char[] buffer = new char[200000];
                    int read = reader.Read(buffer, 0, buffer.Length);
                    filePreview.Text = new string(buffer, 0, read) + (reader.Peek() >= 0 ? "\r\n\r\n[Preview truncated]" : "");
                    filePreview.SelectionStart = 0;
                    filePreview.ScrollToCaret();
                }
            }
            catch (Exception ex) { filePreview.Text = "Preview unavailable: " + ex.Message; }
        }

        private void OpenWorkingDirectory()
        {
            string directory = workingDirectoryBox.Text.Trim();
            if (directory.Length == 0)
                return;
            Directory.CreateDirectory(directory);
            Process.Start("explorer.exe", Quote(directory));
        }

        private void NewProject()
        {
            using (NewProjectWizardForm wizard = new NewProjectWizardForm(workingDirectoryBox.Text.Trim()))
            {
                if (wizard.ShowDialog(this) != DialogResult.OK) return;
                Directory.CreateDirectory(wizard.ProjectDirectory);
                workingDirectoryBox.Text = wizard.ProjectDirectory;
                stemBox.Text = wizard.ProjectStem;
                queue.Clear();
                string mapInput;
                string crossInput;
                if (wizard.UseExample)
                {
                    string example = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "example");
                    mapInput = Path.Combine(wizard.ProjectDirectory, "sample.mps");
                    crossInput = Path.Combine(wizard.ProjectDirectory, "sample.raw");
                    File.Copy(Path.Combine(example, "sample.mps"), mapInput, true);
                    File.Copy(Path.Combine(example, "sample.raw"), crossInput, true);
                }
                else
                {
                    mapInput = Path.Combine(wizard.ProjectDirectory, Path.GetFileName(wizard.MapFile));
                    crossInput = Path.Combine(wizard.ProjectDirectory, Path.GetFileName(wizard.CrossFile));
                    if (!string.Equals(wizard.MapFile, mapInput, StringComparison.OrdinalIgnoreCase)) File.Copy(wizard.MapFile, mapInput, true);
                    if (!string.Equals(wizard.CrossFile, crossInput, StringComparison.OrdinalIgnoreCase)) File.Copy(wizard.CrossFile, crossInput, true);
                }
                AddSampleStep("Rmap", "-X " + Quote(wizard.ProjectStem) + " -i " + Quote(Path.GetFileName(mapInput)) + " -A", wizard.ProjectDirectory);
                AddSampleStep("Rcross", "-X " + Quote(wizard.ProjectStem) + " -i " + Quote(Path.GetFileName(crossInput)) + " -A", wizard.ProjectDirectory);
                AddSampleStep("Qstats", "-X " + Quote(wizard.ProjectStem) + " -A", wizard.ProjectDirectory);
                AddSampleStep("LRmapqtl", "-X " + Quote(wizard.ProjectStem) + " -A", wizard.ProjectDirectory);
                AddSampleStep("SRmapqtl", "-X " + Quote(wizard.ProjectStem) + " -A", wizard.ProjectDirectory);
                AddSampleStep("Zmapqtl", "-X " + Quote(wizard.ProjectStem) + " -A", wizard.ProjectDirectory);
                AddSampleStep("MImapqtl", "-X " + Quote(wizard.ProjectStem) + " -A", wizard.ProjectDirectory);
                AddSampleStep("Eqtl", "-X " + Quote(wizard.ProjectStem) + " -A", wizard.ProjectDirectory);
                AddSampleStep("Preplot", "-X " + Quote(wizard.ProjectStem) + " -A", wizard.ProjectDirectory);
                currentProjectFile = Path.Combine(wizard.ProjectDirectory, wizard.ProjectStem + ".qtlproject");
                RefreshQueue();
                RefreshFiles();
                SaveProject(false);
                tabs.SelectedIndex = 2;
                statusLabel.Text = "New project created and validated";
            }
        }

        private void SaveProject(bool saveAs)
        {
            if (saveAs || string.IsNullOrEmpty(currentProjectFile))
            {
                using (SaveFileDialog dialog = new SaveFileDialog
                {
                    Filter = "QTL Cartographer project (*.qtlproject)|*.qtlproject",
                    InitialDirectory = Directory.Exists(workingDirectoryBox.Text) ? workingDirectoryBox.Text : "",
                    FileName = (string.IsNullOrEmpty(stemBox.Text) ? "qtlcart" : stemBox.Text) + ".qtlproject"
                })
                {
                    if (dialog.ShowDialog(this) != DialogResult.OK) return;
                    currentProjectFile = dialog.FileName;
                }
            }
            ProjectDocument document = new ProjectDocument
            {
                WorkingDirectory = workingDirectoryBox.Text.Trim(),
                Stem = stemBox.Text.Trim(),
                ResourceFile = resourceBox.Text.Trim(),
                SelectedTool = selectedTool == null ? "" : selectedTool.Name,
                SelectedArguments = BuildCurrentRequest() == null ? "" : BuildCurrentRequest().Arguments
            };
            document.InputHashes = ProjectStore.HashInputs(document.WorkingDirectory);
            foreach (CommandRequest request in queue)
                document.Queue.Add(new ProjectStep
                {
                    Tool = request.Tool.Name, Arguments = request.Arguments,
                    Status = request.Status, ElapsedSeconds = request.Elapsed.TotalSeconds
                });
            if (Directory.Exists(document.WorkingDirectory))
                document.Results.AddRange(Directory.GetFiles(document.WorkingDirectory)
                    .Where(path => new[] { ".z", ".eqt", ".lr", ".sr", ".qst", ".mim", ".plt" }.Contains(Path.GetExtension(path).ToLowerInvariant()))
                    .Select(Path.GetFileName));
            ProjectStore.Save(currentProjectFile, document);
            statusLabel.Text = "Project saved: " + Path.GetFileName(currentProjectFile);
        }

        private void OpenProject()
        {
            using (OpenFileDialog dialog = new OpenFileDialog { Filter = "QTL Cartographer project (*.qtlproject)|*.qtlproject" })
            {
                if (dialog.ShowDialog(this) != DialogResult.OK) return;
                ProjectDocument document;
                try { document = ProjectStore.Load(dialog.FileName); }
                catch (Exception ex)
                {
                    MessageBox.Show(this, "The project could not be opened:\r\n" + ex.Message, Text, MessageBoxButtons.OK, MessageBoxIcon.Error);
                    return;
                }
                currentProjectFile = dialog.FileName;
                List<string> changedInputs = ProjectStore.DetectChangedInputs(document);
                workingDirectoryBox.Text = document.WorkingDirectory;
                stemBox.Text = document.Stem;
                resourceBox.Text = document.ResourceFile;
                queue.Clear();
                foreach (ProjectStep step in document.Queue)
                    queue.Add(new CommandRequest
                    {
                        Tool = FindTool(step.Tool), Arguments = step.Arguments,
                        WorkingDirectory = document.WorkingDirectory, Status = step.Status,
                        Elapsed = TimeSpan.FromSeconds(step.ElapsedSeconds)
                    });
                RefreshQueue();
                RefreshFiles();
                statusLabel.Text = "Project opened (created with version " + document.Version + ")";
                if (changedInputs.Count > 0)
                    MessageBox.Show(this,
                        "These inputs changed after the project was saved:\r\n" + string.Join("\r\n", changedInputs.ToArray()) +
                        "\r\n\r\nExisting results may no longer be reproducible. Rerun the affected workflow stages.",
                        "Input data changed", MessageBoxButtons.OK, MessageBoxIcon.Warning);
            }
        }

        private void ShowResultsDashboard()
        {
            string directory = workingDirectoryBox.Text.Trim();
            if (!Directory.Exists(directory))
            {
                MessageBox.Show(this, "Open or create a project first.", Text, MessageBoxButtons.OK, MessageBoxIcon.Information);
                return;
            }
            using (ResultsDashboardForm dashboard = new ResultsDashboardForm(directory, stemBox.Text.Trim()))
                dashboard.ShowDialog(this);
        }

        private CrossData LoadCrossData()
        {
            string directory = workingDirectoryBox.Text.Trim();
            string stem = stemBox.Text.Trim();
            string map = Path.Combine(directory, (stem.Length == 0 ? "qtlcart" : stem) + ".map");
            string cross = Path.Combine(directory, (stem.Length == 0 ? "qtlcart" : stem) + ".cro");
            return CrossDataParser.Load(cross, map);
        }

        private void ShowDiagnostics()
        {
            CrossData data = LoadCrossData();
            if (data.Individuals.Count == 0)
            {
                MessageBox.Show(this, "No native cross data is available. Run Rmap and Rcross first.", Text, MessageBoxButtons.OK, MessageBoxIcon.Information);
                return;
            }
            using (DiagnosticsForm diagnostics = new DiagnosticsForm(data)) diagnostics.ShowDialog(this);
        }

        private void ShowModernFormats()
        {
            using (ModernImportForm form = new ModernImportForm(workingDirectoryBox.Text.Trim(), stemBox.Text.Trim(), LoadCrossData()))
                form.ShowDialog(this);
            RefreshFiles();
        }

        private void ShowBenchmarks()
        {
            using (BenchmarkForm form = new BenchmarkForm(workingDirectoryBox.Text.Trim(), stemBox.Text.Trim())) form.ShowDialog(this);
        }

        private void ShowAnalysisDesign()
        {
            string directory = workingDirectoryBox.Text.Trim();
            if (!Directory.Exists(directory)) return;
            using (AnalysisDesignForm form = new AnalysisDesignForm(directory, stemBox.Text.Trim()))
                if (form.ShowDialog(this) == DialogResult.OK)
                {
                    string stemArgument = stemBox.Text.Trim().Length == 0 ? "" : "-X " + Quote(stemBox.Text.Trim()) + " ";
                    AddSampleStep(form.RecommendedTool, stemArgument + "-A", directory);
                    RefreshQueue();
                    statusLabel.Text = "Validated " + form.RecommendedTool + " design added to workflow";
                }
        }

        private void ShowTutorial()
        {
            using (TutorialForm tutorial = new TutorialForm(delegate { LoadSampleWorkflow(); })) tutorial.ShowDialog(this);
        }

        private void ExportReproducibilityBundle()
        {
            string directory = workingDirectoryBox.Text.Trim();
            if (!Directory.Exists(directory)) return;
            using (SaveFileDialog dialog = new SaveFileDialog
            {
                Filter = "Reproducibility ZIP (*.zip)|*.zip",
                FileName = (stemBox.Text.Trim().Length == 0 ? "qtlcart" : stemBox.Text.Trim()) + "-reproducibility.zip"
            })
                if (dialog.ShowDialog(this) == DialogResult.OK)
                {
                    List<string> included = ReproducibilityBundle.GetIncludedFiles(directory, stemBox.Text.Trim(), currentProjectFile, dialog.FileName);
                    string[] preview = included.Take(25).Select(Path.GetFileName).ToArray();
                    string remaining = included.Count > preview.Length ? "\r\n...and " + (included.Count - preview.Length) + " more files." : "";
                    DialogResult confirmation = MessageBox.Show(this,
                        "The bundle will contain " + included.Count + " project files:\r\n\r\n" +
                        string.Join("\r\n", preview) + remaining +
                        "\r\n\r\nGenotype and phenotype files may contain sensitive data. Review this list before sharing. Continue?",
                        "Review reproducibility bundle", MessageBoxButtons.OKCancel, MessageBoxIcon.Warning);
                    if (confirmation != DialogResult.OK) return;
                    ReproducibilityBundle.Create(directory, stemBox.Text.Trim(), currentProjectFile, queue, dialog.FileName);
                    statusLabel.Text = "Reproducibility bundle created";
                }
        }

        private void SaveConsoleLog()
        {
            using (SaveFileDialog dialog = new SaveFileDialog())
            {
                dialog.Filter = "Text log (*.txt)|*.txt|All files (*.*)|*.*";
                dialog.FileName = "qtl-cartographer-console.txt";
                if (dialog.ShowDialog(this) == DialogResult.OK)
                    File.WriteAllText(dialog.FileName, outputBox.Text, Encoding.UTF8);
            }
        }

        private void ShowProgramHelp()
        {
            Form dialog = new Form
            {
                Text = selectedTool == null ? "Program help" : selectedTool.Name + " help",
                Size = new Size(850, 650),
                StartPosition = FormStartPosition.CenterParent,
                Icon = Icon
            };
            RichTextBox text = new RichTextBox
            {
                Dock = DockStyle.Fill,
                ReadOnly = true,
                Font = new Font("Consolas", 9F),
                Text = lastHelpText,
                BackColor = Color.White,
                WordWrap = false
            };
            dialog.Controls.Add(text);
            dialog.ShowDialog(this);
        }
    }
}
