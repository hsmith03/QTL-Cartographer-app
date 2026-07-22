using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Text;
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
        private ToolStripStatusLabel statusLabel;
        private ToolDefinition selectedTool;
        private Process currentProcess;
        private bool runningQueue;
        private int queueIndex;
        private string lastHelpText = "";

        private string ToolsDirectory
        {
            get { return Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "tools"); }
        }

        public MainForm()
        {
            Text = "QTL Cartographer";
            Icon = SystemIcons.Application;
            MinimumSize = new Size(1050, 700);
            Size = new Size(1280, 820);
            StartPosition = FormStartPosition.CenterScreen;
            BackColor = Pale;
            Font = new Font("Segoe UI", 9F);

            BuildInterface();
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
        }

        private void BuildInterface()
        {
            Panel banner = new Panel { Dock = DockStyle.Top, Height = 76, BackColor = Navy };
            Label appName = new Label
            {
                Text = "QTL Cartographer",
                ForeColor = Color.White,
                Font = new Font("Segoe UI Semibold", 21F),
                AutoSize = true,
                Location = new Point(24, 12)
            };
            Label subtitle = new Label
            {
                Text = "Quantitative trait locus analysis workspace",
                ForeColor = Color.FromArgb(190, 211, 227),
                AutoSize = true,
                Location = new Point(27, 50)
            };
            banner.Controls.Add(appName);
            banner.Controls.Add(subtitle);
            Controls.Add(banner);

            StatusStrip status = new StatusStrip { SizingGrip = false };
            statusLabel = new ToolStripStatusLabel("Ready") { Spring = true, TextAlign = ContentAlignment.MiddleLeft };
            status.Items.Add(statusLabel);
            Controls.Add(status);

            SplitContainer split = new SplitContainer
            {
                Dock = DockStyle.Fill,
                SplitterDistance = 255,
                FixedPanel = FixedPanel.Panel1,
                BackColor = Border
            };
            Controls.Add(split);
            split.BringToFront();
            banner.BringToFront();

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

            TabControl tabs = new TabControl { Dock = DockStyle.Fill };
            tabs.TabPages.Add(BuildConfigurePage());
            tabs.TabPages.Add(BuildOutputPage());
            tabs.TabPages.Add(BuildQueuePage());
            tabs.TabPages.Add(BuildFilesPage());
            tabs.TabPages.Add(BuildAboutPage());
            content.Controls.Add(heading, 0, 0);
            content.Controls.Add(tabs, 0, 1);
            split.Panel2.Controls.Add(content);
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
            workingDirectoryBox.TextChanged += update;
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
            Button sample = MakeButton("Load sample workflow", false);
            sample.Click += delegate { LoadSampleWorkflow(); };
            Button remove = MakeButton("Remove selected", false);
            remove.Click += delegate { RemoveSelectedQueueItems(); };
            Button clear = MakeButton("Clear", false);
            clear.Click += delegate { queue.Clear(); RefreshQueue(); };
            toolbar.Controls.Add(run);
            toolbar.Controls.Add(sample);
            toolbar.Controls.Add(remove);
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
            queueView.Columns.Add("Arguments", 510);
            queueView.Columns.Add("Working directory", 260);
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
            toolbar.Controls.Add(refresh);
            toolbar.Controls.Add(open);
            toolbar.Controls.Add(folder);

            filesView = new ListView
            {
                Dock = DockStyle.Fill,
                View = View.Details,
                FullRowSelect = true,
                GridLines = true
            };
            filesView.Columns.Add("Name", 350);
            filesView.Columns.Add("Size", 100);
            filesView.Columns.Add("Modified", 170);
            filesView.DoubleClick += delegate { OpenSelectedFile(); };
            page.Controls.Add(filesView);
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
                Text = "QTL Cartographer for Windows\r\n\r\n" +
                       "A native Windows desktop interface for QTL Cartographer 1.17. " +
                       "All statistical calculations are performed by the original GPL-licensed C engine.\r\n\r\n" +
                       "The GUI exposes every program option, supports queued workflows, streams console output, " +
                       "and keeps generated files organized in a project directory.\r\n\r\n" +
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
                    TextBox value = new TextBox
                    {
                        Text = option.DefaultValue,
                        Dock = DockStyle.Fill,
                        Enabled = false,
                        Margin = new Padding(3, 5, 3, 5)
                    };
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
                    value.TextChanged += delegate { UpdateCommandPreview(); };
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
                UpdateCommandPreview();
            }
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
            foreach (OptionEditor editor in editors)
            {
                if (!editor.Enabled.Checked)
                    continue;
                arguments.Add(editor.Definition.Flag);
                if (editor.Value.Text.Trim().Length > 0)
                    arguments.Add(Quote(editor.Value.Text.Trim()));
            }
            if (stemBox.Text.Trim().Length > 0)
            {
                arguments.Add("-X");
                arguments.Add(Quote(stemBox.Text.Trim()));
            }
            if (resourceBox.Text.Trim().Length > 0)
            {
                arguments.Add("-R");
                arguments.Add(Quote(resourceBox.Text.Trim()));
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
                WorkingDirectory = workingDirectoryBox.Text.Trim()
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
            string executable = Path.Combine(ToolsDirectory, request.Tool.Name + ".exe");
            if (!File.Exists(executable))
            {
                MessageBox.Show(this, "Missing analysis program:\r\n" + executable, Text, MessageBoxButtons.OK, MessageBoxIcon.Error);
                return;
            }

            AppendOutput("\r\n> " + request.DisplayCommand + "\r\n", Color.FromArgb(100, 190, 235));
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
                if (e.Data != null) AppendOutput(e.Data + Environment.NewLine, Color.FromArgb(255, 150, 140));
            };
            process.Exited += delegate
            {
                process.WaitForExit();
                int code = process.ExitCode;
                BeginInvoke((MethodInvoker)delegate
                {
                    AppendOutput("[Exited with code " + code + "]\r\n", code == 0 ? Color.FromArgb(125, 210, 145) : Color.FromArgb(255, 150, 140));
                    currentProcess = null;
                    runButton.Enabled = selectedTool != null;
                    cancelButton.Enabled = false;
                    statusLabel.Text = code == 0 ? "Completed " + request.Tool.Name : request.Tool.Name + " failed (exit " + code + ")";
                    RefreshFiles();
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
                            MessageBox.Show(this, "The workflow stopped because " + request.Tool.Name + " failed.", Text, MessageBoxButtons.OK, MessageBoxIcon.Error);
                        }
                    }
                });
            };

            try
            {
                currentProcess = process;
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
                runButton.Enabled = true;
                cancelButton.Enabled = false;
                MessageBox.Show(this, ex.Message, Text, MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
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
                row.SubItems.Add(item.Arguments);
                row.SubItems.Add(item.WorkingDirectory);
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
            outputBox.Clear();
            RunNextQueueItem();
        }

        private void RunNextQueueItem()
        {
            if (queueIndex >= queue.Count)
            {
                runningQueue = false;
                statusLabel.Text = "Workflow completed";
                MessageBox.Show(this, "All workflow steps completed successfully.", Text, MessageBoxButtons.OK, MessageBoxIcon.Information);
                return;
            }
            StartRequest(queue[queueIndex], true);
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
            statusLabel.Text = "Sample workflow loaded";
        }

        private void AddSampleStep(string name, string arguments, string work)
        {
            queue.Add(new CommandRequest { Tool = FindTool(name), Arguments = arguments, WorkingDirectory = work });
        }

        private void RefreshFiles()
        {
            filesView.Items.Clear();
            string directory = workingDirectoryBox.Text.Trim();
            if (!Directory.Exists(directory))
                return;
            foreach (FileInfo file in new DirectoryInfo(directory).GetFiles())
            {
                ListViewItem item = new ListViewItem(file.Name) { Tag = file.FullName };
                item.SubItems.Add(FormatSize(file.Length));
                item.SubItems.Add(file.LastWriteTime.ToString("g"));
                filesView.Items.Add(item);
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

        private void OpenWorkingDirectory()
        {
            string directory = workingDirectoryBox.Text.Trim();
            if (directory.Length == 0)
                return;
            Directory.CreateDirectory(directory);
            Process.Start("explorer.exe", Quote(directory));
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
