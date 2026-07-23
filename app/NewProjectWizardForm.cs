using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Windows.Forms;

namespace QTLCartographer.Gui
{
    internal sealed class ImportInspection
    {
        public int Chromosomes { get; set; }
        public int Markers { get; set; }
        public int Individuals { get; set; }
        public int Traits { get; set; }
        public int MissingValues { get; set; }
        public string CrossType { get; set; }
        public List<string> Messages { get; private set; }
        public bool IsValid { get; set; }

        public ImportInspection() { CrossType = ""; Messages = new List<string>(); }
    }

    internal static class DataImportValidator
    {
        public static ImportInspection Inspect(string mapFile, string crossFile, string selectedCross)
        {
            ImportInspection result = new ImportInspection { IsValid = true, CrossType = selectedCross };
            if (!File.Exists(mapFile)) { result.Messages.Add("Map file does not exist."); result.IsValid = false; }
            if (!File.Exists(crossFile)) { result.Messages.Add("Genotype/phenotype file does not exist."); result.IsValid = false; }
            if (!result.IsValid) return result;
            string map = File.ReadAllText(mapFile);
            Match chrom = Regex.Match(map, @"(?im)(?:-chromosomes|\*Chromosomes:)\s*(\d+)");
            if (chrom.Success) result.Chromosomes = int.Parse(chrom.Groups[1].Value);
            else result.Chromosomes = Regex.Matches(map, @"(?im)^\s*(?:-Chromosome|\*c\d+)").Count;
            result.Markers = Regex.Matches(map, @"(?im)^\s*(?:-marker|\*[A-Za-z0-9])").Count;

            string cross = File.ReadAllText(crossFile);
            Match type = Regex.Match(cross, @"(?im)(?:-Cross\s+|data\s+type\s+)([^\s]+(?:\s+intercross)?)");
            if (type.Success) result.CrossType = type.Groups[1].Value.Trim();
            Match sample = Regex.Match(cross, @"(?im)(?:-SampleSize\s+|^\s*)(\d+)\s+\d+\s+\d+\s*$");
            if (sample.Success) result.Individuals = int.Parse(sample.Groups[1].Value);
            Match traits = Regex.Match(cross, @"(?im)-traits\s+(\d+)");
            if (traits.Success) result.Traits = int.Parse(traits.Groups[1].Value);
            else
            {
                Match rawHeader = Regex.Match(cross, @"(?im)^\s*\d+\s+\d+\s+(\d+)\s*$");
                if (rawHeader.Success) result.Traits = int.Parse(rawHeader.Groups[1].Value);
            }
            result.MissingValues = Regex.Matches(cross, @"(?:--|(?<=\s)-(?=\s)|(?<=\s)\.(?=\s))").Count;
            if (result.Chromosomes <= 0) { result.Messages.Add("No chromosome assignments were detected in the map."); result.IsValid = false; }
            if (result.Markers <= 0) result.Messages.Add("Marker count could not be determined; the native importer will perform the final format check.");
            if (result.Individuals <= 0) result.Messages.Add("Sample size could not be determined; verify the cross header.");
            if (result.Traits <= 0) { result.Messages.Add("No phenotype traits were detected."); result.IsValid = false; }
            if (!string.IsNullOrEmpty(selectedCross) && !string.IsNullOrEmpty(result.CrossType) &&
                result.CrossType.IndexOf(selectedCross, StringComparison.OrdinalIgnoreCase) < 0)
                result.Messages.Add("Selected cross type differs from the type declared by the input file.");
            if (result.Messages.Count == 0) result.Messages.Add("Map, chromosome assignments, genotype/phenotype input, and cross metadata passed validation.");
            return result;
        }
    }

    internal sealed class NewProjectWizardForm : Form
    {
        private RadioButton exampleChoice;
        private TextBox folderBox;
        private TextBox stemBox;
        private TextBox mapBox;
        private TextBox crossBox;
        private ComboBox crossType;
        private TextBox validation;
        private Button finish;

        public bool UseExample { get { return exampleChoice.Checked; } }
        public string ProjectDirectory { get { return folderBox.Text.Trim(); } }
        public string ProjectStem { get { return stemBox.Text.Trim(); } }
        public string MapFile { get { return mapBox.Text.Trim(); } }
        public string CrossFile { get { return crossBox.Text.Trim(); } }

        public NewProjectWizardForm(string initialDirectory)
        {
            Text = "New project";
            Size = new Size(720, 590);
            StartPosition = FormStartPosition.CenterParent;
            AutoScaleMode = AutoScaleMode.Dpi;
            FormBorderStyle = FormBorderStyle.FixedDialog;
            MaximizeBox = false;
            TableLayoutPanel grid = new TableLayoutPanel { Dock = DockStyle.Fill, Padding = new Padding(18), ColumnCount = 3, RowCount = 9 };
            grid.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 145));
            grid.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 100));
            grid.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 90));
            Label intro = new Label { Text = "Choose example data for a guided demonstration or provide your own map and genotype/phenotype files.", Dock = DockStyle.Fill, AutoSize = true };
            grid.Controls.Add(intro, 0, 0); grid.SetColumnSpan(intro, 3);
            FlowLayoutPanel choices = new FlowLayoutPanel { Dock = DockStyle.Fill };
            exampleChoice = new RadioButton { Text = "Example data", Checked = true, AutoSize = true };
            RadioButton mine = new RadioButton { Text = "My data", AutoSize = true };
            choices.Controls.Add(exampleChoice); choices.Controls.Add(mine);
            grid.Controls.Add(new Label { Text = "Project type", TextAlign = ContentAlignment.MiddleLeft, Dock = DockStyle.Fill }, 0, 1);
            grid.Controls.Add(choices, 1, 1); grid.SetColumnSpan(choices, 2);
            folderBox = AddPath(grid, 2, "Project folder", initialDirectory, true);
            stemBox = AddText(grid, 3, "Filename stem", "qtlcart");
            mapBox = AddPath(grid, 4, "Linkage map input", "", false);
            crossBox = AddPath(grid, 5, "Genotype + phenotype", "", false);
            crossType = new ComboBox { Dock = DockStyle.Fill, DropDownStyle = ComboBoxStyle.DropDownList };
            crossType.Items.AddRange(new object[] { "Auto-detect", "B1", "B2", "F2", "RI", "SF", "RF" });
            crossType.SelectedIndex = 0;
            grid.Controls.Add(new Label { Text = "Cross type", Dock = DockStyle.Fill, TextAlign = ContentAlignment.MiddleLeft }, 0, 6);
            grid.Controls.Add(crossType, 1, 6); grid.SetColumnSpan(crossType, 2);
            validation = new TextBox { Dock = DockStyle.Fill, Multiline = true, ReadOnly = true, ScrollBars = ScrollBars.Vertical };
            grid.Controls.Add(validation, 0, 7); grid.SetColumnSpan(validation, 3);
            FlowLayoutPanel actions = new FlowLayoutPanel { Dock = DockStyle.Fill, FlowDirection = FlowDirection.RightToLeft };
            finish = new Button { Text = "Create project", AutoSize = true };
            Button validate = new Button { Text = "Validate inputs", AutoSize = true };
            Button cancel = new Button { Text = "Cancel", AutoSize = true, DialogResult = DialogResult.Cancel };
            finish.Click += Finish;
            validate.Click += delegate { ValidateInputs(); };
            actions.Controls.Add(finish); actions.Controls.Add(validate); actions.Controls.Add(cancel);
            grid.Controls.Add(actions, 0, 8); grid.SetColumnSpan(actions, 3);
            Controls.Add(grid);
            exampleChoice.CheckedChanged += delegate { SetMode(); };
            SetMode();
            AcceptButton = finish;
            CancelButton = cancel;
        }

        private void SetMode()
        {
            mapBox.Enabled = crossBox.Enabled = !UseExample;
            crossType.Enabled = !UseExample;
            validation.Text = UseExample
                ? "Example analysis copies sample.mps and sample.raw into the project folder, then prepares the complete analysis workflow."
                : "Select both input files and validate them before creating the project.";
        }

        private bool ValidateInputs()
        {
            if (UseExample) return true;
            ImportInspection inspection = DataImportValidator.Inspect(MapFile, CrossFile, crossType.Text == "Auto-detect" ? "" : crossType.Text);
            validation.Text = "Chromosomes: " + inspection.Chromosomes + "\r\nMarkers detected: " + inspection.Markers +
                "\r\nIndividuals: " + inspection.Individuals + "\r\nCross: " + inspection.CrossType +
                "\r\nPhenotype traits: " + inspection.Traits +
                "\r\nMissing-value tokens: " + inspection.MissingValues + "\r\n\r\n" + string.Join("\r\n", inspection.Messages.ToArray());
            return inspection.IsValid;
        }

        private void Finish(object sender, EventArgs e)
        {
            if (ProjectDirectory.Length == 0 || ProjectStem.Length == 0)
            {
                MessageBox.Show(this, "Choose a project folder and filename stem.", Text, MessageBoxButtons.OK, MessageBoxIcon.Warning);
                return;
            }
            if (!ValidateInputs()) return;
            DialogResult = DialogResult.OK;
            Close();
        }

        private TextBox AddPath(TableLayoutPanel grid, int row, string label, string value, bool folder)
        {
            TextBox box = AddText(grid, row, label, value);
            Button browse = new Button { Text = "Browse…", Dock = DockStyle.Fill };
            browse.Click += delegate
            {
                if (folder)
                {
                    using (FolderBrowserDialog dialog = new FolderBrowserDialog())
                        if (dialog.ShowDialog(this) == DialogResult.OK) box.Text = dialog.SelectedPath;
                }
                else
                {
                    using (OpenFileDialog dialog = new OpenFileDialog())
                        if (dialog.ShowDialog(this) == DialogResult.OK) box.Text = dialog.FileName;
                }
            };
            grid.Controls.Add(browse, 2, row);
            return box;
        }

        private TextBox AddText(TableLayoutPanel grid, int row, string label, string value)
        {
            TextBox box = new TextBox { Text = value, Dock = DockStyle.Fill };
            grid.Controls.Add(new Label { Text = label, Dock = DockStyle.Fill, TextAlign = ContentAlignment.MiddleLeft }, 0, row);
            grid.Controls.Add(box, 1, row);
            grid.SetColumnSpan(box, 2);
            return box;
        }
    }
}
