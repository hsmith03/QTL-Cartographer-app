using System.Collections.Generic;

namespace QTLCartographer.Gui
{
    internal static class ToolCatalog
    {
        public static List<ToolDefinition> Create()
        {
            return new List<ToolDefinition>
            {
                new ToolDefinition("Rmap", "Data preparation", "Create, simulate, or translate a genetic linkage map."),
                new ToolDefinition("Rqtl", "Data preparation", "Create, simulate, or translate a genetic model."),
                new ToolDefinition("Rcross", "Data preparation", "Create, simulate, or translate a mapping population."),
                new ToolDefinition("Prune", "Data preparation", "Prune data sets and prepare bootstrap samples."),
                new ToolDefinition("Qstats", "Exploration", "Calculate marker and quantitative-trait summary statistics."),
                new ToolDefinition("Emap", "Exploration", "Estimate a linkage map from marker data."),
                new ToolDefinition("LRmapqtl", "QTL analysis", "Map QTLs with simple linear regression."),
                new ToolDefinition("SRmapqtl", "QTL analysis", "Map QTLs with stepwise regression."),
                new ToolDefinition("Zmapqtl", "QTL analysis", "Perform interval and composite interval mapping."),
                new ToolDefinition("JZmapqtl", "QTL analysis", "Perform joint or multitrait composite interval mapping."),
                new ToolDefinition("MImapqtl", "QTL analysis", "Build and refine multiple-interval mapping models."),
                new ToolDefinition("BTmapqtl", "QTL analysis", "Map QTLs for binary traits."),
                new ToolDefinition("Bmapqtl", "QTL analysis", "Perform Bayesian and reversible-jump QTL mapping."),
                new ToolDefinition("MultiRegress", "QTL analysis", "Perform multiple-regression analysis."),
                new ToolDefinition("Eqtl", "Results", "Estimate QTL positions and effects from mapping results."),
                new ToolDefinition("Preplot", "Results", "Generate plotting data and gnuplot command files.")
            };
        }
    }
}
