
#Required packages


# Install if needed:
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("limma", "edgeR", "clusterProfiler", "Rsubread"))

library(limma)           # voom, lmFit, eBayes, topTable, topTreat, plotMDS
library(edgeR)           # DGEList-related RNA-seq objects
library(clusterProfiler) # enrichKEGG
library(Rsubread)        # read.columns (if annotation steps are later added)


#PART 1:  Summary plot- mapped vs unmapped reads

barplot(
  as.matrix(
    t(
      cbind(
        propmap[, "NumMapped"] * 1e-6,
        (propmap[, "NumTotal"] - propmap[, "NumMapped"]) * 1e-6
      )
    )
  ),
  ylab = "Number of Reads (millions)",
  names = Targets[1:6, 1],
  las = 3,
  cex.names = 0.3,
  cex.axis = 0.7,
  col = c("violet", "mediumpurple"),
  legend = c("Number of Mapped Reads", "Number of Unmapped Reads"),
  args.legend = list(x = "topright", cex = 0.5, bty = "n")
)


#PART 2: Summary plot-library size

barplot(
  colSums(Counts) * 1e-6,
  names = Targets[1:6, 1],
  ylab = "Library size (millions)",
  las = 2,
  cex.names = 0.3,
  col = "violet"
)

legend(
  "topleft",
  legend = c("1.5", "2.5", "3.5"),
  pch = 15,
  cex = 0.7,
  col = "purple4"
)


#PART 3: Multidimensional scaling (MDS) plot

plotMDS(
  y,
  labels = Targets[1:6, 1],
  gene.selection = "common",
  prior.count = 5
)



#PART 4: Differential expression analysis setup

design <- model.matrix(~ 0 + factor(c(rep(1, 3), rep(2, 3))))
colnames(design) <- c("APR246", "DMSO")

# Define the contrast of interest
contrast.matrix <- makeContrasts(APR246 - DMSO, levels = design)
contrast.matrix

#PART 5: limma-voom workflow

# Apply voom transformation
v <- voom(y, design)

# Fit linear model
fit <- lmFit(v, design)

# Apply contrast
fit2 <- contrasts.fit(fit, contrast.matrix)

# Empirical Bayes moderation
fit3 <- eBayes(fit2)

# View top results
topTable(fit3)

# Check dimensions of filtered top table with fold-change and p-value cutoffs
dim(topTable(fit3, number = 1000, lfc = 1.2, p.value = 0.05))


#PART 6: Extract significant genes and top 5 up/down regulated

# Convert results to a data frame
degdf <- as.data.frame(topTreat(fit3, number = Inf))

# Filter for significant DEGs using adjusted p-value
sig_deg <- subset(degdf, adj.P.Val < 0.05)

# Fallback in case no genes pass threshold
if (nrow(sig_deg) == 0) sig_deg <- degdf

# Top 5 upregulated genes (highest logFC)
top5_up <- head(sig_deg[order(-sig_deg$logFC), ], 5)

# Top 5 downregulated genes (lowest logFC)
top5_down <- head(sig_deg[order(sig_deg$logFC), ], 5)

# Display selected columns
top5_up[, c("GeneID", "Symbol", "logFC", "adj.P.Val")]
top5_down[, c("GeneID", "Symbol", "logFC", "adj.P.Val")]


#PART 7: Volcano plot

# Recreate DEG results table
degdf <- as.data.frame(topTreat(fit3, number = Inf))

# Thresholds used in the report
lfc_thr <- 1.2
fdr_thr <- 0.05

# Extract all genes for volcano plot
DEG <- topTable(fit3, coef = "APR246 - DMSO", number = Inf)
head(DEG)

# Base volcano plot
plot(
  DEG$logFC,
  -log10(DEG$adj.P.Val),
  pch = 16,
  cex = 0.6,
  col = "purple4",
  xlab = "log2 Fold Change",
  ylab = "-log10(adjusted p-value)",
  main = "Volcano Plot (APR246 vs DMSO)"
)

# Add threshold lines
abline(v = c(-lfc_thr, lfc_thr), lty = 2, col = "black")
abline(h = -log10(fdr_thr), lty = 2, col = "black")

# Highlight top 5 upregulated genes
points(
  top5_up$logFC,
  -log10(top5_up$adj.P.Val),
  col = "violet",
  pch = 16,
  cex = 1.2
)

# Highlight top 5 downregulated genes
points(
  top5_down$logFC,
  -log10(top5_down$adj.P.Val),
  col = "mediumpurple",
  pch = 16,
  cex = 1.2
)

# Label top 5 upregulated genes
text(
  top5_up$logFC,
  -log10(top5_up$adj.P.Val),
  labels = top5_up$Symbol,
  pos = 3,
  cex = 0.7,
  col = "violet"
)

# Label top 5 downregulated genes
text(
  top5_down$logFC,
  -log10(top5_down$adj.P.Val),
  labels = top5_down$Symbol,
  pos = 3,
  cex = 0.7,
  col = "mediumpurple"
)



#PART 8: KEGG pathway enrichment analysis

de <- topTreat(fit3, number = 1000, p.value = 0.05)[, 1]

yy <- enrichKEGG(
  de,
  organism = "hsa",
  pvalueCutoff = 0.5
)

# View enrichment results
head(yy)

# Barplot of enriched KEGG pathways
barplot(yy)

# Convert enrichment output to data frame
yydf <- as.data.frame(yy)

# Extract top rows as shown in the report
sig_yy <- yydf[0:3, ]

# Keep selected columns
sig_yy_df <- sig_yy[, c("category", "ID", "Description", "p.adjust", "Count")]

sig_yy_df
