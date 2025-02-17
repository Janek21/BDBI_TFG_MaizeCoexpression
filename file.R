# Load necessary libraries
install.packages("WGCNA")  # Install if not already installed
library(WGCNA)

# Set options for WGCNA
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# 1. Load Gene Expression Data ----------------------------------------
# Replace with your actual gene expression data file (rows = genes, cols = samples)
# Ensure data is properly formatted (genes in rows, samples in columns)
expression_data <- read.csv("gene_expression.csv", row.names = 1)

# Transpose data: WGCNA requires genes as columns and samples as rows
datExpr <- as.data.frame(t(expression_data))

# Check for missing values
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# 2. Choose Soft Threshold -------------------------------------------
powers <- c(1:20)  # Power range for scale-free topology
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Select optimal power (based on scale-free topology criterion)
softPower <- sft$powerEstimate
if (is.na(softPower)) softPower <- 6  # Default to 6 if no estimate

# 3. Construct Adjacency & TOM --------------------------------------
adjacency <- adjacency(datExpr, power = softPower)  # Build adjacency matrix
TOM <- TOMsimilarity(adjacency)  # Transform into Topological Overlap Matrix
dissTOM <- 1 - TOM  # Convert to dissimilarity

# 4. Identify Modules using Hierarchical Clustering -----------------
geneTree <- hclust(as.dist(dissTOM), method = "average")  # Cluster genes

# Dynamic tree cut to detect modules
minModuleSize <- 30  # Minimum module size
dynamicModules <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize)

# Convert labels to colors
moduleColors <- labels2colors(dynamicModules)

# 5. Relate Modules to Eigengenes -----------------------------------
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs)

# 6. Visualize Network & Modules ------------------------------------
# Plot dendrogram with module colors
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, moduleColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Save module assignments
moduleAssignments <- data.frame(Gene = colnames(datExpr), Module = moduleColors)
write.csv(moduleAssignments, "module_assignments.csv", row.names = FALSE)

# 7. Export Network for Cytoscape (Optional) ------------------------
cyt <- exportNetworkToCytoscape(TOM, 
                                edgeFile = "Cytoscape_edges.txt", 
                                nodeFile = "Cytoscape_nodes.txt",
                                weighted = TRUE, 
                                threshold = 0.1, 
                                nodeAttr = moduleColors)

# 8. Save Results ---------------------------------------------------
save(datExpr, adjacency, TOM, moduleColors, MEs, file = "WGCNA_results.RData")

print("WGCNA Network Construction Completed!")

