# Loading required Libraries
library("airway")
library("dplyr")
library("tidyr")
library("tximeta")
library("DESeq2")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("ggbeeswarm")
library("ggplot2")
library("genefilter")
library("pheatmap")
library("plyr")
library("tibble")
library("reshape2")

# Get the path of directory of the project and the sample files located there
dir <- getwd()
sample_files <- list.files(path = paste(dir, "/Samples", sep=""))

# Read the Metadata file, rename the columns them and add a column with the path of samples
metadata_csv <- read.table(file.path(dir, "Samples/samplelist.txt"),sep=' ')
colnames(metadata_csv) <- c("time", "sample")
metadata_csv$files <- file.path(dir, "Samples", metadata_csv$sample, "quant.sf")
# save result
write.csv(metadata_csv,paste(dir, "/Samples/metadata.csv", sep=""), row.names = FALSE)

# Check if the quant.sf files of samples exist
file.exists(metadata_csv$files)

# tximeta() requires a data frame of two columns - file paths and names of experiments
# First we set the names in correct order
codata_names <- factor(metadata_csv[,1], levels = c("TOT_0h", "TOT_1h", "TOT_4h", "TOT_24h"))
coldata <- data.frame(files = metadata_csv[,3], names = codata_names, sample_names = as.factor(metadata_csv[,2]))
# locate and download the relevant annotation data
### Row names are transcripts
se <- tximeta(coldata)
# summarize the transcript-level quantifications to the gene level
### Row Names are genes
gse <- summarizeToGene(se)

# Setting up the DESeq data structure. We know that control is TOT_0h, and the rest are
# experiments, so in design we can just input ~names.
dds <- DESeqDataSet(gse, design = ~names)

# Filtering - keep only genes that are expressed
keep <- rowSums(assay(dds)) >= 10
dds <- dds[keep,]

# Performing the gene differential expression analysis
dds <- DESeq(dds)

# List for result file names
result_files <- list()

# For loop to compare all observations with control
for(i in 2:length(levels(dds$names))){
  res <- results(dds, contrast=c("names",levels(dds$names)[i], levels(dds$names)[1]))
  
  # Filter the results by p-value and by log2 fold change
  # Remove genes with NA for p-value, as, anyway, counts for this gene are zero
  res.filtered <- res[complete.cases(res[ , 5:6]),]
  res.filtered <- res.filtered[(res.filtered$pvalue<0.05),]
  res.filtered <- res.filtered[(res.filtered$log2FoldChange > 1) | (res.filtered$log2FoldChange < -1), ]
  
  summary(res.filtered)
  
  # Annotation of genes
  ens.str <- rownames(res.filtered)
  res.filtered$symbol <- mapIds(org.Hs.eg.db,
                                keys=ens.str,
                                column="SYMBOL",
                                keytype="ENSEMBL",
                                multiVals="first")
  
  # Exporting the results - most significant are on the top
  res.filtered.ordered <- res.filtered[order(res.filtered$pvalue),]
  head(res.filtered.ordered)
  write.csv(res.filtered.ordered, file = paste(dir,"/results",levels(dds$names)[i],"_",levels(dds$names)[1],".csv", sep=""))
  # Log the names of files
  result_files <- append(result_files, paste("results",levels(dds$names)[i],"_",levels(dds$names)[1],".csv", sep=""))
}

# graph for MYCN gene
geneCounts <- plotCounts(dds, gene = "ENSG00000134323", intgroup = c("names"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = names, y = count, color = names)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

# Normalization of the data for visualization via variance shrinkage algorithm.
vsd <- vst(dds, blind = TRUE)
# Select the top 20 most variable genes for visualization
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)

### QUALITY CONTROL

## clustering analysis with heatmap to see how the expression top 20 genes changes.
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
# Rename the columns/rows/levels for better look
rownames(colData(vsd)) <- c("24 Hours", "0 Hours", "1 Hour", "4 Hours", "24 Hours", "0 Hours", "1 Hour", "4 Hours")
colnames(colData(vsd)) <- c("Time", "Sample")
levels(colData(vsd)$Time) <- c("0 Hours", "1 Hour", "4 Hours", "24 Hours")

# Create heatmap
anno <- as.data.frame(colData(vsd)[, c("Time", "Sample")])
pheatmap(mat, annotation_col = anno, annotation_legend = T, annotation_names_col= T)

## Create Principal Components Plot
colData(vsd) <- colData(vsd)[1]
plotPCA(vsd, intgroup = c("Time"))

### Create a graph for expression changes of the genes, present in all of 1/4/24 Hours estimations
# Create a summary data frame
log2_results <- data.frame()
# Loop for iterating through result .csv files
for (p in 1:length(result_files)){
  res_csv <- read.csv(file.path(dir, result_files[[p]]),sep=',')
  res_csv <- subset(res_csv, select=c(X, log2FoldChange, symbol))
  # We are intrested only in Gene ID, its short name and effect size, rest can be dropped
  colnames(res_csv) <- c("Gene", paste("log2FoldChange", p, sep="_"), "symbol")
  # Merge the tables
  if(ncol(log2_results)==0){
    log2_results <- res_csv
  }
  log2_results <- merge(log2_results, res_csv, by=c("Gene"))
}
# Remove duplicated columns
log2_results <- log2_results[,c(2,6,8,9)]
# Transpose the data frame for creating the graph and fix the consequneces of this action
log2_results <- t(log2_results)
colnames(log2_results) <- log2_results[4,]
rownames(log2_results) <- c("1", "4", "24", "x")
log2_results <- log2_results[-4,]

# Convert transposed data to data frame
log2_results <- data.frame(log2_results)
str(log2_results)

# Convert all char values back to floats
for(q in 1:ncol(log2_results)){
  log2_results[,q] <- as.numeric(log2_results[,q])
  print(str(log2_results[,q]))
}
# Add index of hours for graph
log2_results$index <- c(1, 4, 24)

#Create the graph
ggplot(data=log2_results, aes(x=index)) +
  geom_line(aes(y = MYCN), color = "red", show.legend =T)+geom_point(aes(y = MYCN), color = "red")+
  geom_line(aes(y = RIMBP3C), color = "green")+geom_point(aes(y = RIMBP3C), color = "green")+
  geom_line(aes(y = SMIM11), color = "blue")+geom_point(aes(y = SMIM11), color = "blue")+
  geom_line(aes(y = PPP1R10), color = "magenta1")+geom_point(aes(y = PPP1R10), color = "magenta1")+
  geom_line(aes(y = VARS1), color = "salmon")+geom_point(aes(y = VARS1), color = "salmon")+
  geom_line(aes(y = CTSO), color = "yellow")+geom_point(aes(y = CTSO), color = "yellow")+
  geom_line(aes(y = LOC102723553), color = "pink")+geom_point(aes(y = LOC102723553), color = "pink")+
  labs(x="Hours", y="Gene Expression Chanage (log2)")+labs(title="Gene Expression Change during the Experiments")+
  theme()

