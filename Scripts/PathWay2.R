#####################################################################################

library("clusterProfiler")
library("org.Hs.eg.db")
library("DOSE")
library("DESeq2")
library("RCy3")
library("ggplot2")
library("ReactomePA")

# Loading the results of differential gene expression analysis
# Pathway analysis is performed by comparing the results of 1 Hours with control
path_dds <- results(dds, contrast=c("names","TOT_1h", "TOT_0h"))
# Drop the genes that are not expressed and order them by p value
path_dds <- path_dds[complete.cases(path_dds[ , 5:6]),]
path_dds <- path_dds[order(path_dds$pvalue),]

# Subset the differentialy expressed genes
genes <- rownames(path_dds[path_dds$pvalue < 0.05,])
bkgd.genes <- rownames(path_dds)

# Define Entrez ID 
# Used bitr package for convertion ensemble ID into Entrez ID
genes.entrez <- bitr(genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez <- bitr(bkgd.genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

# Perform Gene Ontology search
egobp <- clusterProfiler::enrichGO(
  gene     = genes.entrez[[2]],
  universe = bkgd.genes.entrez[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "ALL",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  readable = TRUE)

# Draw Barplot for Gene Ontology
ggplot(egobp[1:25], aes(x=reorder(Description, -pvalue), y=Count, fill=-p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low="blue", high="red") +
  labs(x = "", y = "Count", fill = "p.adjust") +
  theme(axis.text=element_text(size=11))