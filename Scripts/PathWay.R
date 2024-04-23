#####################################################################################

library("clusterProfiler")
library("org.Hs.eg.db")
library("DOSE")
library("DESeq2")
library("RCy3")
library("ggplot2")
library("ReactomePA")

# Loading the results of differential gene expression analysis
# Pathway analysis is performed by comparing the results of 24 Hours with control
path_dds <- results(dds, contrast=c("names","TOT_24h", "TOT_0h"))
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

###########################################################
# Because preparing the data for Gene Set Enrichment Analysis takes a lot of time
# The results of the preparation were written in a standalone file, read in the next segment.
# Thus this code should be run ONLY if dds_enterzid.csv file is not present.
# Prepare data for Pathway Analysis
gsea_genes <- path_dds
# Add Ensemble ID as own column and convert the data to data frame format
gsea_genes$ENSEMBL <- rownames(gsea_genes)
gsea_genes <- data.frame(gsea_genes)

# Create Dataframe: ENSEMBLE - ENTREZ
gsea_genes_entrez <- bitr(gsea_genes$ENSEMBL,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
# Keep only genes, that have ENTREZ ID
gsea_genes_entrez <- gsea_genes_entrez[complete.cases(gsea_genes_entrez[ , 2]),]

# In the data for Pathway analysis we drop all the genes, that do not have ENTREZ ID
genes_to_keep <- unname(unlist(gsea_genes_entrez[,1]))
gsea_genes <- gsea_genes[gsea_genes$ENSEMBL %in% genes_to_keep,]

# Create data frame for GSEA
gsea_genes_df <- gsea_genes[1,]

# This loop is for merging two data frames:
# 1) Has ENSEMBL ID and p-values
# 2) Has ENSEMBL ID and ENTREZ ID
# The loop was created because ENSEMBL ID can have several ENTREZ ID and not every ENSEMBL ID has own ENTREZ ID,
# which made it not possible to simply merge the data frames.
# The list for keeping all ENTREZ ID, that ENSEMBL ID has
zlist <- c()
# Basically, we recreate the gsea_genes data frame, interating through each gene
# Takes 30 minutes to run
for(gsea_gene in 1:nrow(gsea_genes))
{
  # Find ENTREZ ID individually
  zid <- bitr(gsea_genes[gsea_gene,7],fromType = "ENSEMBL",
              toType = "ENTREZID",OrgDb = org.Hs.eg.db)[,2]
  # Duplicate the ENSEMBL ID entry to maitain the structure: 1 Ensemble ID - 1 p-value - 1 Entrz ID
  for (z in 1:length(zid)){
    gsea_genes_df[nrow(gsea_genes_df)+1,] <- unname(unlist(gsea_genes[gsea_gene,]))
    zlist <- append(zlist, zid[z])
  }
}
# Drop the first line created to define columns
gsea_genes_df <- gsea_genes_df[-1,]
# Now we add ENTREZ ID as a new column
gsea_genes_df$ENTREZID <- zlist

# And finally save this monstrosity, as I am not going to re run it again
write.csv(gsea_genes_df, file = paste(dir,"/dds_enterzid.csv", sep=""))

#########################################################################

# Read the saved data frame for Gene Set Enrichment Analysis
gsea_genes_csv <- read.csv(paste(dir,"/dds_enterzid.csv", sep=""))
# Select only the columns: p value and ENTREZ ID
gsea_genes_csv <- gsea_genes_csv[,c(6, 9)]
# Conversion of p value to scores
gsea_genes_csv$pvalue <- -10*log(gsea_genes_csv$pvalue)

# Create a named ranked genelist for GSEA() function
genelist <- setNames(gsea_genes_csv$pvalue, gsea_genes_csv$ENTREZID)

# Get all human gene sets and subset the Hallmark collection
m_df <- msigdbr(species = "Homo sapiens")
m_df_H <- m_df[m_df$gs_cat=="H",]
# Create a data frame of two columns: gene set name and entrez id
msigdbr_t2g = m_df_H %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()

# Finally, perform the GSE Analysis
res_GSEA <- GSEA(genelist, TERM2GENE = msigdbr_t2g)

# Convert gseaResult to data frame
res_GSEA_df <- as.data.frame(res_GSEA)

# Create dot plot
dotplot(res_GSEA, showCategory=22)

################################################################################

## KEGG analysis using the gseKEGG function
genes_kk <- unname(genes.entrez[,2])
kk <- enrichKEGG(gene = genes_kk, organism = 'hsa', pvalueCutoff = 0.05)
head(kk)
kk[6,]
dotplot(kk, showCategory=100)

## Reactome analysis:

rt <- enrichPathway(gene=genes_kk, pvalueCutoff = 0.05, readable=TRUE)
head(rt)
rt[c(1:3),c(2,3:4, 6)]
rt[,c(2)]

genes_no_duplication <- genelist[!duplicated(names(genelist))]
head(genes_no_duplication)

# visualise a pathway
viewPathway("Cellular response to heat stress", 
            readable = TRUE, 
            foldChange = genes_no_duplication)

