# Gene Expression Analysis Project

# Introduction
Gene expression analysis is an extremely important step in many biological and clinical studies since it 
allows scientists to find similarities and differences in the transcriptome of cells, as well as to some extent,
in their genome and proteome. Those findings open up opportunities to study a variety of cellular 
processes and diseases that may arise from their disorders. Therefore, this project aims to
analyze RNA expression data in the cloned subline of neuroblastoma cells and interpret the results.

# Methods
The Gene Expression analysis was performed in the IDE RStudio using the following libraries:
1) To work with the imported data: dplyr, tidyr, plyr, tibble, reshape2
2) To perform the gene expression analysis itself: airway, DESeq2, tximeta
3) To annotate the results of gene expression analysis: AnnotationDbi, org.Hs.eg.db
4) To perform pathway analysis: clusterProfiler, RCy3, DOSE, ReactomePA. 
5) To visualize the results: ggplot2, ggbeeswarm, pheatmap, genefilter

The gene expression analysis was based on the results of the experiment featuring the treatment of SHSY5Y/6TR(EU)/pTrex-Dest-30/MYCN (SY5Y-MYCN) cell line with Doxycycline in the concentration of 1μg/ml 
in order to induce overexpression of MYCN. Total RNA of cells was extracted at 0 hours, 1 hour, 4 hours 
and 24 hours after the Doxycycline injection using TRI Reagent according to the manufacturer’s protocol and 
sequenced on an Illumina GA IIx following the Illumina protocols. Gene expression quantification was 
performed using the Salmon tool and GRCh38 reference genome. The experiment was performed in two 
iterations, the results of which were used for this analysis.

First of all, the metadata file was read and quant.sf files for a gene-expression time-course in response to
MYCN induction of samples was loaded. Then, the transcript quantification was performed via tximeta
library, and those results were summarized at the gene level. Next, the DESeq data structure was set up in 
order to be used for the following differential gene analysis. Genes with low expression counts (<10 on the 
sum of all measurements of the two experiments) were dropped from the analysis, performed using the 
DESeq () function. The results of expression analysis were then filtered by padj and log2FoldChange variables 
to keep only the genes that have a statistically significant change of expression (p adjusted value <0.05) with a big enough 
effect size (log fold change value >1 or <-1). Then, annotation of genes was performed, and the short symbol "name" was added to 
the end of each table row. Finally, the three tables comparing differential gene expression (1 Hour vs 
Control, 4 Hours vs Control and 24 Hours vs Control) were saved in .csv format.

Next, quality control of the given samples was performed to check if the experimental impact on cells 
was successful and to see how the expression pattern of genes changed through the course of 
measurements in the two experimental iterations. To do this, a graph for counts of MYCN gene was created 
using the plotCounts() and ggplot() functions by looking for the Ensemble ID of ENSG00000134323. Then the data for visualization was normalized via variance shrinkage algorithm, and the top 20 most variable 
genes were selected to create a Principal Components Plot via the plotPCA() function of the BiocGenerics library and 
clustered Heatmap via the pheatmap() function of pheatmap library.

Afterwards, the three created tables of performed differential gene expression analysis were loaded and 
seven genes, present in all three time points, were selected to draw a graph illustrating how their expression was changing 
using the ggplot() function. 

Then, the pathway analysis was performed via the Gene Ontology search[1], Gene Set Enrichment 
Analysis[2], Kyoto Encyclopedia of Genes and Genomes Analysis[3] and Reactome Pathway analysis[4]. The Gene 
Ontology search was performed on all three tables to find out in what process the differentially expressed 
genes are involved. First, all the genes that were not differentially expressed were dropped, and only the genes that had 
a p-value<0.05 were selected. Then, the bitr() function was used to find the ENTREZ ID for each present ENSEMBL ID and inputted 
the result into the enrichGO() function. Afterwards, the ggplot() function was used to draw bar plots for Gene 
Ontology.

Gene Set Enrichment Analysis was the hardest part of Pathway analysis as the GSEA() function used to 
perform it required an ordered named list of Entrez ID and their p-values, which were in two different tables 
of different sizes with duplicated observations. So, first of all, the created during the Gene Ontology search 
ENSEMBL ID – ENTREZ ID data frame was used to drop all the genes from the p-values table that do not have their 
own Entrez ID, and then for each gene, its Entrez ID was obtained via the already mentioned bitr() function. If the 
function returned more than one ID for a gene, the Ensemble IDs and p-values were duplicated to 
maintain the structure: 1 Ensemble ID – 1 p-value – 1 Entrez ID. The result was saved in the dds_enterzid.csv 
file. Then the Homo sapiens Hallmarks collection of genes was loaded, and all p-values were converted to p-scores 
using the formula: p-score = -10*log(p-value). Finally, the list for the GSEA() function was created, and its results 
were visualized using the dotplot() function.

Next, was the Kyoto Encyclopedia of Genes and Genomes Analysis performed via the enrichKEGG() function, 
using the generated Entrez ID data frame and visualized via the dotplot() function. In the end, the Reactome 
analysis was performed using the enrichPathway() function. The returned names of pathways were then 
visualized using the viewPathway() function.

# Results

The performed quality control has highlighted the fact that some genes in the same experiment time were 
expressing differently, as can be seen, looking at the clustered Heatmap below for gene 234487 when 
comparing two controls, gene 280518 when comparing cells at 1 Hour, gene 274442 when comparing cells 
at 24 Hours, et cetera.

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/d95bbab9-e942-4ac6-aa1f-9fa40250ebf1)
