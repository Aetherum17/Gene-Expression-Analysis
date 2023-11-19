# Gene Expression Analysis Project

# Introduction
Gene expression analysis is an extremely important step in many biological and clinical studies since it 
allows scientists to find similarities and differences in the transcriptome of cells, as well as, to some extent,
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

The gene expression analysis was based on the results of the experiment featuring the treatment of SH-SY5Y/6TR(EU)/pTrex-Dest-30/MYCN (SY5Y-MYCN) cell line with Doxycycline in the concentration of 1μg/ml 
in order to induce overexpression of MYCN. The total RNA of cells was extracted at 0 hours, 1 hour, 4 hours 
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
using the plotCounts() and ggplot() functions by looking for the Ensemble ID of ENSG00000134323. Then, the data for visualization was normalized via variance shrinkage algorithm, and the top 20 most variable 
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

The results of the PCA also confirm the difference in gene expression, as only the cells at 
1 and 4 hours are grouped nearly together. This can indicate that the experiment was not repeated in the same 
conditions.

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/e43b3f2a-ca5e-4782-bc79-485eac72d372)

However, the differential expression results were obtained and are included in 
resultsTOT_1h_TOT_0h.csv, resultsTOT_4h_TOT_0h.csv, resultsTOT_24h_TOT_0h.csv files. As they contain 
a lot of genes, the top 5 for each comparison is presented.

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/c91f1515-73c1-4e68-b538-9cd632ff7f4a)

As can be seen, in all comparisons, the top three genes are the same: MYCN, SMIM11 and LOC102723553. 
The presence of the MYCN gene is not surprising, as the cells were induced by Doxycycline in the 
experiment. The success of the intervention can also be confirmed by looking at the graph below, which shows 
the increase in the number of counts of the MYCN gene in the samples through the course of observations.

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/8ab69215-10b1-4332-810b-e14ba7ae40c3)

The other two genes that are present in the top 3 list of differential expression encode Small Integral 
Membrane proteins. The latter genes in top-5 of all tables are VARS1, PPP1R10, GUSBP9, TIGD5 and GIMAP2. 
They are Valine-TRNA Ligase, Protein Phosphatase 1 Regulatory Subunit 10, Tigger Transposable Element 
and GTP-binding/immuno-associated nucleotide-binding protein respectively.

Among the 3 comparisons, not all the genes were present in each table. However, there are 7 of 
them, differentially expressed at 1, 4 and 24 Hours time points. The graph showing how their expression changed compared to the control group after the treatment with Doxycycline is presented below. As can be seen from 
this graph, PPP1R10 and LOC102723553 were firstly expressed at bigger rates, but after 4th hours 
they declined, SMIM11 and MYCN expression after an increase remained the same, while expression of 
RIMBP3C was constantly increasing, and finally, the expression of CTSO and VARS1 does not seem to change 
after induction by the drug.

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/15697ec8-585c-4a02-bc15-c638b6f4a34a)

Another interesting finding is the decrease in expression of the TP53TG3B gene at the 24-hour time point, 
that may play a significant role in the p53/TP53-mediating signalling pathway, which is crucial for preventing the formation of cancer.

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/a412d779-26fd-4387-ba5e-8a4e1497aaed)

After obtaining the results of differential gene expression analysis, the Gene Ontology 
analysis was performed to look at the biological processes in which the differentially expressed genes participate. The 
analysis was performed with the enrichGO() function of the clusterProfiler library.
As can be seen from the bar plot below, after one hour, most genes are involved in protein synthesis 
(ribosomes) and transport (Golgi apparatus). 

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/5caa3c59-6862-4db0-8d89-03b4c0de80aa)

After 4 hours, only the change of activity of the Golgi apparatus can be seen.

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/e28cf687-23cb-4625-8f9d-c060cc4ed055)

While after 24 hours, we can see an increase in the activity inside the nucleus (Cajal bodies, telomerases), 
the synthesis of nucleotide phosphates, that could be used for DNA replication / RNA production or as an 
energy source. 

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/f94d0caf-3450-420e-99a9-5fee978b37ac)

The next step was performing the Gene Set Enrichment Analysis with the GSEA() function of the 
clusterProfiler library to associate the hallmarks of diseases with certain genes. As the results of the Gene 
Ontology search has not shown any signs of abnormal cellular activity at 1 and 4 Hours time points, GSEA 
analysis was performed using the differential gene expression data of 24 Hours vs Control, where an increase 
in processes related to Cajal bodies was observed. The results of Gene Set Enrichment Analysis are 
summarized in the dot plot below, showing the number of genes associated with human Hallmarks by gene 
ratio - the number of genes related to the Gene Ontology process divided by the total number of expected 
genes. Counts of such genes and p-adjusted values are also shown.

As can be seen on the graph below, at 24 hours after exposure of cells to Doxycycline, hallmarks of normal 
cell processes such as oxidative phosphorylation or fatty acid metabolism are present along with the 
hallmarks of MYC and E2F targets, associated with the genesis of cancer. Also, a hallmark of the G2M
checkpoint is present, which can indicate the proximity of the Mitosis stage, though in some tumours it is 
also connected with mutations leading to malignant transformation.

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/04d68813-4fc0-4bca-87b4-667757cd6948)

Afterwards, the analysis using the Kyoto Encyclopedia of Genes and Genomes was performed via enrichKEGG() 
function on the Entrez ID of genes, differentially expressed at 24 Hour time points, since GSEA has shown 
the presence of cancer hallmarks there. As can be seen on the dot plot below, KEGG analysis has found 
groups of genes, associated with various diseases, among which we can also see the process connected with 
tumorogenesis: Chemical Carcinogenesis – reactive oxygen species, Thyroid cancer and ErbB signalling 
pathway

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/caeda737-581c-4cc0-8220-ffa14393ba1f)

Finally, the Reactome analysis was performed to see if it could find any pathways connected to 
cancer in which the differentially expressed genes on 24 Hour time point are involved. However, the results 
of the enrichPathway() function have not found any, as can be seen from the tables below.

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/ac9246fc-ad9f-4600-a226-15958f285618)

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/e28de73c-ffca-4631-8499-b13ad892e5c9)
![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/a65c02df-91d9-40a7-8fc7-3c3037468028)








