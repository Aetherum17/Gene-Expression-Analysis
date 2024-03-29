# Gene Expression Analysis Project

# Introduction
Gene expression analysis is an extremely important step in many biological and clinical studies since it 
allows scientists to find similarities and differences in the transcriptome of cells, as well as, to some extent,
in their genome and proteome. Those findings open up opportunities to study a variety of cellular 
processes and diseases that may arise from their disorders. Therefore, this project aims to
analyze RNA expression data in the cloned subline of neuroblastoma cells and interpret the results.

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
and GTP-binding/immuno-associated nucleotide-binding protein, respectively.

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
the synthesis of nucleotide phosphates that could be used for DNA replication / RNA production or as an 
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
hallmarks of MYC and E2F targets associated with the genesis of cancer. Also, a hallmark of the G2M
checkpoint is present, which can indicate the proximity of the Mitosis stage, though in some tumours, it is 
also connected with mutations leading to malignant transformation.

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/04d68813-4fc0-4bca-87b4-667757cd6948)

Afterwards, the analysis using the Kyoto Encyclopedia of Genes and Genomes was performed via enrichKEGG() 
function on the Entrez ID of genes, differentially expressed at 24-hour time points, since GSEA has shown 
the presence of cancer hallmarks there. As can be seen on the dot plot below, KEGG analysis has found 
groups of genes associated with various diseases, among which we can also see the process connected with 
tumorogenesis: Chemical Carcinogenesis – reactive oxygen species, Thyroid cancer and ErbB signalling 
pathway

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/caeda737-581c-4cc0-8220-ffa14393ba1f)

Finally, the Reactome analysis was performed to see if it could find any pathways connected to 
cancer in which the differentially expressed genes on 24-hour time points are involved. However, the results 
of the enrichPathway() function have not found any, as can be seen from the tables below.

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/ac9246fc-ad9f-4600-a226-15958f285618)

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/e28de73c-ffca-4631-8499-b13ad892e5c9)
![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/a65c02df-91d9-40a7-8fc7-3c3037468028)

# Discussion

From the results section, it can be seen that overall, the usage of Doxycycline has caused the desired effect, 
as in Figure 3, we see an increase in counts of the MYCN gene. This gene is a proto-oncogene, deregulation 
of which occurs in various types of cancer, including Neuroblastoma - a malignant tumour of the sympathetic 
nervous system[5,6]. However, since this disease looks to have a complex nature, it usually involves the 
mutation of other genes as well, such as p53, which plays an important role in the p53/TP53-mediating 
signalling pathway[5]. Indeed, in our case, we can also observe a decrease in the expression of the TP53TG3B 
gene, as mentioned in the results, which is also part of the p53/TP53 pathway7
, as can be seen in Figure 5.

Another possible sign indicating that at 24 24-hour time point, the cell starts forming a tumour can be seen in 
figure 8, where there is a high number of processes associated with Cajal Bodies – present in neuronal and 
cancer cells, but not in other normal diploid cells[8]. Though this also can be caused by the neuronal origin of 
the tumour, despite the fact that such activity was not observed in 1 and 4-hour time points (Figures 6 and 
7).

The later performed Gene Set Enrichment Analysis further highlights the cancerogenic nature of studied 
cells, as shows the presence of hallmarks MYC Targets V1, MYC Targets V2, G2 Checkpoint and E2F Targets 
(Figure 9). MYC is a group of C-MYC, L-MYC and N-MYC proteins, which are part of the MYC oncogenic family, 
consisting of MYCC, MYCL and already mentioned MYCN genes[9]. MYC proteins act as oncogenic 
transcriptional factors, as they promote cell growth and proliferation, alter apoptosis and activate 
telomerases. The latter could also be observed in the Gene Ontology results of Figure 8, where a number of 
processes related to Telomerase RNA, used as a template for telomere replication by telomerases, are 
present. G2M checkpoint is another detected cancer hallmark that controls the transition from G2 to M 
stages of the cell cycle and marks an increased aggressiveness of a tumour10. The same applies to the last 
unveiled cancer hallmark - E2F Targets – which mark cancer cells with increased expression of proliferation-related genes[11].

After confirming that the cells in the experiment had obtained multiple hallmarks of cancer, KEGG analysis 
was performed (Figure 10). It has found genes associated with various neurodegeneration diseases, such as 
Alzheimer, Parkinson or Huntington, as well as a couple of groups for malignant tumour development, but 
not for Neuroblastoma. Apparently, the https://www.genome.jp/kegg/pathway.html website just does not 
have any pathway map designed specifically for it, as can be seen in Figure 13, most likely due to mentioned 
earlier complex biology of the disease.

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/e0d457b1-e5b1-4c27-a681-6ff7dac43679)

Thus, a closer look was taken at the Chemical carcinogenesis - reactive oxygen species pathway. This 
is a complex pathway[12], as can be seen at Figure 15, so let’s take a closer look at the ENSG00000134184
and ENSG00000277897 genes, as they show the biggest effect size change among all differentially 
expressed genes related to Chemical carcinogenesis reactive oxygen species pathway, if judged by log2 fold 
change metrics. This pair of genes encode glutathione S-transferase enzyme[13,14], that in the Kegg database 
have the symbols of GSTT215 and GSTM116 and are shown in Figure 16 under the symbol GSTO1.

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/0e558bb4-ab66-4ed7-aa47-4c26a4e3b25c)

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/84887763-bfe0-4ad1-90c7-6ae5b3b1a6fc)

![image](https://github.com/Aetherum17/Gene-Expression-Analysis/assets/46795020/4c78fd59-df18-48e7-8fc6-2ef83a63ca31)

Glutathione transferases play an important role in the processes of cell detoxification, modification and 
synthesis of leukotrienes as well as prostaglandins. Those enzymes also help to avoid oxidative stress. 
Glutathione transferases are capable of catalyzing the conjuration of glutathione (GSH) to a variety of 
hydrophobic and electrophilic molecules, making them less toxic and prone to further modifications to later 
be discharged from cells[17, 18].

However, glutathione transferases can also participate in the signalling pathways, as they control the activity 
of mitogen-activated protein kinases. This feature is used by some tumour cells, for example, to bind to the 
c-Jun N-terminal kinase of the MAPK signalling pathway, deactivate the protein and therefore inhibit 
apoptotic signal[17, 19]. In the same way, glutathione transferases can interfere with the work of tumour 
necrosis factor receptor-associated factor 2 and p3819.

Moreover, glutathione transferases can use their conjugating function to provide anti-cancer drug resistance 
to malignant cells by detoxification of medicine. Indeed, in many cancers, an increased expression of this 
enzyme’s gene is shown compared to normal cells17. The same result is obtained in our observations (Figure 
14).

Thus, it is possible to assume that increased levels of Glutathione S-transferase will contribute to the 
survival, proliferation and drug resistance of Neuroblastoma cells, as this enzyme will be able to disrupt the 
work of multiple pathways, including MAPK, and make the cancer cells more chemoresistant. 

Finally, it is worth noting that based on the obtained results, no more additional experiments should be conducted with this cell line,
until more repeats for the current experiment are completed since the results of quality control imply that two 
iterations of the experiments were conducted under different conditions, as at the same time point in two 
repeats, we see a different expression of the same genes on the heatmap (Figure 1). Additionally, the principal 
component analysis has put only one pair of samples out of four close to each other (Figure 2). Moreover, 
the results of Reactome analysis show the activity of paths related to viruses: Export of Viral 
Ribonucleoproteins from the Nucleus, Rev-mediated nuclear export of HIV RNA, SARS-CoV-2 
activates/modulates innate and adaptive immune responses, SARS-CoV-2-host interactions, Viral Messenger 
RNA Synthesis, Host Interactions of HIV factors, HIV Infection (Figure 12). Hence, maybe it is also worth
testing the cell culture on viral contamination, as its presence can nullify the results of all future conducted 
work with this material.

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

# References:
1. Carlson M (2019). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.8.2.
2. https://www.gsea-msigdb.org/gsea/index.jsp
3. https://www.genome.jp/kegg/kegg1.html
4. https://reactome.org/PathwayBrowser/
5. Johnsen JI, Dyberg C and Wickström M (2019) Neuroblastoma—A Neural Crest Derived Embryonal 
Malignancy. Front. Mol. Neurosci. 12:9. doi: 10.3389/fnmol.2019.00009
6. Liu Z, Chen SS, Clarke S, Veschi V and Thiele CJ (2021) Targeting MYCN in Pediatric and Adult 
Cancers. Front. Oncol. 10:623679. doi: 10.3389/fonc.2020.623679
7. https://www.genecards.org/cgi-bin/carddisp.pl?gene=TP53TG3B
8. Sawyer, Iain A. (2018). Nuclear Architecture and Dynamics || Nuclear Bodies. , (), 235–
256. doi:10.1016/B978-0-12-803480-4.00010-7
9. Liu R, Shi P, Wang Z, Yuan C and Cui H (2021) Molecular Mechanisms of MYCN Dysregulation in 
Cancers. Front. Oncol. 10:625332. doi: 10.3389/fonc.2020.625332
10. Oshi M., Takahashi H., Tokumaru Y., Yan L., Rashid O.M., Matsuyama R., Endo I., Takabe K. G2M Cell 
Cycle Pathway Score as a Prognostic Biomarker of Metastasis in Estrogen Receptor (ER)-Positive 
Breast Cancer. Int. J. Mol. Sci. 2020;21:2921. doi: 10.3390/ijms21082921.
11. Oshi M., Takahashi H., Tokumaru Y., Yan L., Rashid O.M., Nagahashi M., Matsuyama R., Endo I., 
Takabe K. The E2F Pathway Score as a Predictive Biomarker of Response to Neoadjuvant Therapy in 
ER+/HER2- Breast Cancer. Cells. 2020;9:1643. doi: 10.3390/cells9071643.
12. https://www.genome.jp/pathway/hsa05208
13. https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000134184;r=1:109687814-
109709039
14. https://en.wikipedia.org/wiki/GSTT2
15. https://www.genome.jp/dbget-bin/www_bget?hsa:2953
16. https://www.genome.jp/entry/hsa:2944
17. Allocati, N., Masulli, M., Di Ilio, C. et al. Glutathione transferases: substrates, inihibitors and prodrugs in cancer and neurodegenerative diseases. Oncogenesis 7, 8 (2018). 
https://doi.org/10.1038/s41389-017-0025-3
18. Laborde, E. Glutathione transferases as mediators of signaling pathways involved in cell 
proliferation and cell death. Cell Death Differ 17, 1373–1380 (2010). 
https://doi.org/10.1038/cdd.2010.80
19. Singh, R.R.; Reindl, K.M. Glutathione S-Transferases in Cancer. Antioxidants 2021, 10, 701. 
https://doi.org/10.3390/antiox10050701
