# Description
APAtizer is a tool designed to analyse alternative polyadenylation of RNA-Seq data. Besides standard RNA-Seq data, APAtizer is also capable of performing analysis of 3'mRNA-Seq data and can also be used to analyse data obtained by other methods of sequencing than Illumina, for example, Ion Torrent. APAtizer has the option to analyse not only human but also mouse RNA-Seq data. Additionally, it is capable of performing differential gene expression, gene ontology analysis, visualizing Venn diagram intersections and Pearson correlation analysis. APAtizer is equipped to handle BAM files, DaPars files and htseq files. It is a user-friendly interface, that allows users to generate informative visualizations, including volcano plots, heatmaps, Venn intersections and gene lists. The APAtizer tool also provides the functionality to download the aforementioned plots and gene lists for further analysis and exploration. 

# Workflow
<img src="https://github.com/user-attachments/assets/7ae428b3-9d76-4ea4-b0ef-3619cf3a497b" alt="workflow">


# Installing dependencies
To install the required command line tools for the creation of the input files necessary to use APAtizer, the user must run the [install_dependencies_linux.sh](install_dependencies_linux.sh) script for linux or run the [install_dependencies_macos.sh](install_dependencies_macos.sh) script for macOS.
```shell
./install_dependencies_linux.sh #Linux
```
```shell
./install_dependencies_macos.sh #MacOS
```

# Creating the input files
The script to create the input files requires the raw BAM files to be placed in a folder called **RAW_BAM**. To start, clone the repository in the same directory of the **RAW_BAM** and enter the folder with the following commands.

```shell
git clone https://github.com/brss12/APAtizer.git && cd APAtizer
```

After this, you will find a file called [create_inputs.sh](create_inputs.sh). Run it using the following command.

```shell
chmod +x create_inputs.sh && ./create_inputs.sh
```

This script will prompt the user to select the number (1-4) corresponding to the genome version used in the creation of the BAM files. This will be essencial for the creation of the input files necessary for APAtizer.

![image](https://github.com/user-attachments/assets/a0386d9d-9767-4bd5-be45-e16aaca687ad)

Upon selecting the number, the snakemake workflow scripts for the genome version chosen by the user will automatically run and create the input files necessary for the analysis in APAtizer. This script automatically sorts and removes the duplicates from the raw BAM files required for the APA analysis using the APAlyzer algorithm, creates the DaPars txt files required for APA analysis employing the DaPars algorithm and it also creates the htseq files required for the DGE analysis using the DESeq2 package. All of these downstream analysis take place in the APAtizer's user interface.

After finishing running, the folders **TRIMMED_READS**, **TRIMMED_QC**, **TRIMMED_htseq** and **DaPars_data** are created. The **TRIMMED_READS** folder is where the de-duplicated BAM files along with the corresponding BAI index files are located. The **TRIMMED_QC** folder is where the fastqc reports of the de-duplicated BAM files are located for the user to obtain information about the number of reads, length of reads and many more parameters about the BAM files. The **TRIMMED_htseq** folder is where the htseq files for the DGE analysis are located. Finally, the **DaPars_data** folder is where the txt files necessary for the DaPars analysis are located.

Also, it is important to mention that depending on the size and ammount of BAM files, we recommend performing the aforementioned steps in a High Performance Computing (HPC) environment.

***All done! Now you are ready to use APAtizer!***

To run the tool, the user can run the R script [APAtizer.R](APAtizer.R) using the following command or open the file on RStudio and press "Run App". Upon running the script, the R packages required for the utilization of the tool are all automatically installed.

```shell
Rscript APAtizer.R
```

Now, to showcase the capabilities of APAtizer, we performed three case studies using our tool. 

The first case was done on standard RNA-Seq data from 8 samples (4 samples from Heart and 4 samples from Testis) of Mouse retrieved from GEO (GSM900199 and GSM900193 accession numbers). The FASTQ files were obtained via Illumina sequencing technologies and were aligned to the mm9 reference genome to obtain raw BAM files that were processed via our snakemake workflow to create the inputs necessary for APAtizer. Our aim with this case study was to show to the users that, in addition to Human RNA-Seq data, APAtizer can also work with RNA-Seq data derived from Mouse.

The second case was done on 3'mRNA-Seq data from 4 samples from M1 macrophages, published in https://doi.org/10.3389/fimmu.2023.1182525. The FASTQ files were obtained via Illumina sequencing technologies and were previously aligned to the hg38 reference genome to obtain raw BAM files that were put in our snakemake workflow with [create_inputs.sh](create_inputs.sh). Our aim with this case study was to showcase that APAtizer can not only work with data from standard RNA-Seq but also with 3'mRNA-Seq data.

The third case study was done on standard RNA-Seq data from 4 samples (2 samples from N-Diethylnitrosamine (DEN) treated WT and 2 samples from WT) of Mouse also retrieved from NCBI (PRJNA214241 BioProject). The FASTQ files were obtained via Ion Torrent sequencing technologies and were once again aligned to the mm9 reference genome to obtain raw BAM files that were processed using our snakemake workflow script to create the inputs necessary for APAtizer. Our aim with this final case study was to showcase that the APAtizer tool, while working with Illumina sequencing technologies, can also work with Ion Torrent sequencing technologies.

Below, we showcase a walkthrough of the APAtizer tool showing the different tabs, inputs and outputs that can be obtained by the user.


# APAtizer walkthrough case study 1 (Illumina 3'mRNA-Seq samples from M1 Macrophages)

## Sample Sheet
### Creating the sample sheet
In this section, the user may start by creating the sample sheet by clicking on the **Add row** button to add the necessary number of rows to construct the sample sheet. This sample sheet consists of two columns called **File.Name** and **Sample.Type**. The first column indicated the name of the BAM files and the second column indicated the name of the corresponding condition. An example of a sample sheet for this case study is shown below.

<img src="https://github.com/user-attachments/assets/ec0553a0-fa0a-4b07-8b8c-9d4cbfe7487a" alt="sample_sheet_case1">


## DaPars2
### 3'UTR-APA gene list
<img src="https://github.com/user-attachments/assets/5daaca60-d5a4-47d5-9394-7c99ae9d1c88" alt="dapars_len_case1">

In this section, in the input space the user can select all the output files originated by the DaPars2 analysis that are located in the folder **DaPars_data**.

In the output space, we can observe the lists of genes that go through 3'UTR APA lengthening events (*Len genes*) and 3'UTR APA shortening events (*Short genes*). The user can also search the lists for a specific gene of interest and download the lists using the download button below the search box.

## APA_APALYZER
### 3'UTR-APA gene list
<img src="https://github.com/user-attachments/assets/fb9d43ee-1303-486e-acc0-494821155e78" alt="apa_up_case1">


### 3'UTR-APA top-40 Volcano plot
<img src="https://github.com/user-attachments/assets/4e70f8f6-4aef-47fa-afd3-d44aeebc9fa4" alt="apa_top40_volcano_case1">


### 3'UTR-APA Volcano plot
<img src="https://github.com/user-attachments/assets/89b184db-648d-4a56-9989-7eb0a40e4268" alt="apa_volcano_case1">


### 3'UTR-APA Box plot
<img src="https://github.com/user-attachments/assets/aa734dad-1e2e-4f01-afec-eab25ffd786a" alt="apa_box_case1">


In APA_APALYZER, in the input space the user may paste the full path of the folder **TRIMMED_READS** that contains all the de-duplicated BAM files and the index files. The user can also select the reference PAS for the analysis (hg19, hg38, mm9 and mm10), select the sequencing method (paired-end and single-end), select the strandedness of the BAM files (forward stranded, reverse stranded or non-stranded) and also select the statistical analysis to perform (unpaired t-test, paired t-test and ANOVA). The user may also select the output types of the analysis such as the lists that are displayed and the plots. For the lists, the user can choose between 3'UTR APA lengthening (*NvsT_APA_UP*), 3'UTR APA shortening (*NvsT_APA_DN*) and non-significant (*NvsT_APA_NC*). In the case of the plots, the choice is between a Volcano plot with the top 40 significant genes highlighted (*APA Volcano Plot (top40)*), the same plot but with no highlights (*APA Volcano Plot*) and a box plot (*APA Box*).

In the output space, in the tab called *Number of APA events* one can see a small table where the number of non significant, lengthening and shortening genes is present. In *NvsT_APA* the full lists chosen in the input space are presented to the user and he can search the list for a gene of interest and download it. Finally, in *Plots* the user can visualize the plots selected in the input space and download them as well. Below are some examples of how the chosen statistical analysis may affect the number of 3'UTR-APA events detected by the APAlyzer algorithm.

## Unpaired t-test
<img src="https://github.com/user-attachments/assets/18859a11-c48a-457f-84ce-8469b70f0d0e" alt="apa_events_unpaired_case1">


## Paired t-test
<img src="https://github.com/user-attachments/assets/da1f0cb0-c2dd-4c6f-bef7-7f5d77a540b9" alt="apa_events_paired_case1">


## ANOVA
<img src="https://github.com/user-attachments/assets/944f68b7-c304-4fde-9b12-0b0079aa6445" alt="apa_events_anova_case1">


## IPA APALYZER
### IPA events list
<img src="https://github.com/user-attachments/assets/5e97da11-0994-4ab5-ae93-c370488f1306" alt="ipa_events_up_case1">


### IPA gene list
<img src="https://github.com/user-attachments/assets/4895c011-32bf-46e0-a447-13cc6f1c4bff" alt="ipa_genes_up_case1">


### IPA top-40 Volcano plot
<img src="https://github.com/user-attachments/assets/7a9d088f-83c6-4830-9090-4cadf156a342" alt="ipa_top40_volcano_case1">


### IPA Volcano plot
<img src="https://github.com/user-attachments/assets/071f362e-987a-4ed1-8c75-5205616b0e0c" alt="ipa_volcano_case1">


### IPA Box plot
<img src="https://github.com/user-attachments/assets/decfc250-bd4b-4fdb-8535-43233f9fb22e" alt="ipa_box_case1">


In IPA_APALYZER, in the input space the user may paste the full path of the folder **TRIMMED_READS** that contains all the de-duplicated BAM files and the index files. The user can also select the reference PAS for the analysis (hg19, hg38, mm9 and mm10), select the sequencing method (paired-end and single-end), select the strandedness of the BAM files (forward stranded, reverse stranded or non-stranded) and also, specifically for this analysis, select the number of threads to be used in the analysis. The user then selects the output types of the analysis such as the lists that are displayed and the plots. For the lists, the user can choose between IPA upregulated events (*NvsT_IPA_events_UP*), IPA downregulated events (*NvsT_IPA_events_DN*) and non-significant events (*NvsT_IPA_events_NC*). Another output type is the gene lists with the unique genes such as IPA upregulated genes (*NvsT_IPA_genes_UP*), IPA downregulated genes (*NvsT_IPA_genes_DN*) and non-significant genes (*NvsT_IPA_genes_NC*). In the case of the plots, the choice is between a Volcano plot with the top 40 significant genes highlighted (*IPA Volcano Plot (top40)*), the same plot but with no highlights (*IPA Volcano Plot*) and a box plot (*IPA Box*).

In the output space, in the tab called *Number of IPA events* one can see a small table where the number of non significant, lengthening and shortening IPA events are present. In *NvsT_IPA_events* the full lists of the IPA events are presented to the user and he can search the list for a gene of interest and download it. In *NvsT_IPA_genes* the full list for the unique genes is presented to the user and he can, also, search for a gene of interest and download it. Finally, in *Plots* the user can visualize the plots selected in the input space and download them as well.

## DGE
### DGE gene list
<img src="https://github.com/user-attachments/assets/d4fa8c84-5ba3-4c53-a10e-d3c0efbb40cf" alt="dge_up_case1">


### DGE PCA plot
<img src="https://github.com/user-attachments/assets/4546af7f-8efe-406a-a2af-eccc494c47ec" alt="dge_pca_case1">


### DGE Volcano plot
<img src="https://github.com/user-attachments/assets/f8028476-114f-4688-b580-8bd2d5763f61" alt="dge_volcano_case1">


### DGE Heatmap
<img src="https://github.com/user-attachments/assets/fb263765-c198-478a-b711-3c9fd9e375dc" alt="dge_heatmap_case1">


For the differential gene expression analysis, in the input space the user may paste the full path for the folder **TRIMMED_htseq/FILTERED/** that has all of the filtered htseq files. Once again, the user can also select the lists and the plots that will be displayed in the outputs. For the lists, the user can choose between DGE upregulated (*DGE_Genes_UP*), DGE downregulated (*DGE_Genes_DN*) and non-significant (*DGE_Genes_NC*). For the plots, the user can select a PCA plot (*PCA Plot*), a Volcano plot (*DGE Volcano Plot*) and a heatmap to evaluate the pattern of gene expression between conditions(*DGE Heatmap*).

In the output space, in the tab named *Number of DGE genes* the user can see a table showing the total number of upregulated, downregulated and non-significant genes. In *DGE_Genes* is where the full lists chosen before in the input space are displayed and the user can search those lists for a gene of interest and download them. In *Plots*, the same as in APAlyzer, is where one can visualize the selected plots and download them.

## GO_TERMS
### Biological Process (BP)
<img src="https://github.com/user-attachments/assets/8e338810-fb2d-41fe-963a-f08cde2ec5b5" alt="go_bp_case1">


### Molecular Function (MF)
<img src="https://github.com/user-attachments/assets/d35ecfd5-680d-4983-9c33-2555ffd273c9" alt="go_mf_case1">


In this section, the user needs to select the list of genes in which he wants to perform gene ontology exploration, the organism database to be used (Human and Mouse) and the type of analysis to be performed, *Biological Process (BP)* and *Molecular Function (MF)*. The provided list should be one of the gene lists obtained in the previous steps (DaPars, APAlyzer or DGE analysis).

The output space only has one tab called *GO Plots* in which the resulting plots are displayed and can be downloaded.

## VENN DIAGRAMS
<img src="https://github.com/user-attachments/assets/e06d0b77-7509-451f-b504-83377452a0ef" alt="venn_case1">


For the Venn diagram intersections, the user can provide from 2 up to 5 gene lists obtained in the previous steps to execute the analysis. In the output section, the Venn diagram is displayed and can be downloaded in the tab *Venn Diagram*. Next, the user can obtain and download a list of the common genes between all the gene lists provided in the intersection.

## APA CORRELATION ANALYSIS
<img src="https://github.com/user-attachments/assets/2fd4c20b-bd9a-4aa8-a15e-7f684dd9af09" alt="apa_corr_case1">


In this section we can see the pearson correlation analysis scatter plot between the 3'UTR APA and DGE events. In this case study for M1 macrophages we can see that genes that undergo 3'UTR APA shortening events are being upregulated and the genes that undergo 3'UTR APA lengthening events are being downregulated.

## IPA CORRELATION ANALYSIS
<img src="https://github.com/user-attachments/assets/28027bf2-de76-4398-9a38-16ac5855c82b" alt="ipa_corr_case1">


Now, in this section we have the pearson correlation analysis scatter plot between IPA and DGE events. In this case study for M1 macrophages we can also see that we have a significant negative correlation between IPA and DGE events. Genes that undergo IPA downregulation are being more expressed, whereas genes that undergo IPA upregulation are being less expressed.


# Final remarks
In this README, three case studies were used to demonstrate and explain the features and capabilities of the APAtizer tool. With this tool, the user can analyze RNA-Seq data from various sources and retrieve many plots and useful information regarding 3'UTR-APA & IPA events via DaPars2 and APAlyzer analysis, DGE via DESeq2, the function of those genes via GO analysis, the common genes between gene lists using Venn diagram intersections and the correlation between 3'UTR-APA & IPA events and DGE via Pearson correlation analysis scatter plots.
