# Description
APAtizer is a tool designed to analyse alternative polyadenylation of RNA-Seq data. Besides standard RNA-Seq data, APAtizer is also capable of performing analysis of 3'mRNA-Seq data and can also be used to analyse data obtained by other methods of sequencing than Illumina, for example, Ion Torrent. APAtizer has the option to analyse not only human but also mouse RNA-Seq data. Additionally, it is capable of performing differential gene expression, gene ontology analysis, visualizing Venn diagram intersections and Pearson correlation analysis. APAtizer is equipped to handle BAM files, DaPars files and htseq files. It is a user-friendly interface, that allows users to generate informative visualizations, including volcano plots, heatmaps, Venn intersections and gene lists. The APAtizer tool also provides the functionality to download the aforementioned plots and gene lists for further analysis and exploration. 

# Workflow
<img src="https://github.com/user-attachments/assets/7ae428b3-9d76-4ea4-b0ef-3619cf3a497b" alt="workflow">

# Required dependencies
1. samtools
2. snakemake
3. HTSeq
4. python3
5. gatk

# Installing dependencies
To install the required command line tools for the creation of the input files necessary to use APAtizer, the user must run the [install_dependencies_linux.sh](install_dependencies_linux.sh) script for linux or run the [install_dependencies_macos.sh](install_dependencies_macos.sh) script for macOS.
```shell
./install_dependencies_linux.sh #Linux
```
```shell
./install_dependencies_macos.sh #MacOS
```

# Required input files for APAtizer
1. Mapped reads .bam files (3'UTR-APA and IPA analysis with APAlyzer)
2. DaPars bedgraph files (3'UTR-APA analysis with DaPars)
3. HTSeq files (DGE analysis with DESeq2)


# Creating input files for APAtizer
The script to create the input files requires the raw BAM files to be placed in a folder called **RAW_BAM**. To start, clone the repository in the same directory of the **RAW_BAM** folder and enter the repository with the following command.

```shell
git clone https://github.com/brss12/APAtizer.git && cd APAtizer
```

After this, you will find a file called [create_inputs.sh](create_inputs.sh). Run it using the following command.

```shell
chmod +x create_inputs.sh && ./create_inputs.sh
```

This script will prompt the user to select the number (1-4) corresponding to the genome version used in the creation of the BAM files. This will be essencial for the creation of the input files necessary for APAtizer.

![image](https://github.com/user-attachments/assets/a0386d9d-9767-4bd5-be45-e16aaca687ad)

Upon selecting the number, the snakemake workflow scripts for the genome version chosen by the user will automatically run and create the input files necessary for the analysis in APAtizer. This script automatically sorts and removes the duplicates from the raw BAM files required for the APA analysis using the APAlyzer algorithm, and, using the BAM files, creates the DaPars bedgraph files required for 3'UTR-APA analysis employing the DaPars algorithm and it also creates the HTSeq files required for the DGE analysis using the DESeq2 package. All of these downstream analysis take place in the APAtizer's user interface.

Also, it is important to mention that depending on the size and ammount of BAM files, we recommend performing the aforementioned steps in a High Performance Computing (HPC) environment.

After finishing running, the following folder and files are created:

```plaintext
.
├── TRIMMED_READS/                           Folder containing the de-duplicated BAM files and the index files
│   ├── X.trim.bam                           De-duplicated BAM files
│   ├── X.trim.bam.bai                       Index files
├── TRIMMED_htseq/                           Folder containing the HTSeq files
│   ├── X.htseq.txt                          HTSeq files
├── DaPars_data/                             Folder containing the DaPars bedgraph files
│   ├── DaPars_data_result_temp.chrX.txt     DaPars bedgraph files
```

Also, it is important to mention that depending on the size and ammount of BAM files, we recommend performing the aforementioned steps in a High Performance Computing (HPC) environment.

***All done! Now you are ready to use APAtizer!***

To run the tool, the user can run the R script [APAtizer.R](APAtizer.R) using the following command or open the file on RStudio and press "Run App". Upon running the script, the R packages required for the utilization of the tool are all automatically installed.

```shell
Rscript APAtizer.R
```

Now, to showcase the capabilities of APAtizer, we performed two case studies using our tool. 

The first case was done on standard RNA-Seq data from 10 samples (5 sample pairs) from COAD (colon adenocarcinoma) from TCGA. The raw BAM files were obtained on the TCGA repository and were put in our snakemake workflow with [create_inputs.sh](create_inputs.sh). With this case study, our aim was to showcase that APAtizer is able to analyse standard RNA-Seq data.

The second case was done on 3'mRNA-Seq data from 4 samples from M1 macrophages, published in https://doi.org/10.3389/fimmu.2023.1182525, and uploaded in GEO (GSE163726). The FASTQ files were obtained and were previously aligned to the hg38 reference genome to obtain raw BAM files that were then put in our snakemake workflow with [create_inputs.sh](create_inputs.sh). Our aim with this case study was to showcase that APAtizer can not only work with data from standard RNA-Seq but also with 3'mRNA-Seq data.

Below, we showcase a walkthrough of the APAtizer tool showing the different tabs, inputs and outputs that can be obtained by the user.


# APAtizer walkthrough case study 1 (Illumina standard RNA-Seq samples from TCGA COAD)

## Sample Sheet
### Creating the sample sheet
In this section, the user may start by creating the sample sheet by clicking on the **Add row** button to add the necessary number of rows to construct the sample sheet. This sample sheet consists of two columns called **File.Name** and **Sample.Type**. The first column indicated the name of the BAM files and the second column indicated the name of the corresponding condition. An example of a sample sheet for this case study is shown below.

<img src="https://github.com/user-attachments/assets/071b9413-f2b9-41b9-b838-110fc5dd872a" alt="sample_sheet_case1">


## DaPars2
#### Inputs:
```plaintext
1. Select Multiple .txt files ---> Select all of the .txt files located in the **DaPars_data** folder
```
#### Outputs:
```plaintext
1. Len genes ---> List of genes undergoing 3'UTR-APA lengthening
2. Short genes ---> List of genes undergoing 3'UTR-APA shortening
```
The user can also search the gene lists for a specific gene of interest and download the said list using the download button below the search box.

### 3'UTR-APA gene list
<img src="https://github.com/user-attachments/assets/78e7e96d-4d46-4bcf-a47d-7213dd022b5a" alt="dapars_len_case1">


## APA_APALYZER
#### Inputs:
```plaintext
1. TRIMMED BAM files directory path ---> Paste the full path of the **TRIMMED_READS** folder
2. Select reference PAS ---> Select the reference PAS to use in the 3'UTR-APA analysis with APAlyzer
   2.1. hg19 ---> Reference PAS for hg19
   2.2. hg38 ---> Reference PAS for hg38
   2.3. mm9 ---> Reference PAS for mm9
   2.4. mm10 ---> Reference PAS for mm10
3. Select sequencing method ---> Select the sequening method used to obtain the reads
   3.1. Paired-end ---> For BAM files with paired-end reads
   3.2. Single-end ---> For BAM files with single-end reads
4. Select strandedness ---> Select the strandedness for the BAM files
   4.1. Forward stranded ---> For forward strand-specific BAM files  
   4.2. Reverse stranded ---> For reverse strand-specific BAM files
   4.3. Non-stranded ---> For non strand-specific BAM files
5. Select statistical test ---> Select the statistical test to employ in the analysis
   5.1. Unpaired t-test
   5.2. Paired t-test
   5.3. ANOVA
```
#### Outputs:
```plaintext
1. Select Output Type ---> Select the gene list to be shown on the output space
   1.1. NvsT_APA_UP ---> List of genes undergoing 3'UTR-APA lengthening
   1.2. NvsT_APA_DN ---> List of genes undergoing 3'UTR-APA shortening
   1.3. NvsT_APA_NC ---> List of non-significant genes
2. Selec Plot Type ---> Select the plot to be shown on the output space
   2.1. APA Volcano plot (top40) ---> Volcano plot for the 3'UTR-APA events with the top 40 most significant genes highlighted
   2.2. APA Volcano plot ---> Volcano plot for the 3'UTR-APA events
   2.3. APA Box ---> Box plot for the 3'UTR-APA events  
```

### 3'UTR-APA gene list
<img src="https://github.com/user-attachments/assets/e1754d6a-5c5b-4d5f-ba52-9dd071b99422" alt="apa_up_case1">


### 3'UTR-APA top-40 Volcano plot
<img src="https://github.com/user-attachments/assets/d80b2b47-7ca3-4338-bcb7-9b65d56874e1" alt="apa_top40_volcano_case1">


### 3'UTR-APA Volcano plot
<img src="https://github.com/user-attachments/assets/e09b6a32-e460-4892-92aa-fd7f680c9400" alt="apa_volcano_case1">


### 3'UTR-APA Box plot
<img src="https://github.com/user-attachments/assets/5d9b3122-5c95-48a8-8ecd-929a495e4c69" alt="apa_box_case1">


## IPA APALYZER
#### Inputs:
```plaintext
1. TRIMMED BAM files directory path ---> Paste the full path of the **TRIMMED_READS** folder
2. Select reference PAS ---> Select the reference PAS to use in the 3'UTR-APA analysis with APAlyzer
  2.1. hg19 ---> Reference PAS for hg19
  2.2. hg38 ---> Reference PAS for hg38
  2.3. mm9 ---> Reference PAS for mm9
  2.4. mm10 ---> Reference PAS for mm10
3. Select sequencing method ---> Select the sequening method used to obtain the reads
  3.1. Paired-end ---> For BAM files with paired-end reads
  3.2. Single-end ---> For BAM files with single-end reads
4. Select strandedness ---> Select the strandedness for the BAM files
  4.1. Forward stranded ---> For forward strand-specific BAM files  
  4.2. Reverse stranded ---> For reverse strand-specific BAM files
  4.3. Non-stranded ---> For non strand-specific BAM files
5. Select statistical test ---> Select the statistical test to employ in the analysis
  5.1. Unpaired t-test
  5.2. Paired t-test
  5.3. ANOVA
6. Number of threads ---> Select the number of threads to use for parallelization (Only available in IPA analysis with APAlyzer)
```
#### Outputs:
```plaintext
1. Select Output Type ---> Select the events list to be shown on the output space
  1.1. NvsT_IPA_events_UP ---> List of IPA upregulation events
  1.2. NvsT_IPA_events_DN ---> List of IPA downregulation events
  1.3. NvsT_IPA_events_NC ---> List of non-significant events
2. Select Output Type ---> Select the gene list to be shown on the output space
  2.1. NvsT_IPA_genes_UP ---> List of genes undergoing IPA upregulation
  2.2. NvsT_IPA_genes_DN ---> List of genes undergoing IPA downregulation
  2.3. NvsT_IPA_genes_NC ---> List of non-significant genes
3. Selec Plot Type ---> Select the plot to be shown on the output space
  3.1. IPA Volcano plot (top40) ---> Volcano plot for the IPA events with the top 40 most significant genes highlighted
  3.2. IPA Volcano plot ---> Volcano plot for the IPA events
  3.3. IPA Box ---> Box plot for the IPA events
```

### IPA events list
<img src="https://github.com/user-attachments/assets/a7a39299-ab10-4a02-b735-5b49d8d906c9" alt="ipa_events_up_case1">


### IPA gene list
<img src="https://github.com/user-attachments/assets/587766d4-ddb4-4469-9435-d40d4fa6e171" alt="ipa_genes_up_case1">


### IPA top-40 Volcano plot
<img src="https://github.com/user-attachments/assets/f9c699ae-5265-4e6e-a067-3d71207c9547" alt="ipa_top40_volcano_case1">


### IPA Volcano plot
<img src="https://github.com/user-attachments/assets/16f2a151-74d1-4529-a877-2810dc3843cc" alt="ipa_volcano_case1">


### IPA Box plot
<img src="https://github.com/user-attachments/assets/61cb91af-64c4-4a72-8c1b-7a5a06c8b82b" alt="ipa_box_case1">


## DGE

### DGE gene list
<img src="https://github.com/user-attachments/assets/d73702cd-b69e-42a9-84c8-5ba7859bfd83" alt="dge_up_case1">


### DGE PCA plot
<img src="https://github.com/user-attachments/assets/bf7696c4-40b7-4088-b5ee-423f370fe7e8" alt="dge_pca_case1">


### DGE Volcano plot
<img src="https://github.com/user-attachments/assets/74127a62-c414-4615-8018-2a52f1fb483e" alt="dge_volcano_case1">


### DGE Heatmap
<img src="https://github.com/user-attachments/assets/2e80c57c-b7fd-4a0b-add7-753e51b19133" alt="dge_heatmap_case1">


For the differential gene expression analysis, in the input space the user may paste the full path for the folder **TRIMMED_htseq/FILTERED/** that has all of the filtered htseq files. Once again, the user can also select the lists and the plots that will be displayed in the outputs. For the lists, the user can choose between DGE upregulated (*DGE_Genes_UP*), DGE downregulated (*DGE_Genes_DN*) and non-significant (*DGE_Genes_NC*). For the plots, the user can select a PCA plot (*PCA Plot*), a Volcano plot (*DGE Volcano Plot*) and a heatmap to evaluate the pattern of gene expression between conditions(*DGE Heatmap*).

In the output space, in the tab named *Number of DGE genes* the user can see a table showing the total number of upregulated, downregulated and non-significant genes. In *DGE_Genes* is where the full lists chosen before in the input space are displayed and the user can search those lists for a gene of interest and download them. In *Plots*, the same as in APAlyzer, is where one can visualize the selected plots and download them.

## GO_TERMS
### Biological Process (BP)
<img src="https://github.com/user-attachments/assets/3be1420e-4610-49f7-9e03-9622c06e1c0b" alt="go_bp_case1">


### Molecular Function (MF)
<img src="https://github.com/user-attachments/assets/40bdc6fc-7002-406f-a2d2-85f306d5f730" alt="go_mf_case1">


In this section, the user needs to select the list of genes in which he wants to perform gene ontology exploration, the organism database to be used (Human and Mouse) and the type of analysis to be performed, *Biological Process (BP)* and *Molecular Function (MF)*. The provided list should be one of the gene lists obtained in the previous steps (DaPars, APAlyzer or DGE analysis).

The output space only has one tab called *GO Plots* in which the resulting plots are displayed and can be downloaded.

## VENN DIAGRAMS
<img src="https://github.com/user-attachments/assets/52579063-3498-451f-8f6f-b637099fc220" alt="venn_case1">


For the Venn diagram intersections, the user can provide from 2 up to 5 gene lists obtained in the previous steps to execute the analysis. In the output section, the Venn diagram is displayed and can be downloaded in the tab *Venn Diagram*. Next, the user can obtain and download a list of the common genes between all the gene lists provided in the intersection.

## APA CORRELATION ANALYSIS
<img src="https://github.com/user-attachments/assets/d513b991-7b4c-4215-9fa0-de03abb87302" alt="apa_corr_case1">


In this section we can see the pearson correlation analysis scatter plot between the 3'UTR APA and DGE events. In this case study for M1 macrophages we can see that genes that undergo 3'UTR APA shortening events are being upregulated and the genes that undergo 3'UTR APA lengthening events are being downregulated.

## IPA CORRELATION ANALYSIS
<img src="https://github.com/user-attachments/assets/33a13a1b-938d-45b2-81b8-4ed8178e1383" alt="ipa_corr_case1">


# APAtizer walkthrough case study 2 (Illumina 3'mRNA-Seq samples from M1 Macrophages)

## Sample Sheet
### Creating the sample sheet
In this section, the user may start by creating the sample sheet by clicking on the **Add row** button to add the necessary number of rows to construct the sample sheet. This sample sheet consists of two columns called **File.Name** and **Sample.Type**. The first column indicated the name of the BAM files and the second column indicated the name of the corresponding condition. An example of a sample sheet for this case study is shown below.

<img src="https://github.com/user-attachments/assets/ec0553a0-fa0a-4b07-8b8c-9d4cbfe7487a" alt="sample_sheet_case2">


## DaPars2
### 3'UTR-APA gene list
<img src="https://github.com/user-attachments/assets/5daaca60-d5a4-47d5-9394-7c99ae9d1c88" alt="dapars_len_case2">

In this section, in the input space the user can select all the output files originated by the DaPars2 analysis that are located in the folder **DaPars_data**.

In the output space, we can observe the lists of genes that go through 3'UTR APA lengthening events (*Len genes*) and 3'UTR APA shortening events (*Short genes*). The user can also search the lists for a specific gene of interest and download the lists using the download button below the search box.

## APA_APALYZER
### 3'UTR-APA gene list
<img src="https://github.com/user-attachments/assets/fb9d43ee-1303-486e-acc0-494821155e78" alt="apa_up_case2">


### 3'UTR-APA top-40 Volcano plot
<img src="https://github.com/user-attachments/assets/4e70f8f6-4aef-47fa-afd3-d44aeebc9fa4" alt="apa_top40_volcano_case2">


### 3'UTR-APA Volcano plot
<img src="https://github.com/user-attachments/assets/89b184db-648d-4a56-9989-7eb0a40e4268" alt="apa_volcano_case2">


### 3'UTR-APA Box plot
<img src="https://github.com/user-attachments/assets/aa734dad-1e2e-4f01-afec-eab25ffd786a" alt="apa_box_case2">


In APA_APALYZER, in the input space the user may paste the full path of the folder **TRIMMED_READS** that contains all the de-duplicated BAM files and the index files. The user can also select the reference PAS for the analysis (hg19, hg38, mm9 and mm10), select the sequencing method (paired-end and single-end), select the strandedness of the BAM files (forward stranded, reverse stranded or non-stranded) and also select the statistical analysis to perform (unpaired t-test, paired t-test and ANOVA). The user may also select the output types of the analysis such as the lists that are displayed and the plots. For the lists, the user can choose between 3'UTR APA lengthening (*NvsT_APA_UP*), 3'UTR APA shortening (*NvsT_APA_DN*) and non-significant (*NvsT_APA_NC*). In the case of the plots, the choice is between a Volcano plot with the top 40 significant genes highlighted (*APA Volcano Plot (top40)*), the same plot but with no highlights (*APA Volcano Plot*) and a box plot (*APA Box*).

In the output space, in the tab called *Number of APA events* one can see a small table where the number of non significant, lengthening and shortening genes is present. In *NvsT_APA* the full lists chosen in the input space are presented to the user and he can search the list for a gene of interest and download it. Finally, in *Plots* the user can visualize the plots selected in the input space and download them as well. Below are some examples of how the chosen statistical analysis may affect the number of 3'UTR-APA events detected by the APAlyzer algorithm.


## IPA APALYZER
### IPA events list
<img src="https://github.com/user-attachments/assets/5e97da11-0994-4ab5-ae93-c370488f1306" alt="ipa_events_up_case2">


### IPA gene list
<img src="https://github.com/user-attachments/assets/4895c011-32bf-46e0-a447-13cc6f1c4bff" alt="ipa_genes_up_case2">


### IPA top-40 Volcano plot
<img src="https://github.com/user-attachments/assets/7a9d088f-83c6-4830-9090-4cadf156a342" alt="ipa_top40_volcano_case2">


### IPA Volcano plot
<img src="https://github.com/user-attachments/assets/071f362e-987a-4ed1-8c75-5205616b0e0c" alt="ipa_volcano_case2">


### IPA Box plot
<img src="https://github.com/user-attachments/assets/decfc250-bd4b-4fdb-8535-43233f9fb22e" alt="ipa_box_case2">


In IPA_APALYZER, in the input space the user may paste the full path of the folder **TRIMMED_READS** that contains all the de-duplicated BAM files and the index files. The user can also select the reference PAS for the analysis (hg19, hg38, mm9 and mm10), select the sequencing method (paired-end and single-end), select the strandedness of the BAM files (forward stranded, reverse stranded or non-stranded) and also, specifically for this analysis, select the number of threads to be used in the analysis. The user then selects the output types of the analysis such as the lists that are displayed and the plots. For the lists, the user can choose between IPA upregulated events (*NvsT_IPA_events_UP*), IPA downregulated events (*NvsT_IPA_events_DN*) and non-significant events (*NvsT_IPA_events_NC*). Another output type is the gene lists with the unique genes such as IPA upregulated genes (*NvsT_IPA_genes_UP*), IPA downregulated genes (*NvsT_IPA_genes_DN*) and non-significant genes (*NvsT_IPA_genes_NC*). In the case of the plots, the choice is between a Volcano plot with the top 40 significant genes highlighted (*IPA Volcano Plot (top40)*), the same plot but with no highlights (*IPA Volcano Plot*) and a box plot (*IPA Box*).

In the output space, in the tab called *Number of IPA events* one can see a small table where the number of non significant, lengthening and shortening IPA events are present. In *NvsT_IPA_events* the full lists of the IPA events are presented to the user and he can search the list for a gene of interest and download it. In *NvsT_IPA_genes* the full list for the unique genes is presented to the user and he can, also, search for a gene of interest and download it. Finally, in *Plots* the user can visualize the plots selected in the input space and download them as well.

## DGE
### DGE gene list
<img src="https://github.com/user-attachments/assets/d4fa8c84-5ba3-4c53-a10e-d3c0efbb40cf" alt="dge_up_case2">


### DGE PCA plot
<img src="https://github.com/user-attachments/assets/4546af7f-8efe-406a-a2af-eccc494c47ec" alt="dge_pca_case2">


### DGE Volcano plot
<img src="https://github.com/user-attachments/assets/f8028476-114f-4688-b580-8bd2d5763f61" alt="dge_volcano_case2">


### DGE Heatmap
<img src="https://github.com/user-attachments/assets/fb263765-c198-478a-b711-3c9fd9e375dc" alt="dge_heatmap_case2">


For the differential gene expression analysis, in the input space the user may paste the full path for the folder **TRIMMED_htseq/FILTERED/** that has all of the filtered htseq files. Once again, the user can also select the lists and the plots that will be displayed in the outputs. For the lists, the user can choose between DGE upregulated (*DGE_Genes_UP*), DGE downregulated (*DGE_Genes_DN*) and non-significant (*DGE_Genes_NC*). For the plots, the user can select a PCA plot (*PCA Plot*), a Volcano plot (*DGE Volcano Plot*) and a heatmap to evaluate the pattern of gene expression between conditions(*DGE Heatmap*).

In the output space, in the tab named *Number of DGE genes* the user can see a table showing the total number of upregulated, downregulated and non-significant genes. In *DGE_Genes* is where the full lists chosen before in the input space are displayed and the user can search those lists for a gene of interest and download them. In *Plots*, the same as in APAlyzer, is where one can visualize the selected plots and download them.

## GO_TERMS
### Biological Process (BP)
<img src="https://github.com/user-attachments/assets/8e338810-fb2d-41fe-963a-f08cde2ec5b5" alt="go_bp_case2">


### Molecular Function (MF)
<img src="https://github.com/user-attachments/assets/d35ecfd5-680d-4983-9c33-2555ffd273c9" alt="go_mf_case2">


In this section, the user needs to select the list of genes in which he wants to perform gene ontology exploration, the organism database to be used (Human and Mouse) and the type of analysis to be performed, *Biological Process (BP)* and *Molecular Function (MF)*. The provided list should be one of the gene lists obtained in the previous steps (DaPars, APAlyzer or DGE analysis).

The output space only has one tab called *GO Plots* in which the resulting plots are displayed and can be downloaded.

## VENN DIAGRAMS
<img src="https://github.com/user-attachments/assets/e06d0b77-7509-451f-b504-83377452a0ef" alt="venn_case2">


For the Venn diagram intersections, the user can provide from 2 up to 5 gene lists obtained in the previous steps to execute the analysis. In the output section, the Venn diagram is displayed and can be downloaded in the tab *Venn Diagram*. Next, the user can obtain and download a list of the common genes between all the gene lists provided in the intersection.

## APA CORRELATION ANALYSIS
<img src="https://github.com/user-attachments/assets/2fd4c20b-bd9a-4aa8-a15e-7f684dd9af09" alt="apa_corr_case2">


In this section we can see the pearson correlation analysis scatter plot between the 3'UTR-APA and DGE events. In this case study for M1 macrophages we can see that genes that undergo 3'UTR-APA shortening events are being upregulated and the genes that undergo 3'UTR APA lengthening events are being downregulated.

## IPA CORRELATION ANALYSIS
<img src="https://github.com/user-attachments/assets/28027bf2-de76-4398-9a38-16ac5855c82b" alt="ipa_corr_case2">


Now, in this section we have the pearson correlation analysis scatter plot between IPA and DGE events. In this case study for M1 macrophages we can also see that we have a significant negative correlation between IPA and DGE events. Genes that undergo IPA downregulation are being more expressed, whereas genes that undergo IPA upregulation are being less expressed.


# Final remarks
In this README, three case studies were used to demonstrate and explain the features and capabilities of the APAtizer tool. With this tool, the user can analyze RNA-Seq data from various sources and retrieve many plots and useful information regarding 3'UTR-APA & IPA events via DaPars2 and APAlyzer analysis, DGE via DESeq2, the function of those genes via GO analysis, the common genes between gene lists using Venn diagram intersections and the correlation between 3'UTR-APA & IPA events and DGE via Pearson correlation analysis scatter plots.
