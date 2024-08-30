# Description
APAtizer is a tool designed to analyse alternative polyadenylation of RNA-Seq data. Additionally, it is capable of performing differential gene expression, gene ontology analysis, visualizing Venn diagram intersections and Pearson correlation analysis. APAtizer is equipped to handle BAM files, DaPars txt files and htseq files. It is a user-friendly interface, that allows users to generate informative visualizations, including volcano plots, heatmaps, Venn intersections and gene lists. The APAtizer tool also provides the functionality to download the aforementioned plots and gene lists for further analysis and exploration. 

# Workflow
![image](https://github.com/brss12/APAtizer/assets/121204829/02d6eb3a-6bd1-47f1-9c40-b0a9ec19af1d)


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

This script will prompt the user to select the number (1-4) corresponding to the genome version using in the alignment of the FASTQ files to create the BAM files. This will be essencial for the creation of the input files necessary for APAtizer.

![image](https://github.com/user-attachments/assets/a0386d9d-9767-4bd5-be45-e16aaca687ad)

Upon selecting the number, the snakemake workflow scripts for the genome version chosen by the user, will automatically run and create the input files necessary for the analysis in APAtizer. This script automatically sorts and removes the duplicates from the raw BAM files required for the APA analysis using the APAlyzer algorithm, creates the DaPars txt files required for APA analysis employing the DaPars algorithm and it also created the htseq files required for the differential gene expression analysis feature of the tool. All of these downstream analysis take place in the APAtizer's user interface.

After finishing running, creates two new folders called **SORTED_BAM** and **TRIMMED_READS**. In the **SORTED_BAM** folder, is where the sorted bam files along with fastqc reports are placed and in the **TRIMMED_READS** is where the de-duplicated bam files along with fastqc reports are placed. After this, a folder called **DaPars_data** is created with the final DaPars txt files.

Depending on the size and ammount of BAM files, we recommend running both scripts in a HPC environment. With the following command, run the first script on the same directory of the **RAW_BAM** folder to create the sorted and de-duplicated BAM files and the htseq files.
```shell
nohup bash -c "snakemake -s APAtizer_inputs_1.smk --cores <n_cores> --wait-for-files" & #The first script can be run using parallelization
```
After this, to create the required DaPars2 files, run the second script on the same directory of the **RAW_BAM** folder using the following command.
```shell
nohup bash -c "snakemake -s APAtizer_inputs_2.smk --cores 1 --wait-for-files" & #The second script must be run sequentially
```

***All done! Now you are ready to use APAtizer!***

# APAtizer walkthrough
In the walkthrough of the APAtizer tool, a case study using 46 samples (23 Normal Tissue samples and 23 Primary Tumor samples) from colorectal cancer from TCGA was performed to showcase the interface and the capabilities of the tool all while retrieving useful insights from the data.

To run the tool, the user must run the R script [APAtizer.R](APAtizer.R) using the following command.
```shell
Rscript APAtizer.R
```
Upon running the command, all the R packages needed for APAtizer will be installed automatically and a link to the tool will appear in the end. Click said link to open APAtizer in your default web browser!
The user can also open the R script in Rstudio and then click on the ***Run App*** button to launch the tool in a separate window.

## DaPars2 interface
### 3'UTR APA lengthening genes
<img src="https://github.com/brss12/APAtizer/assets/121204829/4d32cb27-7bfb-4add-894a-a64b9b2a148d" alt="DAPARS_LEN">

### 3'UTR APA shortening genes
<img src="https://github.com/brss12/APAtizer/assets/121204829/7c57fff7-d9f0-4f80-afe6-cb03c3dd9081" alt="DAPARS_SHORT">

In this section, in the input space the user can select the 24 output files originated by the DaPars2 analysis that are located in the folder **DaPars_data** and the TCGA sample sheet.

In the output space, we can observe the lists of genes that go through 3'UTR APA lengthening events (*Len genes*) and 3'UTR APA shortening events (*Short genes*). The user can also search the lists for a specific gene of interest and download the lists using the download button below the search box.

## APA_APALYZER interface
### 3'UTR APA lengthening genes
<img src="https://github.com/brss12/APAtizer/assets/121204829/9f48a0a1-d5d6-497b-a841-b4165536a1ff" alt="APA_UP">

### 3'UTR APA non-significant genes
<img src="https://github.com/brss12/APAtizer/assets/121204829/ff179f5e-22ac-4bed-97d6-5d7dbfce2aa6" alt="APA_NC">

### 3'UTR APA shortening genes
<img src="https://github.com/brss12/APAtizer/assets/121204829/702be557-57b8-4adb-8b5f-ba2b251d6c34" alt="APA_DN">

### 3'UTR APA top-40 Volcano plot
<img src="https://github.com/brss12/APAtizer/assets/121204829/418a9875-556e-4222-adcf-e043665ec355" alt="APA_top40_Volcano">

### 3'UTR APA Volcano plot
<img src="https://github.com/brss12/APAtizer/assets/121204829/45510c74-f9e9-4a39-a182-044f836d3a84" alt="APA_Volcano">

### 3'UTR APA Box plot
<img src="https://github.com/brss12/APAtizer/assets/121204829/aba89b3d-a60f-4fb3-ae20-583cf9de528b" alt="APA_Box">

In APA_APALYZER, in the input space the user may paste the full path of the folder **TRIMMED_READS** that contains all the de-duplicated BAM files and the index files and select the TCGA sample sheet. The user may also select the output types of the analysis such as the lists that are displayed and the plots. For the lists, the user can choose between 3'UTR APA lengthening (*NvsT_APA_UP*), 3'UTR APA shortening (*NvsT_APA_DN*) and non-significant (*NvsT_APA_NC*). In the case of the plots, the choice is between a Volcano plot with the top 40 significant genes highlighted (*APA Volcano Plot (top40)*), the same plot but with no highlights (*APA Volcano Plot*) and a box plot (*APA Box*).

In the output space, in the tab called *Number of APA events* one can see a small table where the number of non significant, lengthening and shortening genes is present. In *NvsT_APA* the full lists chosen in the input space are presented to the user and he can search the list for a gene of interest and download it. Finally, in *Plots* the user can visualize the plots selected in the input space and download them as well.

## IPA APALYZER interface
### IPA upregulated events
<img src="https://github.com/brss12/APAtizer/assets/121204829/cc32470c-97f6-41ab-b98f-56098dc17a87" alt="IPA_events_UP">

### IPA non-significant events
<img src="https://github.com/brss12/APAtizer/assets/121204829/b96a6383-a894-4216-84f4-dc333420d806" alt="IPA_events_NC">

### IPA downregulated events
<img src="https://github.com/brss12/APAtizer/assets/121204829/2c84e83c-8d76-42ac-9a92-0a81ba4f4383" alt="IPA_events_DN">

### IPA upregulated genes
<img src="https://github.com/brss12/APAtizer/assets/121204829/382adc6e-2b15-49a3-a41c-3a54f6947833" alt="IPA_genes_UP">

### IPA non-significant genes
<img src="https://github.com/brss12/APAtizer/assets/121204829/fd2026f3-493b-4039-925c-397ac2edcf2e" alt="IPA_genes_NC">

### IPA downregulated genes
<img src="https://github.com/brss12/APAtizer/assets/121204829/1c164835-eb12-41e2-a382-9bcaec2aa992" alt="IPA_genes_DN">

### IPA top-40 Volcano plot
<img src="https://github.com/brss12/APAtizer/assets/121204829/77235954-264f-48a0-9818-76e529e531c5" alt="IPA_top40_Volcano">

### IPA Volcano plot
<img src="https://github.com/brss12/APAtizer/assets/121204829/62c236c0-499c-43ef-b980-224a720f50af" alt="IPA_Volcano">

### IPA Box plot
<img src="https://github.com/brss12/APAtizer/assets/121204829/9dec095b-0927-48d9-90fc-eee29456b0f6" alt="IPA_Box">

In IPA_APALYZER, in the input space the user may paste the full path of the folder **TRIMMED_READS** that contains all the de-duplicated BAM files and the index files and select the TCGA sample sheet. The user then selects the output types of the analysis such as the lists that are displayed and the plots. For the lists, the user can choose between IPA upregulated events (*NvsT_IPA_events_UP*), IPA downregulated events (*NvsT_IPA_events_DN*) and non-significant events (*NvsT_IPA_events_NC*). Another output type is the gene lists with the unique genes such as IPA upregulated genes (*NvsT_IPA_genes_UP*), IPA downregulated genes (*NvsT_IPA_genes_DN*) and non-significant genes (*NvsT_IPA_genes_NC*). In the case of the plots, the choice is between a Volcano plot with the top 40 significant genes highlighted (*IPA Volcano Plot (top40)*), the same plot but with no highlights (*IPA Volcano Plot*) and a box plot (*IPA Box*).

In the output space, in the tab called *Number of IPA events* one can see a small table where the number of non significant, lengthening and shortening IPA events are present. In *NvsT_IPA_events* the full lists of the IPA events are presented to the user and he can search the list for a gene of interest and download it. In *NvsT_IPA_genes* the full list for the unique genes is presented to the user and he can, also, search for a gene of interest and download it. Finally, in *Plots* the user can visualize the plots selected in the input space and download them as well.

## DGE interface
### DGE upregulated genes
<img src="https://github.com/brss12/APAtizer/assets/121204829/cc59faa6-3ce2-45f8-b89c-7515a2339f5e" alt="DGE_UP">

### DGE non-significant genes
<img src="https://github.com/brss12/APAtizer/assets/121204829/1ea5a5b8-6fd0-42fd-9ebd-72815daeebc4" alt="DGE_NC">

### DGE downregulated genes
<img src="https://github.com/brss12/APAtizer/assets/121204829/fae6ce61-6d3a-4829-9f62-7361f297b7ec" alt="DGE_DN">

### DGE PCA plot
<img src="https://github.com/brss12/APAtizer/assets/121204829/91d89676-9c49-4c98-8131-7b9ea4b92aa0" alt="DGE_PCA">

### DGE Volcano plot
<img src="https://github.com/brss12/APAtizer/assets/121204829/78fe38b1-1cf9-4e0d-a368-fc45d2f1f435" alt="DGE_Volcano">

### DGE Heatmap
<img src="https://github.com/brss12/APAtizer/assets/121204829/27bec740-51a8-4604-bbe9-6920129ca7d9" alt="DGE_Heatmap">

For the differential gene expression analysis, in the input space the user may paste the full path for the folder **FILTERED** that has all of the htseq files and select the TCGA sample sheet. Once again, the user can also select the lists and the plots that will be displayed in the outputs. For the lists, the user can choose between DGE upregulated (*DGE_Genes_UP*), DGE downregulated (*DGE_Genes_DN*) and non-significant (*DGE_Genes_NC*). For the plots, the user can select a PCA plot (*PCA Plot*), a Volcano plot (*DGE Volcano Plot*) and a heatmap to evaluate the pattern of gene expression between conditions(*DGE Heatmap*).

In the output space, in the tab named *Number of DGE genes* the user can see a table showing the total number of upregulated, downregulated and non-significant genes. In *DGE_Genes* is where the full lists chosen before in the input space are displayed and the user can search those lists for a gene of interest and download them. In *Plots*, the same as in APAlyzer, is where one can visualize the selected plots and download them.

## GO_TERMS interface
### Biological Process (BP)
<img src="https://github.com/brss12/APAtizer/assets/121204829/92c7946b-add8-4b7a-aae4-3b18e71d4cf9" alt="GO_BP">

### Molecular Function (MF)
<img src="https://github.com/brss12/APAtizer/assets/121204829/716fa7b8-3ffd-4a0e-a5d0-0dffbc1d795c" alt="GO_MF">

In this section, the user only needs to select the list of genes in which he wants to perform gene ontology exploration and the type of analysis to be performed, *Biological Process (BP)* and *Molecular Function (MF)*. The provided list should be one of the gene lists obtained in the previous steps (DaPars, APAlyzer or DGE analysis).

The output space only has one tab called *GO Plots* in which the resulting plots are displayed and can be downloaded.

## VENN DIAGRAMS interface
<img src="https://github.com/brss12/APAtizer/assets/121204829/0faebdf0-1baa-4be4-9d0b-75aab8cf99d3" alt="Venn_2">

For the Venn diagram intersections, the user can provide up to 5 gene lists obtained in the previous steps to execute the analysis. In the output section, the Venn diagram is displayed and can be downloaded in the tab *Venn Diagram*. Next, the user can obtain and download a list of the common genes between all the gene lists provided in the intersection.

## APA CORRELATION ANALYSIS interface
<img src="https://github.com/brss12/APAtizer/assets/121204829/66f1f08b-3795-454b-855a-652ce2080572" alt="APA_CORR">

In this section we can see the pearson correlation analysis scatter plot between the 3'UTR APA and DGE events. In this case study for colon cancer we can see that genes that undergo 3'UTR APA shortening events are being upregulated and the genes that undergo 3'UTR APA lengthening events are being downregulated.

## IPA CORRELATION ANALYSIS interface
<img src="https://github.com/brss12/APAtizer/assets/121204829/dc65aaae-8635-4438-a228-0af8bf70a613" alt="IPA_CORR">

Now, in this section we have the pearson correlation analysis scatter plot between IPA and DGE events. For colon cancer, we can also see that we have a significant negative correlation between IPA and DGE events. Genes that undergo IPA downregulation are being more expressed, whereas genes that undergo IPA upregulation are being less expressed.

# Final remarks
In this README, a case study for colon cancer (COAD) from the TCGA repository was used to demonstrate and explain the features and capabilities of the APAtizer tool. With this tool, the user can analyze RNA-Seq data from TCGA and retrieve many plots and useful information regarding 3'UTR APA & IPA events via DaPars2 and APAlyzer analysis, the function of those genes via gene ontology analysis, the common genes between cancers using venn diagram intersections and the correlation between 3'UTR APA & IPA events and DGE via pearson correlation analysis scatter plots.
