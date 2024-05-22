# Description
APAtizer is a tool designed to analyse alternative polyadenylation of TCGA RNA-Seq data. Additionally, it is capable of performing differential gene expression, gene ontology analysis, visualizing Venn diagram intersections and survival analysis. APAtizer is equipped to handle TCGA files, such as BAM files, sample sheets, clinical data from patients and downstream analysis files such as htseq and DaPars txt files. It is a user-friendly interface, that allows users with limited knowledge on bioinformatics to generate informative visualizations, including Volcano plots, heatmaps and gene lists. The APAtizer tool also provides the functionality to download the aforementioned plots and gene lists for further analysis and exploration. 

# Workflow
![image](https://github.com/brss12/APAtizer/assets/121204829/ea7c789d-907a-42bc-b331-94387a9b4325)

# Installing dependencies
To install the required dependencies for the creation of the input files for APAtizer, the user must run the [install_dependencies.sh](install_dependencies.sh) script.
```shell
./install_dependencies.sh
```

# Creating the input files
To start, clone the repository in the same directory where your raw BAM files are located. The script [APAtizer_inputs.smk](APAtizer_inputs.smk) employs the snakemake workflow to create fastqc reports, sort and remove the duplicates from the raw TCGA BAM files required for the APA analysis step using the APAlyzer algorithm and creates the Dapars files required for APA analysis employing the DaPars2 algorithm. In addition, it also creates the htseq files required for differential gene expression analysis using the DESeq2 package.

The script requires the raw BAM files to be placed in a folder called **RAW_BAM** and, upon running, creates two new folders called **SORTED_BAM** and **TRIMMED_READS**. In the **SORTED_BAM** folder, is where the sorted bam files along with fastqc reports are placed and in the **TRIMMED_READS** is where the de-duplicated bam files along with fastqc reports are placed. After this, a folder called **DaPars_data** is created with the final DaPars txt files.

Depending on the size and ammount of BAM files, we recommend running the script in a HPC environment. Run the script on the same directory of the **RAW_BAM** folder with the following command. 
```shell
nohup bash -c "snakemake -s APAtizer_inputs.smk --cores 1 --jobs 1 --wait-for-files" &
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

## VENN DIAGRAMS interface   TO DO
<img src="https://github.com/brss12/APAtizer/assets/121204829/0faebdf0-1baa-4be4-9d0b-75aab8cf99d3" alt="Venn_2">

For the Venn diagram intersections, the user can provide up to 5 gene lists obtained in the previous steps to execute the analysis. In the output section, the Venn diagram is displayed and can be downloaded in the tab *Venn Diagram*. Next, the user can obtain and download a list of the common genes between all the gene lists provided in the intersection.

## SURVIVAL ANALYSIS interface TO DO
<img src="https://github.com/brss12/APAtizer/assets/121204829/e661e322-0d8a-4914-8222-404d7d7f7e30" alt="Survival">

For the final section, a TCGA sample sheet should be provided along with the correspondent TCGA clinical data. Also, the full path for the folder **FILTERED** that contains all the htseq files must be provided along with a specific gene of interest to run the survival analysis. In the output side, in the tab *Plot* is where the resulting plot is displayed and can be downloaded.
