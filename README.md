# Description
APAtizer is a tool designed to analyse alternative polyadenylation and differential gene expression from TCGA RNA-Seq data. Additionally, it is capable of performing gene ontology analysis, visualizing Venn diagram intersections and survival analysis. APAtizer is equipped to handle TCGA files, such as BAM files, sample sheets, clinical data from patients and downstream analysis files such as htseq and DaPars txt files. It is a user-friendly interface, that allows users with limited knowledge on bioinformatics to generate informative visualizations, including Volcano plots, heatmaps and gene lists. The APAtizer tool also provides the functionality to download the aforementioned plots and gene lists for further analysis and exploration. 

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

## DaPars2 interface
### 3'UTR APA lengthening genes
<img src="https://github.com/brss12/APAtizer/assets/121204829/4d32cb27-7bfb-4add-894a-a64b9b2a148d" alt="DAPARS_LEN">

### 3'UTR APA shortening genes
<img src="https://github.com/brss12/APAtizer/assets/121204829/7c57fff7-d9f0-4f80-afe6-cb03c3dd9081" alt="DAPARS_SHORT">

In this section, in the input space the user can select the 24 output files originated by the DaPars2 analysis that are located in the folder **DaPars_data_final** and the TCGA sample sheet.

In the output space, we can observe the lists of genes that go through 3'UTR APA lengthening events (*Len genes*) and 3'UTR APA shortening events (*Short genes*). The user can also search the lists for a specific gene of interest and download the lists using the download button below the search box.

## APAlyzer interface
<img src="https://github.com/brss12/APAtizer/assets/121204829/f96b5fde-b5c6-4827-9a4d-542281aeca4f" alt="APA_List">
<img src="https://github.com/brss12/APAtizer/assets/121204829/a9ed2b83-dc63-462a-9bd2-43a32168f60f" alt="APA_Volcano">

Here, in the input space the user may paste the full path of the folder **TRIMMED_READS** that contains all the de-duplicated BAM files and the index files and select the TCGA sample sheet. The user may also select the output types of the analysis such as the lists that are displayed and the plots. For the lists, the user can choose between 3'UTR APA lengthening (*NvsT_APA_UP*), 3'UTR APA shortening (*NvsT_APA_DN*) and non-significant (*NvsT_APA_NC*). In the case of the plots, the choice is between a Volcano plot with the top 40 significant genes highlighted (*APA Volcano Plot (top40)*), the same plot but with no highlights (*APA Volcano Plot*) and a box plot (APA Box).

In the output space, in the tab called *Number of APA events* one can see a small table where the number of non significant, lengthening and shortening genes is present. In *NvsT_APA* the full lists chosen in the input space are presented to the user and he can search the list and download it. Finally, in *Plots* the user can visualize the plots selected in the input space and download them as well.

## DGE interface
<img src="https://github.com/brss12/APAtizer/assets/121204829/d69edb1f-7d95-4e78-9da4-9fecd09b7c44" alt="DGE_List">
<img src="https://github.com/brss12/APAtizer/assets/121204829/3ae16b84-87f0-494d-9c08-31b9c048f346" alt="DGE_Pca">
<img src="https://github.com/brss12/APAtizer/assets/121204829/2dd18423-06f7-4285-a528-8dc7b24b1ce2" alt="DGE_Volcano">
<img src="https://github.com/brss12/APAtizer/assets/121204829/4f76e1e6-ff5a-4c8c-8688-f3e9c1e386f9" alt="DGE_Heatmap">

For the differential gene expression analysis, in the input space the user may paste the full path for the folder **FILTERED** that has all of the htseq files and select the TCGA sample sheet. Once again, the user can also select the lists and the plots that will be displayed in the outputs. For the lists, the user can choose between DGE upregulated (*DGE_Genes_UP*), DGE downregulated (*DGE_Genes_DN*) and non-significant (*DGE_Genes_NC*). For the plots, the user can select a PCA plot (*PCA Plot*), a Volcano plot (*DGE Volcano Plot*) and a heatmap to evaluate the pattern of gene expression between conditions(*DGE Heatmap*).

In the output space, in the tab named *Number of DGE genes* the user can see a table showing the total number of upregulated, downregulated and non-significant genes. In *DGE_Genes* is where the full lists chosen before in the input space are displayed and the user can search those lists for a gene of interest and download them. In *Plots*, the same as in APAlyzer, is where one can visualize the selected plots and download them.

## GO TERMS interface
<img src="https://github.com/brss12/APAtizer/assets/121204829/9c48d65a-e3d9-4c5d-8b4f-c4d9a6f2fdee" alt="GO_BP">
<img src="https://github.com/brss12/APAtizer/assets/121204829/45e43e24-c6c1-48a1-8244-ffb25a837756" alt="GO_MF">

In this section, the user only needs to select the list of genes in which he wants to perform gene ontology exploration and the type of analysis to be performed, *Biological Process (BP)* and *Molecular Function (MF)*. The provided list should be one of the gene lists obtained in the previous steps (DaPars, APAlyzer or DGE analysis).

The output space only has one tab called *GO Plots* in which the resulting plots are displayed and can be downloaded.

## VENN DIAGRAMS interface
<img src="https://github.com/brss12/APAtizer/assets/121204829/0faebdf0-1baa-4be4-9d0b-75aab8cf99d3" alt="Venn_2">

For the Venn diagram intersections, the user can provide up to 5 gene lists obtained in the previous steps to execute the analysis. In the output section, the Venn diagram is displayed and can be downloaded in the tab *Venn Diagram*. Next, the user can obtain and download a list of the common genes between all the gene lists provided in the intersection.

## SURVIVAL ANALYSIS interface
<img src="https://github.com/brss12/APAtizer/assets/121204829/e661e322-0d8a-4914-8222-404d7d7f7e30" alt="Survival">

For the final section, a TCGA sample sheet should be provided along with the correspondent TCGA clinical data. Also, the full path for the folder **FILTERED** that contains all the htseq files must be provided along with a specific gene of interest to run the survival analysis. In the output side, in the tab *Plot* is where the resulting plot is displayed and can be downloaded.
