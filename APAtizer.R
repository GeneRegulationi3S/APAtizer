if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages_R <- c(
  "tidyverse", "base", "data.table", "stats", "shinyalert", "shiny", "purrr",
  "dplyr", "splitstackshape", "shinythemes", "readr", "ggplot2", "repmis",
  "pheatmap", "RColorBrewer", "EnhancedVolcano", "ggvenn", "ggpubr", "DT"
)

required_packages_Bioc <- c(
  "TCGAbiolinks", "SummarizedExperiment", "DESeq2", "APAlyzer", "Rsamtools", "apeglm",
  "clusterProfiler", "enrichplot", "org.Hs.eg.db", "org.Mm.eg.db", "VennDiagram", "EnhancedVolcano"
)

install_missing_packages_R <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      install.packages(package)
    }
  }
}

install_missing_packages_Bioc <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      BiocManager::install(package)
    }
  }
}
install_missing_packages_R(required_packages_R)
install_missing_packages_Bioc(required_packages_Bioc)

library("apeglm")
library("TCGAbiolinks")
library("SummarizedExperiment")
library("DESeq2")
library("APAlyzer")
library("data.table")
library("Rsamtools")
library("tidyverse")
library("stats")
library("shinyalert")
library("shiny")
library("DT")
library("purrr")
library("base")
library("dplyr")
library("splitstackshape")
library("shinythemes")
library("readr")
library("ggplot2")
library("repmis")
library("pheatmap")
library("RColorBrewer")
library("EnhancedVolcano")
library("clusterProfiler")
library("enrichplot")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("VennDiagram")
library("ggvenn")
library("ggpubr")

options(shiny.maxRequestSize=10000000*1024^2)

ui <- fluidPage(theme = shinytheme("darkly"),
                navbarPage(
                  "APAtizer",
                  tabPanel("SAMPLE SHEET",
                           titlePanel("Create a sample sheet with BAM file names and conditions"),
                           sidebarLayout(
                             sidebarPanel(
                               actionButton("addRow", "Add row"),
                               actionButton("removeRow", "Remove row"),
                               br(),
                               br(),
                               actionButton("saveData", "Save data")
                             ),
                             mainPanel(
                               tabsetPanel(
                                 tabPanel("Edit Table", DTOutput("editableTable")),
                                 tabPanel("View Data", tableOutput("displayData"))
                               )
                             )
                           ),
                           tags$style(HTML("
                                .dataTables_wrapper .dataTables_scroll {
                                  background-color: white;
                                  border: none !important;
                                }
                                table.dataTable tbody tr {
                                  background-color: white;
                                  border: none !important;
                                }
                                table.dataTable {
                                  background-color: white;
                                  border: none !important;
                                }
                                table.dataTable thead th {
                                  background-color: #f8f9fa;
                                  color: #333;
                                  border: none !important;
                                }
                                table.dataTable tbody td {
                                  color: #333;
                                  border: none !important;
                                }
                                .tab-content {
                                  border: 1px solid #ddd;
                                  padding: 10px;
                                  border: none !important;
                                }"))
                  ),
                  tabPanel("DAPARS",
                           titlePanel("Use DaPars2 to analyse de novo 3'UTR-APA from RNA-Seq data"),
                           sidebarLayout(
                             sidebarPanel(
                               fileInput("txt_files", "Select Multiple .txt Files", multiple = TRUE),
                               #fileInput("txt_file2", "Select sample sheet"),
                               actionButton("run2", "DaPars Analysis")
                             ),
                             mainPanel(
                               tabsetPanel(
                                 tabPanel("Len genes",
                                          br(),
                                          textInput(inputId = "search_term", label = "Search for a gene"),
                                          downloadButton("download_datax", "Download Len genes"),
                                          tableOutput("data_table_x")
                                 ),
                                 tabPanel("Short genes",
                                          br(),
                                          textInput(inputId = "search_term2", label = "Search for a gene"),
                                          downloadButton("download_datay", "Download Short genes"),
                                          tableOutput("data_table_y")
                                 )
                               )
                             )
                           )
                           
                  ),
                  tabPanel("APA APALYZER",
                           titlePanel("Use Apalyzer to analyse 3'UTR-APA from RNA-Seq data"),
                           sidebarLayout(
                             sidebarPanel(
                               textInput("path2", "TRIMMED BAM files directory path:"),
                               #fileInput("txt_file3", "Select sample sheet"),
                               selectInput(inputId = "ref_PAS1", label = "Select reference PAS",
                                           choices = c("hg19", "hg38", "mm9", "mm10")),
                               selectInput(inputId = "seq_method1", label = "Select sequencing method",
                                           choices = c("Paired-end", "Single-end")),
                               selectInput(inputId = "seq_strand1", label = "Select strandedness",
                                           choices = c("Forward stranded", "Reverse stranded", "Non-stranded")),
                               selectInput(inputId = "multi_test", label = "Select statistical test",
                                           choices = c("Unpaired t-test", "Paired t-test", "ANOVA")),
                               actionButton("run3", "APA Analysis"),
                               br(),
                               br(),
                               br(),
                               selectInput(inputId = "output_type4", label = "Select Output Type",
                                           choices = c("NvsT_APA_UP", "NvsT_APA_DN","NvsT_APA_NC" )),
                               br(),
                               selectInput(inputId = "output_type5", label = "Select Plot Type",
                                           choices = c("APA Volcano plot (top 40)", "APA Volcano plot", "APA Box"))
                             ),
                             mainPanel(
                               tabsetPanel(
                                 tabPanel("Number of APA events",
                                          tableOutput("data_table_2")
                                 ),
                                 tabPanel("NvsT_APA",
                                          conditionalPanel(
                                            condition = "input.output_type4 == 'NvsT_APA_UP'",
                                            br(),
                                            textInput(inputId = "search_term3", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data2", label = "Download NvsT_APA_UP"),
                                            tableOutput(outputId = "data_table_3")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type4 == 'NvsT_APA_DN'",
                                            br(),
                                            textInput(inputId = "search_term4", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data3", label = "Download NvsT_APA_DN"),
                                            tableOutput(outputId = "data_table_4")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type4 == 'NvsT_APA_NC'",
                                            br(),
                                            textInput(inputId = "search_term5", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data4", label = "Download NvsT_APA_NC"),
                                            tableOutput(outputId = "data_table_5")
                                          )
                                 ),
                                 tabPanel("Plots",
                                          conditionalPanel(
                                            condition = "input.output_type5 == 'APA Volcano plot (top 40)'",
                                            br(),
                                            fluidRow(
                                              column(3,
                                                     downloadButton(outputId = "download_plot1", label = "Download APA Volcano plot (top 40)")
                                              ),
                                              column(1,
                                                     div(style = "margin-top: 5px;", # Adjusts vertical alignment
                                                         selectInput("file_format", NULL, 
                                                                     choices = c("PNG", "PDF"), 
                                                                     selected = "PNG", 
                                                                     width = "100%")
                                                     )
                                              )
                                            ),
                                            plotOutput(outputId = "plot1")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type5 == 'APA Volcano plot'",
                                            br(),
                                            fluidRow(
                                              column(3,
                                                     downloadButton(outputId = "download_plot2", label = "Download APA Volcano plot")
                                              ),
                                              column(1,
                                                     div(style = "margin-top: 5px;", # Adjusts vertical alignment
                                                         selectInput("file_format", NULL, 
                                                                     choices = c("PNG", "PDF"), 
                                                                     selected = "PNG", 
                                                                     width = "100%")
                                                     )
                                              )
                                            ),
                                            plotOutput(outputId = "plot2")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type5 == 'APA Box'",
                                            br(),
                                            fluidRow(
                                              column(3,
                                                     downloadButton(outputId = "download_plot3", label = "Download APA Box")
                                              ),
                                              column(1,
                                                     div(style = "margin-top: 5px;", # Adjusts vertical alignment
                                                         selectInput("file_format", NULL, 
                                                                     choices = c("PNG", "PDF"), 
                                                                     selected = "PNG", 
                                                                     width = "100%")
                                                     )
                                              )
                                            ),
                                            plotOutput(outputId = "plot3")
                                          )
                                 )
                               )
                             )
                           )
                  ),
                  tabPanel("IPA APALYZER",
                           titlePanel("Use Apalyzer to analyse IPA from RNA-Seq data"),
                           sidebarLayout(
                             sidebarPanel(
                               textInput("path", "TRIMMED BAM files directory path:"),
                               #fileInput("txt_file", "Select sample sheet"),
                               selectInput(inputId = "ref_PAS2", label = "Select reference PAS",
                                           choices = c("hg19", "hg38", "mm9", "mm10")),
                               selectInput(inputId = "seq_method2", label = "Select sequencing method",
                                           choices = c("Paired-end", "Single-end")),
                               selectInput(inputId = "seq_strand2", label = "Select strandedness",
                                           choices = c("Forward stranded", "Reverse stranded", "Non-stranded")),
                               selectInput(inputId = "multi_test2", label = "Select statistical test",
                                           choices = c("Unpaired t-test", "Paired t-test", "ANOVA")),
                               sliderInput("n_threads", "Number of threads:", min = 1, max = 16, value = 1, step = 1),
                               actionButton("run", "IPA Analysis"),
                               br(),
                               br(),
                               br(),
                               selectInput(inputId = "output_type", label = "Select Output Type",
                                           choices = c("NvsT_IPA_events_UP", "NvsT_IPA_events_DN", "NvsT_IPA_events_NC")),
                               br(),
                               selectInput(inputId = "output_type2", label = "Select Output Type",
                                           choices = c("NvsT_IPA_genes_UP", "NvsT_IPA_genes_DN", "NvsT_IPA_genes_NC")),
                               br(),
                               selectInput(inputId = "output_type3", label = "Select Plot Type",
                                           choices = c("IPA Volcano plot (top 40)", "IPA Volcano plot", "IPA Box"))
                             ),
                             mainPanel(
                               tabsetPanel(
                                 tabPanel("Number of IPA events",
                                          tableOutput("data_table_6")
                                 ),
                                 tabPanel("NvsT_IPA_events",
                                          conditionalPanel(
                                            condition = "input.output_type == 'NvsT_IPA_events_UP'",
                                            br(),
                                            textInput(inputId = "search_term6", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data5", label = "Download NvsT_IPA_events_UP"),
                                            tableOutput(outputId = "data_table_7")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type == 'NvsT_IPA_events_DN'",
                                            br(),
                                            textInput(inputId = "search_term7", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data6", label = "Download NvsT_IPA_events_DN"),
                                            tableOutput(outputId = "data_table_8")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type == 'NvsT_IPA_events_NC'",
                                            br(),
                                            textInput(inputId = "search_term8", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data7", label = "Download NvsT_IPA_events_NC"),
                                            tableOutput(outputId = "data_table_9")
                                          )
                                 ),
                                 tabPanel("NvsT_IPA_genes",
                                          conditionalPanel(
                                            condition = "input.output_type2 == 'NvsT_IPA_genes_UP'",
                                            br(),
                                            textInput(inputId = "search_term9", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data8", label = "Download NvsT_IPA_genes_UP"),
                                            tableOutput(outputId = "data_table_10")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type2 == 'NvsT_IPA_genes_DN'",
                                            br(),
                                            textInput(inputId = "search_term10", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data9", label = "Download NvsT_IPA_genes_DN"),
                                            tableOutput(outputId = "data_table_11")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type2 == 'NvsT_IPA_genes_NC'",
                                            br(),
                                            textInput(inputId = "search_term11", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data10", label = "Download NvsT_IPA_genes_NC"),
                                            tableOutput(outputId = "data_table_12")
                                          )
                                 ),
                                 tabPanel("Plots",
                                          conditionalPanel(
                                            condition = "input.output_type3 == 'IPA Volcano plot (top 40)'",
                                            br(),
                                            fluidRow(
                                              column(3,
                                                     downloadButton(outputId = "download_plot4", label = "Download IPA Volcano plot (top 40)")
                                              ),
                                              column(1,
                                                     div(style = "margin-top: 5px;", # Adjusts vertical alignment
                                                         selectInput("file_format", NULL, 
                                                                     choices = c("PNG", "PDF"), 
                                                                     selected = "PNG", 
                                                                     width = "100%")
                                                     )
                                              )
                                            ),
                                            plotOutput(outputId = "plot4")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type3 == 'IPA Volcano plot'",
                                            br(),
                                            fluidRow(
                                              column(3,
                                                     downloadButton(outputId = "download_plot5", label = "Download IPA Volcano plot")
                                              ),
                                              column(1,
                                                     div(style = "margin-top: 5px;", # Adjusts vertical alignment
                                                         selectInput("file_format", NULL, 
                                                                     choices = c("PNG", "PDF"), 
                                                                     selected = "PNG", 
                                                                     width = "100%")
                                                     )
                                              )
                                            ),
                                            plotOutput(outputId = "plot5")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type3 == 'IPA Box'",
                                            br(),
                                            fluidRow(
                                              column(3,
                                                     downloadButton(outputId = "download_plot6", label = "Download IPA Box")
                                              ),
                                              column(1,
                                                     div(style = "margin-top: 5px;", # Adjusts vertical alignment
                                                         selectInput("file_format", NULL, 
                                                                     choices = c("PNG", "PDF"), 
                                                                     selected = "PNG", 
                                                                     width = "100%")
                                                     )
                                              )
                                            ),
                                            plotOutput(outputId = "plot6")
                                          )
                                 )
                               )
                             )
                           )
                  ),
                  tabPanel("DGE",
                           titlePanel("Use DESeq2 to analyse differentially expressed genes"),
                           sidebarLayout(
                             sidebarPanel(
                               textInput("path_dge", "htseq files directory path:"),
                               #fileInput("txt_file_dge", "Select sample sheet"),
                               actionButton("run_dge", "DGE Analysis"),
                               br(),
                               br(),
                               br(),
                               selectInput(inputId = "output_type_dge", label = "Select Output Type",
                                           choices = c("DGE_Genes_UP", "DGE_Genes_DN", "DGE_Genes_NC")),
                               br(),
                               selectInput(inputId = "output_type2_dge", label = "Select Output Type",
                                           choices = c("PCA Plot", "DGE Volcano Plot", "DGE Heatmap")),
                               # ConditionalPanel to show sliders only when "DGE Heatmap" is selected
                               conditionalPanel(
                                 condition = "input.output_type2_dge == 'DGE Heatmap'",
                                 sliderInput("cellwidth", "Heatmap Width:", min = 5, max = 100, value = 5, step = 1),
                                 sliderInput("cellheight", "Heatmap Height:", min = 0.01, max = 10, value = 0.01, step = 0.01)
                               )
                             ),
                             mainPanel(
                               tabsetPanel(
                                 tabPanel("Number of DGE genes",
                                          tableOutput("data_table_13")
                                 ),
                                 tabPanel("DGE_Genes",
                                          conditionalPanel(
                                            condition = "input.output_type_dge == 'DGE_Genes_UP'",
                                            br(),
                                            textInput(inputId = "search_term_dge", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data_dge", label = "Download DGE_Genes_UP"),
                                            tableOutput(outputId = "data_table_14")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type_dge == 'DGE_Genes_DN'",
                                            br(),
                                            textInput(inputId = "search_term2_dge", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data2_dge", label = "Download DGE_Genes_DN"),
                                            tableOutput(outputId = "data_table_15")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type_dge == 'DGE_Genes_NC'",
                                            br(),
                                            textInput(inputId = "search_term3_dge", label = "Search for a gene"),
                                            downloadButton(outputId = "download_data3_dge", label = "Download DGE_Genes_NC"),
                                            tableOutput(outputId = "data_table_16")
                                          )
                                 ),
                                 tabPanel("Plots",
                                          conditionalPanel(
                                            condition = "input.output_type2_dge == 'PCA Plot'",
                                            br(),
                                            fluidRow(
                                              column(3,
                                                     downloadButton(outputId = "download_plot7", label = "Download PCA Plot")
                                              ),
                                              column(1,
                                                     div(style = "margin-top: 5px;", # Adjusts vertical alignment
                                                         selectInput("file_format", NULL, 
                                                                     choices = c("PNG", "PDF"), 
                                                                     selected = "PNG", 
                                                                     width = "100%")
                                                     )
                                              )
                                            ),
                                            plotOutput(outputId = "plot7")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type2_dge == 'DGE Volcano Plot'",
                                            br(),
                                            fluidRow(
                                              column(3,
                                                     downloadButton(outputId = "download_plot8", label = "Download DGE Volcano plot")
                                              ),
                                              column(1,
                                                     div(style = "margin-top: 5px;", # Adjusts vertical alignment
                                                         selectInput("file_format", NULL, 
                                                                     choices = c("PNG", "PDF"), 
                                                                     selected = "PNG", 
                                                                     width = "100%")
                                                     )
                                              )
                                            ),
                                            plotOutput(outputId = "plot8")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type2_dge == 'DGE Heatmap'",
                                            br(),
                                            fluidRow(
                                              column(3,
                                                     downloadButton(outputId = "download_plot9", label = "Download DGE Heatmap")
                                              ),
                                              column(1,
                                                     div(style = "margin-top: 5px;", # Adjusts vertical alignment
                                                         selectInput("file_format", NULL, 
                                                                     choices = c("PNG", "PDF"), 
                                                                     selected = "PNG", 
                                                                     width = "100%")
                                                     )
                                              )
                                            ),
                                            plotOutput(outputId = "plot9")
                                          )
                                 )
                               )
                             )
                           )
                  ),
                  tabPanel("GO TERMS",
                           titlePanel("Perform Gene Ontology analysis"),
                           sidebarLayout(
                             sidebarPanel(
                               fileInput("txt_file_go", "Select Gene List"),
                               selectInput(inputId = "database_go", label = "Select Organism Database",
                                           choices = c("Human", "Mouse")),
                               actionButton("run_go", "GO Analysis"),
                               br(),
                               br(),
                               br(),
                               selectInput(inputId = "output_type_go", label = "Select Output Type",
                                           choices = c("Biological Process (BP)", "Molecular Function (MF)"))
                             ),
                             mainPanel(
                               tabsetPanel(
                                 tabPanel("GO Plots",
                                          conditionalPanel(
                                            condition = "input.output_type_go == 'Biological Process (BP)'",
                                            br(),
                                            fluidRow(
                                              column(3,
                                                     downloadButton(outputId = "download_plot_go", label = "Download GO Plot BP")
                                              ),
                                              column(1,
                                                     div(style = "margin-top: 5px;", # Adjusts vertical alignment
                                                         selectInput("file_format", NULL, 
                                                                     choices = c("PNG", "PDF"), 
                                                                     selected = "PNG", 
                                                                     width = "100%")
                                                     )
                                              )
                                            ),
                                            plotOutput(outputId = "plot_go")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type_go == 'Molecular Function (MF)'",
                                            br(),
                                            fluidRow(
                                              column(3,
                                                     downloadButton(outputId = "download_plot2_go", label = "Download APA GO Plot MF")
                                              ),
                                              column(1,
                                                     div(style = "margin-top: 5px;", # Adjusts vertical alignment
                                                         selectInput("file_format", NULL, 
                                                                     choices = c("PNG", "PDF"), 
                                                                     selected = "PNG", 
                                                                     width = "100%")
                                                     )
                                              )
                                            ),
                                            plotOutput(outputId = "plot2_go")
                                          )
                                 )
                               )
                             )
                           )
                  ),
                  tabPanel("VENN DIAGRAMS",
                           titlePanel("Perform Venn Diagram analysis"),
                           sidebarLayout(
                             sidebarPanel(
                               fileInput("txt_file_venn", "Select Gene Lists", multiple = TRUE),
                               actionButton("run_venn", "Venn Diagram Analysis")
                             ),
                             mainPanel(
                               tabsetPanel(
                                 tabPanel("Venn Diagram",
                                          br(),
                                          fluidRow(
                                            column(3,
                                                   downloadButton(outputId = "download_plot_venn", label = "Download Venn Diagram")
                                            ),
                                            column(1,
                                                   div(style = "margin-top: 5px;", # Adjusts vertical alignment
                                                       selectInput("file_format", NULL, 
                                                                   choices = c("PNG", "PDF"), 
                                                                   selected = "PNG", 
                                                                   width = "100%")
                                                   )
                                            )
                                          ),
                                          plotOutput(outputId = "plot_venn")
                                 ),
                                 tabPanel("Common Genes",
                                          br(),
                                          textInput(inputId = "search_term_common_genes", label = "Search for a gene"),
                                          downloadButton(outputId = "download_common_genes", label = "Download Common Genes"),
                                          tableOutput(outputId = "common_genes")
                                       )
                                    )
                                )
                            )
                  ),
                  tabPanel("APA CORRELATION ANALYSIS",
                           titlePanel("Perform pearson correlation analysis between APA and DGE"),
                           sidebarLayout(
                             sidebarPanel(
                               fileInput("corr_apa", "Select APA LEN and SHORT lists", multiple = TRUE),
                               fileInput("corr_dge", "Select DGE UP, NC and DN lists", multiple = TRUE),
                               br(),
                               actionButton("run_corr", "Correlation Analysis")
                             ),
                             mainPanel(
                               tabsetPanel(
                                 tabPanel("Plot",
                                          br(),
                                          fluidRow(
                                            column(3,
                                                   downloadButton(outputId = "download_plot_corr", label = "Download Correlation Analysis Plot")
                                            ),
                                            column(1,
                                                   div(style = "margin-top: 5px;", # Adjusts vertical alignment
                                                       selectInput("file_format", NULL, 
                                                                   choices = c("PNG", "PDF"), 
                                                                   selected = "PNG", 
                                                                   width = "100%")
                                                   )
                                            )
                                          ),
                                          plotOutput(outputId = "plot_corr")
                                        
                                 )
                               )
                             )
                           )
                  ),
                  tabPanel("IPA CORRELATION ANALYSIS",
                           titlePanel("Perform pearson correlation analysis between IPA and DGE"),
                           sidebarLayout(
                             sidebarPanel(
                               fileInput("corr_ipa", "Select IPA events UP and DN lists", multiple = TRUE),
                               fileInput("corr_dge2", "Select DGE UP, NC and DN lists", multiple = TRUE),
                               br(),
                               actionButton("run_corr2", "Correlation Analysis")
                             ),
                             mainPanel(
                               tabsetPanel(
                                 tabPanel("Plot",
                                          br(),
                                          fluidRow(
                                            column(3,
                                                   downloadButton(outputId = "download_plot_corr2", label = "Download Correlation Analysis Plot")
                                            ),
                                            column(1,
                                                   div(style = "margin-top: 5px;", # Adjusts vertical alignment
                                                       selectInput("file_format", NULL, 
                                                                   choices = c("PNG", "PDF"), 
                                                                   selected = "PNG", 
                                                                   width = "100%")
                                                   )
                                            )
                                          ),
                                          plotOutput(outputId = "plot_corr2")
                                 )
                               )
                             )
                           )
                  )
          )
)
                                            
                                            

server <- function(input, output,session) {
  
  ##### SAMPLE SHEET #####
  
  # Initialize reactive values with data from .GlobalEnv if available
  if (exists("savedData", envir = .GlobalEnv)) {
    initialData <- get("savedData", envir = .GlobalEnv)
  } else {
    initialData <- data.frame(File.Name = character(), Sample.Type = character(), stringsAsFactors = FALSE)
  }
  
  rv <- reactiveValues(data = initialData)
  
  # Render the editable data table
  output$editableTable <- renderDT({
    datatable(rv$data, editable = TRUE, rownames = FALSE,
              options = list(dom = 't', paging = FALSE, ordering = FALSE, searching = FALSE))
  }, server = FALSE)
  
  # Observe add row button
  observeEvent(input$addRow, {
    newRow <- data.frame(File.Name = "", Sample.Type = "", stringsAsFactors = FALSE)
    rv$data <- rbind(rv$data, newRow)
  })
  
  # Observe remove row button
  observeEvent(input$removeRow, {
    if (nrow(rv$data) > 0) {
      rv$data <- rv$data[-nrow(rv$data), ]  # Remove the last row
    }
  })
  
  # Observe changes to the table
  observeEvent(input$editableTable_cell_edit, {
    info <- input$editableTable_cell_edit
    rv$data[info$row, info$col + 1] <- info$value
  })
  
  # Save the data to the R environment
  observeEvent(input$saveData, {
    assign("savedData", rv$data, envir = .GlobalEnv)  # Save data to .GlobalEnv
    showModal(modalDialog(
      title = "Data Saved",
      "The data has been saved for use in the tool.",
      easyClose = TRUE
    ))
  })
  
  # Display the data in another tab
  output$displayData <- renderTable({
    rv$data
  })
  
  
  ##### IPA #####
  
  df_pacientes <- eventReactive(input$run,{
    #file <- input$txt_file
    datapath <- input$path
    #if (is.null(file) || datapath == "") {
    #  shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
    #  return(NULL)
    #}
    
    if (nrow(rv$data) == 0) {
      shinyalert("Error", "Please make sure that the sample sheet is correct before running the analysis.", type = "error")
      return(NULL)
    }
    
    df_pacientes <- rv$data
    
    #df_pacientes <- read.table(file$datapath, sep = "\t", header = TRUE)
    df_pacientes <- dplyr::select(df_pacientes, File.Name,Sample.Type)
    
    df_pacientes$Sample.Type <- gsub(" ", "_", df_pacientes$Sample.Type)
    #df_pacientes$File.Name<- str_replace(df_pacientes$File.Name, "rna_seq.genomic.gdc_realn.bam", "rna_seq.genomic.gdc_realn.trim.bam")
    
    #case_order <- unique(df_pacientes$Case.ID)
    #df_pacientes <- df_pacientes %>% arrange(df_pacientes$Sample.Type, match(df_pacientes$Case.ID, case_order))
    
    #normal <- c("Solid Tissue Normal")
    #df_pacientes$category <- ifelse(df_pacientes$Sample.Type %in% c("Solid Tissue Normal"), "Normal", "Tumor")
    #df_pacientes$category = paste(df_pacientes$Case.ID, df_pacientes$category, sep="_")
    
    return(df_pacientes)
    
    
  })
  
  
  NvsT_IPA <- eventReactive(input$run, {
    #file <- input$txt_file
    datapath <- input$path
    #if (is.null(file) || datapath == "") {
    #  shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
    #  return(NULL)
    #}
    
    shinyalert("Processing", "IPA Analysis has started. Please wait...", type = "info")
    
    datapath <- str_replace_all(datapath, "\\\\", "/")
    flsall <- dir(datapath, pattern = ".*\\.bam$")
    flsall <- paste0(datapath, flsall)
    names(flsall) <- gsub('.bam','',dir(datapath, pattern = ".*\\.bam$"))
    
    #flsall <- paste0(datapath, df_pacientes()$File.Name) #D:/TRIMMED_READS is where the bam files are
    #names(flsall) <- df_pacientes()$category
    #flsall
    
    #Genomic reference
    library("repmis")
    URL="https://github.com/RJWANGbioinfo/PAS_reference_RData/blob/master/"
    file=paste0(input$ref_PAS2, "_REF.RData")
    source_data(paste0(URL,file,"?raw=True"))
    
    if (input$ref_PAS2 == "hg19") {
      refUTRraw=refUTRraw_hg19
      dfIPAraw=dfIPA_hg19
      dfLEraw=dfLE_hg19
    } else if (input$ref_PAS2 == "hg38") {
      refUTRraw=refUTRraw_hg38
      dfIPAraw=dfIPA_hg38
      dfLEraw=dfLE_hg38
    } else if (input$ref_PAS2 == "mm9") {
      refUTRraw=refUTRraw
      dfIPAraw=dfIPA
      dfLEraw=dfLE
    } else if (input$ref_PAS2 == "mm10") {
      refUTRraw=refUTRraw
      dfIPAraw=dfIPA
      dfLEraw=dfLE
    }
    
    PASREF=REF4PAS(refUTRraw,dfIPAraw,dfLEraw)
    UTRdbraw=PASREF$UTRdbraw
    dfIPA=PASREF$dfIPA
    dfLE=PASREF$dfLE   
    dfIPA=dfIPA
    dfLE=dfLE
    
    #IPA_OUTraw <- read_csv("/home/bruno/I3S/COAD/Non-Smoking/TRIMMED_APALYZER/IPA_OUTraw.csv")
    #IPA_OUTraw=PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype="forward", SeqType="ThreeMostPairEnd")
    
    if (input$seq_strand2=="Forward stranded" & input$seq_method2=="Paired-end") {
      IPA_OUTraw=PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype="forward", SeqType="ThreeMostPairEnd", nts=input$n_threads)
    } else if (input$seq_strand2 == "Reverse stranded" & input$seq_method2=="Paired-end") {
      IPA_OUTraw=PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype="invert", SeqType="ThreeMostPairEnd", nts=input$n_threads)
    } else if (input$seq_strand2 == "Non-stranded" & input$seq_method2=="Paired-ends") {
      IPA_OUTraw=PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype="NONE", SeqType="ThreeMostPairEnd", nts=input$n_threads)
    } else if (input$seq_strand2 == "Forward stranded" & input$seq_method2=="Single-end") {
      IPA_OUTraw=PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype="forward", nts=input$n_threads)
    } else if (input$seq_strand2 == "Reverse stranded" & input$seq_method2=="Single-end") {
      IPA_OUTraw=PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype="invert", nts=input$n_threads)
    } else if (input$seq_strand2 == "Non-stranded" & input$seq_method2=="Single-end") {
      IPA_OUTraw=PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype="NONE", nts=input$n_threads)
    }
    
    x=nrow(df_pacientes())/2
    sampleTable2 = data.frame(samplename = gsub(".bam", "", df_pacientes()$File.Name),
                              condition = c(rep(unique(df_pacientes()$Sample.Type)[1],x),rep(unique(df_pacientes()$Sample.Type)[2],x)))
    sampleTable2$samplename <- paste0(sampleTable2$samplename, ".trim")
    
    #x=nrow(df_pacientes())/2
    #sampleTable2 = data.frame(samplename = c(names(flsall)),
    #                          condition = c(rep("KD",x),rep("NT",x))) #OTIMIZAR, se são 2 amostras, 1 pro KD e 1 pro NT
    
    if (input$multi_test2 == "Unpaired t-test") {
      NvsT_IPA=APAdiff(sampleTable2, IPA_OUTraw, 
                       conKET=unique(df_pacientes()$Sample.Type)[2],
                       trtKEY=unique(df_pacientes()$Sample.Type)[1],
                       PAS='IPA', CUTreads=5, p_adjust_methods="fdr", MultiTest = "unpaired t-test")
    } else if (input$multi_test2 == "Paired t-test") {
      NvsT_IPA=APAdiff(sampleTable2, IPA_OUTraw, 
                       conKET=unique(df_pacientes()$Sample.Type)[2],
                       trtKEY=unique(df_pacientes()$Sample.Type)[1],
                       PAS='IPA', CUTreads=5, p_adjust_methods="fdr", MultiTest = "paired t-test")
    } else if (input$multi_test2 == "ANOVA") {
      NvsT_IPA=APAdiff(sampleTable2, IPA_OUTraw, 
                       conKET=unique(df_pacientes()$Sample.Type)[2],
                       trtKEY=unique(df_pacientes()$Sample.Type)[1],
                       PAS='IPA', CUTreads=5, p_adjust_methods="fdr", MultiTest = "ANOVA")
    }
    
    NvsT_IPA <- as.data.frame(NvsT_IPA)
    
    shinyalert("Success", "IPA Analysis has been completed successfully.", type = "success")
    
    return(NvsT_IPA)
  })
  
  Nr_IPA_events <- eventReactive(input$run,{
    Nr_IPA_events<- table(NvsT_IPA()$APAreg)  
    return(Nr_IPA_events)
  })
  # NvsT_IPA_UP
  NvsT_IPA_events_UP <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    NvsT_IPA_events_UP <- NvsT_IPA[ which(NvsT_IPA$APAreg=='UP'),]
    
    return(NvsT_IPA_events_UP)
  })
  
  # NvsT_IPA_DN
  NvsT_IPA_events_DN <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    NvsT_IPA_events_DN <- NvsT_IPA[ which(NvsT_IPA$APAreg=='DN'),]
    
    return(NvsT_IPA_events_DN)
  })
  
  # NvsT_IPA_NC
  NvsT_IPA_events_NC <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    NvsT_IPA_events_NC <- NvsT_IPA[ which(NvsT_IPA$APAreg=='NC'),]
    
    return(NvsT_IPA_events_NC)
  })
  
  NvsT_IPA_genes_UP <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    NvsT_IPA_events_UP <- NvsT_IPA_events_UP()
    NvsT_IPA_genes_UP <- distinct(NvsT_IPA_events_UP,select=c(gene_symbol))
    NvsT_IPA_genes_UP <- NvsT_IPA_genes_UP %>% rename(gene_symbol = select)
    
    return(NvsT_IPA_genes_UP)
  })
  
  # NvsT_IPA_DN
  NvsT_IPA_genes_DN <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    NvsT_IPA_events_DN <- NvsT_IPA_events_DN()
    NvsT_IPA_genes_DN <- distinct(NvsT_IPA_events_DN,select=c(gene_symbol))
    NvsT_IPA_genes_DN <- NvsT_IPA_genes_DN %>% rename(gene_symbol = select)
    
    return(NvsT_IPA_genes_DN)
  })
  
  # NvsT_IPA_NC
  NvsT_IPA_genes_NC <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    NvsT_IPA_events_NC <- NvsT_IPA_events_NC()
    NvsT_IPA_genes_NC <- distinct(NvsT_IPA_events_NC,select=c(gene_symbol))
    NvsT_IPA_genes_NC <- NvsT_IPA_genes_NC %>% rename(gene_symbol = select)
    
    return(NvsT_IPA_genes_NC)
  })
  
  e <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    e <- APAVolcano(NvsT_IPA, PAS='IPA', Pcol = "pvalue", top=40)
    
    return(e)
  })
  f <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    f <- APAVolcano(NvsT_IPA, PAS='IPA', Pcol = "pvalue")
    
    return(f)
  })
  g <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    g <- APABox(NvsT_IPA, xlab = "APAreg", ylab = "RED", plot_title = NULL)
    
    return(g)
  })
  
  
  output$data_table_6 <- renderTable({
    Nr_IPA_events()
  })
  
  output$data_table_7 <- renderTable({
    if (input$search_term6 != "") {
      NvsT_IPA_events_UP() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term6, ignore_case = TRUE))))
    } else {
      NvsT_IPA_events_UP()
    }
  })
  
  output$data_table_8 <- renderTable({
    if (input$search_term7 != "") {
      NvsT_IPA_events_DN() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term7, ignore_case = TRUE))))
    } else {
      NvsT_IPA_events_DN()
    }
  })
  
  output$data_table_9 <- renderTable({
    if (input$search_term8 != "") {
      NvsT_IPA_events_NC() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term8, ignore_case = TRUE))))
    } else {
      NvsT_IPA_events_NC()
    }
  })
  
  output$data_table_10 <- renderTable({
    if (input$search_term9 != "") {
      NvsT_IPA_genes_UP() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term9, ignore_case = TRUE))))
    } else {
      NvsT_IPA_genes_UP()
    }
  })
  
  output$data_table_11 <- renderTable({
    if (input$search_term10 != "") {
      NvsT_IPA_genes_DN() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term10, ignore_case = TRUE))))
    } else {
      NvsT_IPA_genes_DN()
    }
  })
  
  output$data_table_12 <- renderTable({
    if (input$search_term11 != "") {
      NvsT_IPA_genes_NC() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term11, ignore_case = TRUE))))
    } else {
      NvsT_IPA_genes_NC()
    }
  })
  
  output$plot4 <- renderPlot({
    e()
  }, width = 1200, height = 750)
  output$plot5 <- renderPlot({
    f()
  }, width = 1200, height = 750)
  output$plot6 <- renderPlot({
    g()
  }, width = 1200, height = 750)
  
  
  
  output$download_data5 <- downloadHandler(
    filename = function() {
      paste("NvsT_IPA_events_UP", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_IPA_events_UP(), file, row.names = FALSE)
    }
  )
  output$download_data6 <- downloadHandler(
    filename = function() {
      paste("NvsT_IPA_events_DN", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_IPA_events_DN(), file, row.names = FALSE)
    }
  )
  output$download_data7 <- downloadHandler(
    filename = function() {
      paste("NvsT_IPA_events_NC", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_IPA_events_NC(), file, row.names = FALSE)
    }
  )
  output$download_data8 <- downloadHandler(
    filename = function() {
      paste("NvsT_IPA_genes_UP", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_IPA_genes_UP(), file, row.names = FALSE)
    }
  )
  output$download_data9 <- downloadHandler(
    filename = function() {
      paste("NvsT_IPA_genes_DN", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_IPA_genes_DN(), file, row.names = FALSE)
    }
  )
  output$download_data10 <- downloadHandler(
    filename = function() {
      paste("NvsT_IPA_genes_NC", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_IPA_genes_NC(), file, row.names = FALSE)
    }
  )
  
  output$download_plot4 <- downloadHandler(
    filename = function() {
      paste("plot4", switch(input$file_format, "PNG" = ".png", "PDF" = ".pdf"), sep = "")
    },
    content = function(file) {
      if (input$file_format == "PNG") {
        ggsave(file, plot = e(), width = 5000, height = 5000, units = "px", bg = "white", dpi = 600)
      } else if (input$file_format == "PDF") {
        ggsave(file, plot = e(), width = 7, height = 7, device = "pdf", bg = "white")
      }
    }
  )
  
  output$download_plot5 <- downloadHandler(
    filename = function() {
      paste("plot5", switch(input$file_format, "PNG" = ".png", "PDF" = ".pdf"), sep = "")
    },
    content = function(file) {
      if (input$file_format == "PNG") {
        ggsave(file, plot = f(), width = 5000, height = 5000, units = "px", bg = "white", dpi = 600)
      } else if (input$file_format == "PDF") {
        ggsave(file, plot = f(), width = 7, height = 7, device = "pdf", bg = "white")
      }
    }
  )
  
  output$download_plot6 <- downloadHandler(
    filename = function() {
      paste("plot6", switch(input$file_format, "PNG" = ".png", "PDF" = ".pdf"), sep = "")
    },
    content = function(file) {
      if (input$file_format == "PNG") {
        ggsave(file, plot = g(), width = 5000, height = 5000, units = "px", bg = "white", dpi = 600)
      } else if (input$file_format == "PDF") {
        ggsave(file, plot = g(), width = 7, height = 7, device = "pdf", bg = "white")
      }
    }
  )
  
  
  ##### DAPARS #####
  
  df_pacientes2 <- eventReactive(input$run2,{
    files <- input$txt_files
    #file <- input$txt_file2
    if (is.null(files)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    
    if (nrow(rv$data) == 0) {
      shinyalert("Error", "Please make sure that the sample sheet is correct before running the analysis.", type = "error")
      return(NULL)
    }
    
    df_pacientes2 <- rv$data
    
    #df_pacientes2 <- read.table(file$datapath, sep = "\t", header = TRUE)
    df_pacientes2 <- dplyr::select(df_pacientes2, File.Name,Sample.Type)
    
    #df_pacientes2 <- read.table(file$datapath, sep = "\t", header = TRUE)
    #df_pacientes2 <- dplyr::select(df_pacientes2, File.Name,Case.ID,Sample.Type)
    df_pacientes2$File.Name <- paste("../WIG/", df_pacientes2$File.Name, sep = "") #"./WIG/" or "WIG/"
    #case_order <- unique(df_pacientes2$Case.ID)
    #df_pacientes2 <- df_pacientes2 %>% arrange(df_pacientes2$Sample.Type, match(df_pacientes2$Case.ID, case_order))
    #normal <- c("Solid Tissue Normal")
    #df_pacientes2$category <- ifelse(df_pacientes2$Sample.Type %in% normal, "Normal", "Tumor")
    #df_pacientes2$category = paste(df_pacientes2$Case.ID, df_pacientes2$category, sep="_")
    df_pacientes2$File.Name<- str_replace(df_pacientes2$File.Name, ".bam", "_PDUI")
    df_pacientes2$Sample.Type <- gsub(" ", "_", df_pacientes2$Sample.Type)
    
    return(df_pacientes2)
    
    
  })
  
  # DPDUI
  dpdui <- eventReactive(input$run2,{
    files <- input$txt_files
    #file <- input$txt_file2
    if (is.null(files)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    
    df <- rbindlist(sapply(files$datapath, fread, simplify = FALSE), use.names = TRUE)
    df <- subset(df, select = -c(fit_value,Loci,Predicted_Proximal_APA))
    
    idx <- match(df_pacientes2()$File.Name, colnames(df))
    idx <- append(1, idx)
    
    df <- df[, ..idx]
    
    num_cols = ncol(df)
    
    colnames(df)[2:num_cols] <- df_pacientes2()$Sample.Type
    
    num_cols2 = ((ncol(df)-1)/2)+1
    num_cols3 = ((ncol(df)-1)/2)+2
    num_cols4 = ncol(df)+1
    
    for (i in colnames(df)[2:num_cols2]) {
      df[[i]][is.na(df[[i]])] <- rowMeans(df[,2:num_cols2], na.rm = TRUE)[is.na(df[[i]])]
    }
    
    for (i in colnames(df)[num_cols3:num_cols]) {
      df[[i]][is.na(df[[i]])] <- rowMeans(df[,num_cols3:num_cols], na.rm = TRUE)[is.na(df[[i]])]
    }
    
    df <- data.frame(df)
    res <- df[, grepl(unique(df_pacientes2()$Sample.Type)[1], colnames(df))] - df[, grepl(unique(df_pacientes2()$Sample.Type)[2], colnames(df))]
    
    colnames(res) <- paste(colnames(df[, grepl(unique(df_pacientes2()$Sample.Type)[1], colnames(df))]),
                           colnames(df[, grepl(unique(df_pacientes2()$Sample.Type)[2], colnames(df))]), sep = "-")
    
    df <- cbind(df, res)
    
    num_cols5 = ncol(df)
    
    dpdui <- df[,c(1,num_cols4:num_cols5)]

    dpdui$Gene <- sapply(strsplit(as.character(dpdui$Gene), "\\|"), `[`, 2)
    dpdui <- dpdui %>% rename(gene_symbol = Gene)

    return(dpdui)
  })
  
  # SHORT GENES
  short_genes <- eventReactive(input$run2,{
    dpdui <- dpdui()
    dpdui$Mean <- rowMeans(dpdui[,2:length(dpdui())],na.rm=TRUE) 
    dpdui$Len <- dpdui[,c("Mean")] >= 0.2
    dpdui$Short <- dpdui[,c("Mean")] <= -0.2
    gene_len <- dpdui[which(dpdui$Len == 1), ]                             
    gene_short <- dpdui[which(dpdui$Short == 1), ]
    
    return(gene_short)
  })
  
  # LEN GENES
  len_genes <- eventReactive(input$run2,{
    dpdui <- dpdui()
    dpdui$Mean <- rowMeans(dpdui[,2:length(dpdui())],na.rm=TRUE) 
    dpdui$Len <- dpdui[,c("Mean")] >= 0.2
    dpdui$Short <- dpdui[,c("Mean")] <= -0.2
    gene_len <- dpdui[which(dpdui$Len == 1), ]                             
    gene_short <- dpdui[which(dpdui$Short == 1), ]  
    
    return(gene_len)
  })
  

  # Display combined data as a table
  output$data_table_x <- renderTable({
    if (input$search_term != "") {
      len_genes() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term, ignore_case = TRUE))))
    } else {
      len_genes()
    }
  })
  
  output$data_table_y <- renderTable({
    if (input$search_term2 != "") {
      short_genes() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term2, ignore_case = TRUE))))
    } else {
      short_genes()
    }
  })
  
  
  # Download combined data as a .csv file

  output$download_datax <- downloadHandler(
    filename = function() {
      paste("Len_Genes", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(len_genes(), file, row.names = FALSE)
    }
  )
  output$download_datay <- downloadHandler(
    filename = function() {
      paste("Short_Genes", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(short_genes(), file, row.names = FALSE)
    }
  )
  
  
  ##### APA #####
  
  df_pacientes3 <- eventReactive(input$run3,{
    #file <- input$txt_file3
    datapath <- input$path2
    
    if (datapath == "") {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    
    if (nrow(rv$data) == 0) {
      shinyalert("Error", "Please make sure that the sample sheet is correct before running the analysis.", type = "error")
      return(NULL)
    }
    
    df_pacientes3 <- rv$data
    
    #df_pacientes3 <- read.table(file$datapath, sep = "\t", header = TRUE)
    df_pacientes3 <- dplyr::select(df_pacientes3, File.Name,Sample.Type)
    
    df_pacientes3$Sample.Type <- gsub(" ", "_", df_pacientes3$Sample.Type)
    #df_pacientes3$File.Name<- str_replace(df_pacientes3$File.Name, "rna_seq.genomic.gdc_realn.bam", "rna_seq.genomic.gdc_realn.trim.bam")
    
    #case_order <- unique(df_pacientes3$Case.ID)
    #df_pacientes3 <- df_pacientes3 %>% arrange(df_pacientes3$Sample.Type, match(df_pacientes3$Case.ID, case_order))
    
    #normal <- c("Solid Tissue Normal")
    #df_pacientes3$category <- ifelse(df_pacientes3$Sample.Type %in% c("Solid Tissue Normal"), "Normal", "Tumor")
    #df_pacientes3$category = paste(df_pacientes3$Case.ID, df_pacientes3$category, sep="_")
    
    return(df_pacientes3)
    
    
  })
  
  NvsT_APA <- eventReactive(input$run3,{
    #file <- input$txt_file3
    datapath <- input$path2
    
    if (datapath == "") {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    
    shinyalert("Processing", "APA Analysis has started. Please wait...", type = "info",
               closeOnEsc = FALSE, closeOnClickOutside = FALSE, showConfirmButton = TRUE)
    
    datapath <- str_replace_all(datapath, "\\\\", "/")
    flsall <- dir(datapath, pattern = ".*\\.bam$")
    flsall <- paste0(datapath, flsall)
    names(flsall) <- gsub('.bam','',dir(datapath, pattern = ".*\\.bam$"))
    
    #flsall <- paste0(datapath, df_pacientes3()$File.Name) #D:/TRIMMED_READS is where the bam files are
    #names(flsall) <- df_pacientes3()$category
    #flsall
    
    #Genomic reference
    URL="https://github.com/RJWANGbioinfo/PAS_reference_RData/blob/master/"
    file=paste0(input$ref_PAS1, "_REF.RData")
    source_data(paste0(URL,file,"?raw=True"))
    
    if (input$ref_PAS1 == "hg19") {
      refUTRraw=refUTRraw_hg19
      dfIPAraw=dfIPA_hg19
      dfLEraw=dfLE_hg19
    } else if (input$ref_PAS1 == "hg38") {
      refUTRraw=refUTRraw_hg38
      dfIPAraw=dfIPA_hg38
      dfLEraw=dfLE_hg38
    } else if (input$ref_PAS1 == "mm9") {
      refUTRraw=refUTRraw
      dfIPAraw=dfIPA
      dfLEraw=dfLE
    } else if (input$ref_PAS1 == "mm10") {
      refUTRraw=refUTRraw
      dfIPAraw=dfIPA
      dfLEraw=dfLE
    }
    
    refUTRraw=refUTRraw
    dfIPAraw=dfIPAraw
    dfLEraw=dfLEraw
    PASREF=REF4PAS(refUTRraw,dfIPAraw,dfLEraw)
    UTRdbraw=PASREF$UTRdbraw
    dfIPA=PASREF$dfIPA
    dfLE=PASREF$dfLE
    
    #Analysis of APA in 3’UTRs
    refUTRraw=refUTRraw
    UTRdbraw=REF3UTR(refUTRraw)
    #DFUTRraw <- read_csv("/home/bruno/I3S/IPO/Lexogen_IPO_no6/No_6_APALYZER/DFUTRraw.csv")
    
    if (input$seq_strand1 == "Forward stranded") {
      DFUTRraw=PASEXP_3UTR(UTRdbraw, flsall, Strandtype="forward")
    } else if (input$seq_strand1 == "Reverse stranded") {
      DFUTRraw=PASEXP_3UTR(UTRdbraw, flsall, Strandtype="invert")
    } else if (input$seq_strand1 == "Non-stranded") {
      DFUTRraw=PASEXP_3UTR(UTRdbraw, flsall, Strandtype="NONE")
    }
    
    x=nrow(df_pacientes3())/2
    sampleTable1 = data.frame(samplename = gsub(".bam", "", df_pacientes3()$File.Name),
                              condition = c(rep(unique(df_pacientes3()$Sample.Type)[1],x),rep(unique(df_pacientes3()$Sample.Type)[2],x)))
    sampleTable1$samplename <- paste0(sampleTable1$samplename, ".trim")
    
    #x=nrow(df_pacientes3())/2
    #sampleTable1 = data.frame(samplename = c(names(flsall)),
    #                          condition = c(rep("KD",x),rep("NT",x))) #OTIMIZAR, se são 2 amostras, 1 pro KD e 1 pro NT
    
    if (input$multi_test == "Unpaired t-test") {
      NvsT_APA=APAdiff(sampleTable1,DFUTRraw, 
                       conKET=unique(df_pacientes3()$Sample.Type)[2],
                       trtKEY=unique(df_pacientes3()$Sample.Type)[1],
                       PAS='3UTR', CUTreads=5, p_adjust_methods="fdr", MultiTest = "unpaired t-test")
    } else if (input$multi_test == "Paired t-test") {
      NvsT_APA=APAdiff(sampleTable1,DFUTRraw, 
                       conKET=unique(df_pacientes3()$Sample.Type)[2],
                       trtKEY=unique(df_pacientes3()$Sample.Type)[1],
                       PAS='3UTR', CUTreads=5, p_adjust_methods="fdr", MultiTest = "paired t-test")
    } else if (input$multi_test == "ANOVA") {
      NvsT_APA=APAdiff(sampleTable1,DFUTRraw, 
                       conKET=unique(df_pacientes3()$Sample.Type)[2],
                       trtKEY=unique(df_pacientes3()$Sample.Type)[1],
                       PAS='3UTR', CUTreads=5, p_adjust_methods="fdr", MultiTest = "ANOVA")
    }
    
    NvsT_APA <- as.data.frame(NvsT_APA)
    
    shinyalert("Success", "APA Analysis has been completed successfully.", type = "success")
    
    return(NvsT_APA)
  })
  
  Nr_APA_events <- eventReactive(input$run3,{
    #NvsT_APA <- NvsT_APA()
    Nr_APA_events<- table(NvsT_APA()$APAreg)  
    return(Nr_APA_events)
  })
  # NvsT_APA_UP
  NvsT_APA_UP <- eventReactive(input$run3,{
    NvsT_APA <- NvsT_APA()
    NvsT_APA_UP <- NvsT_APA[ which(NvsT_APA$APAreg=='UP'),]
    
    return(NvsT_APA_UP)
  })
  
  # NvsT_APA_DN
  NvsT_APA_DN <- eventReactive(input$run3,{
    NvsT_APA <- NvsT_APA()
    NvsT_APA_DN <- NvsT_APA[ which(NvsT_APA$APAreg=='DN'),]
    
    return(NvsT_APA_DN)
  })
  
  # NvsT_APA_NC
  NvsT_APA_NC <- eventReactive(input$run3,{
    NvsT_APA <- NvsT_APA()
    NvsT_APA_NC <- NvsT_APA[ which(NvsT_APA$APAreg=='NC'),]
    
    return(NvsT_APA_NC)
  })
  
  a <- eventReactive(input$run3,{
    NvsT_APA <- NvsT_APA()
    a <- APAVolcano(NvsT_APA, PAS='3UTR', Pcol = "pvalue", top=40)
    
    return(a)
  })
  b <- eventReactive(input$run3,{
    NvsT_APA <- NvsT_APA()
    b <- APAVolcano(NvsT_APA, PAS='3UTR', Pcol = "pvalue")
    
    return(b)
  })
  d <- eventReactive(input$run3,{
    NvsT_APA <- NvsT_APA()
    d <- APABox(NvsT_APA, xlab = "APAreg", ylab = "RED", plot_title = NULL)
    
    return(d)
  })
  
  
  
  
  output$data_table_2 <- renderTable({
    Nr_APA_events()
  })
  
  output$data_table_3 <- renderTable({
    if (input$search_term3 != "") {
      NvsT_APA_UP() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term3, ignore_case = TRUE))))
    } else {
      NvsT_APA_UP()
    }
  })
  
  output$data_table_4 <- renderTable({
    if (input$search_term4 != "") {
      NvsT_APA_DN() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term4, ignore_case = TRUE))))
    } else {
      NvsT_APA_DN()
    }
  })
  
  output$data_table_5 <- renderTable({
    if (input$search_term5 != "") {
      NvsT_APA_NC() %>%
        filter_all(any_vars(str_detect(., regex(input$search_term5, ignore_case = TRUE))))
    } else {
      NvsT_APA_NC()
    }
  })
  
  
  
  output$plot1 <- renderPlot({
    a()
  }, width = 1200, height = 750)
  output$plot2 <- renderPlot({
    b()
  }, width = 1200, height = 750)
  output$plot3 <- renderPlot({
    d()
  }, width = 1200, height = 750)
  
  
  
  
  output$download_data2 <- downloadHandler(
    filename = function() {
      paste("NvsT_APA_UP", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_APA_UP(), file, row.names = FALSE)
    }
  )
  output$download_data3 <- downloadHandler(
    filename = function() {
      paste("NvsT_APA_DN", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_APA_DN(), file, row.names = FALSE)
    }
  )
  output$download_data4 <- downloadHandler(
    filename = function() {
      paste("NvsT_APA_NC", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NvsT_APA_NC(), file, row.names = FALSE)
    }
  )
  
  output$download_plot1 <- downloadHandler(
    filename = function() {
      paste("plot1", switch(input$file_format, "PNG" = ".png", "PDF" = ".pdf"), sep = "")
    },
    content = function(file) {
      if (input$file_format == "PNG") {
        ggsave(file, plot = a(), width = 5000, height = 5000, units = "px", bg = "white", dpi = 600)
      } else if (input$file_format == "PDF") {
        ggsave(file, plot = a(), width = 7, height = 7, device = "pdf", bg = "white")
      }
    }
  )
  
  output$download_plot2 <- downloadHandler(
    filename = function() {
      paste("plot2", switch(input$file_format, "PNG" = ".png", "PDF" = ".pdf"), sep = "")
    },
    content = function(file) {
      if (input$file_format == "PNG") {
        ggsave(file, plot = b(), width = 5000, height = 5000, units = "px", bg = "white", dpi = 600)
      } else if (input$file_format == "PDF") {
        ggsave(file, plot = b(), width = 7, height = 7, device = "pdf", bg = "white")
      }
    }
  )
  
  output$download_plot3 <- downloadHandler(
    filename = function() {
      paste("plot3", switch(input$file_format, "PNG" = ".png", "PDF" = ".pdf"), sep = "")
    },
    content = function(file) {
      if (input$file_format == "PNG") {
        ggsave(file, plot = d(), width = 5000, height = 5000, units = "px", bg = "white", dpi = 600)
      } else if (input$file_format == "PDF") {
        ggsave(file, plot = d(), width = 7, height = 7, device = "pdf", bg = "white")
      }
    }
  )
  
  
  ##### DGE #####
  
  df_pacientes_dge <- eventReactive(input$run_dge,{
    
    rm(list = ls(all.names = TRUE), envir = .GlobalEnv)
    
    datapath <- input$path_dge
    
    #file <- input$txt_file_dge
    if (datapath == "") {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    
    if (nrow(rv$data) == 0) {
      shinyalert("Error", "Please make sure that the sample sheet is correct before running the analysis.", type = "error")
      return(NULL)
    }
    
    df_pacientes_dge <- rv$data
    
    #df_pacientes_dge <- read.table(file$datapath, sep = "\t", header = TRUE)
    df_pacientes_dge <- dplyr::select(df_pacientes_dge, File.Name,Sample.Type)
    
    #df_pacientes_dge <- read.table(file$datapath, sep = "\t", header = TRUE)
    #df_pacientes_dge <- dplyr::select(df_pacientes_dge, File.Name, Case.ID, Sample.Type)
    df_pacientes_dge$File.Name<- str_replace(df_pacientes_dge$File.Name, ".bam", ".trim.htseq.txt")
    
    #df_pacientes_dge$Sample.Type <- str_replace(df_pacientes_dge$Sample.Type, "Primary Tumor", "PrimaryTumor")
    #df_pacientes_dge$Sample.Type <- str_replace(df_pacientes_dge$Sample.Type, "Solid Tissue Normal", "NormalTissue")
    #normal <- c("NormalTissue")
    #df_pacientes_dge$category <- ifelse(df_pacientes_dge$Sample.Type %in% c("NormalTissue"), "NormalTissue", "PrimaryTumor")
    #df_pacientes_dge$category2 = paste(df_pacientes_dge$Case.ID, df_pacientes_dge$category, sep="_")
    #df_pacientes_dge$category3 = paste(df_pacientes_dge$File.Name, df_pacientes_dge$Sample.Type, sep="_")
    df_pacientes_dge$Sample.Type <- gsub(" ", "_", df_pacientes_dge$Sample.Type)
      
    return(df_pacientes_dge)
  })
  
  dds <- eventReactive(input$run_dge,{
    dir <- input$path_dge
    if (is.null(dir)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    
    shinyalert("Processing", "DGE Analysis has started. Please wait...", type = "info",
               closeOnEsc = FALSE, closeOnClickOutside = FALSE, showConfirmButton = TRUE)
    
    dir <- str_replace_all(dir, "\\\\", "/")
    setwd(dir)
    getwd()
    sampleFiles=grep('.htseq.txt', list.files(dir), value=TRUE)
    sampleNames=names(sampleFiles) <- df_pacientes_dge()$Sample.Type
    sampleCondition=gsub("[a-zA-Z0-9]*.trim.htseq.txt_",'',sampleNames)
    #sampleCondition=gsub("\\d[A-Za-z].sorted.htseq.txt_",'',sampleNames)
    
    sampleTable=data.frame(sampleName = df_pacientes_dge()$File.Name, fileName = sampleFiles, condition = sampleCondition)
    
    reorder_idx <- match(sampleTable$sampleName,sampleTable$fileName) 
    sampleTable$fileName <- sampleTable$fileName[reorder_idx]
    
    sampleTable <- sampleTable[order(sampleTable$condition),]
    #reorder_idx <- match(sampleTable$sampleName, sampleTable$fileName) 
    #sampleTable$fileName <- sampleTable$fileName[reorder_idx]
    #sampleTable <- sampleTable[order(sampleTable$condition), ]
    
    dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                      directory = dir,
                                      design = ~ condition)
    
    
    keep <- rowSums(counts(dds)) >= 5
    dds2 <- dds[keep, ]
    dds2$condition <- factor(dds2$condition, levels = c(unique(df_pacientes_dge()$Sample.Type)[2],unique(df_pacientes_dge()$Sample.Type)[1]))
    
    return(dds2)
  })
  
  vst_dds <- eventReactive(input$run_dge,{
    dds2 <- dds()
    vst_dds <- vst(dds2)
    
    return(vst_dds)
  })
  
  pca_plot <- eventReactive(input$run_dge,{
    vst_dds <- vst_dds()
    pca_plot <- plotPCA(vst_dds, intgroup = "condition",)
    
    return(pca_plot)
  })
  
  ddx <- eventReactive(input$run_dge,{
    dds2 <- dds()
    
    ddx <- DESeq(dds2)
    ddx <- estimateSizeFactors(ddx)
    
    return(ddx)
  })
  
  res <- reactive({
    ddx <- ddx()
    
    res <- results(ddx, contrast=c("condition",unique(df_pacientes_dge()$Sample.Type)[2],unique(df_pacientes_dge()$Sample.Type)[1]))
    res <- as.data.frame(res)
    
    return(res)
  })
  
  normalized_counts <- reactive({
    ddx <- ddx()
    
    normalized_counts <- counts(ddx, normalized=TRUE)
    
    return(normalized_counts)
  })
  
  resLFC <- eventReactive(input$run_dge,{
    ddx <- ddx()
    
    cond <- paste0("condition_", unique(df_pacientes_dge()$Sample.Type)[1], "_vs_", unique(df_pacientes_dge()$Sample.Type)[2])  
    
    resLFC <- lfcShrink(ddx, coef=cond, type="apeglm")
    
    return(resLFC)
  })
  
  
  vsd <- eventReactive(input$run_dge,{
    ddx <- ddx()
    
    vst(ddx, blind=FALSE)
    
    return(vst)
  })
  
  debounced_cellwidth <- debounce(reactive(input$cellwidth), 1000)  # 500 ms delay
  debounced_cellheight <- debounce(reactive(input$cellheight), 1000)  # 500 ms delay
  
  htmap <- reactive({
    req(input$run_dge)
    
    res <- res()
    normalized_counts <- normalized_counts()
    
    signi <- subset(res, (padj <= 0.05))
    allSig <- merge(normalized_counts, signi, by = 0)
    sigCounts <- allSig[, 2:(ncol(allSig) - 6)]
    row.names(sigCounts) <- allSig$Row.names
    
    htmap <- pheatmap(log2(sigCounts+1), scale = "row", cluster_rows=TRUE, show_rownames=FALSE, show_colnames=TRUE,
             cluster_cols=TRUE, treeheight_row = 0, treeheight_col = 50,  display_numbers=FALSE,
             color = colorRampPalette(c("blue", "white", "red"))(100), cellwidth = debounced_cellwidth(), cellheight = debounced_cellheight())
    
    return(htmap)
  })
  
  Volcano_dge <- eventReactive(input$run_dge,{
    resLFC <- resLFC()
    
    keyvals <- ifelse(
      (resLFC$log2FoldChange < -2 & resLFC$pvalue < 0.05), 'purple',
      ifelse((resLFC$log2FoldChange > 2 & resLFC$pvalue < 0.05), 'darkgreen',
             'black'))
    keyvals[is.na(keyvals)] <- 'black'
    names(keyvals)[keyvals == 'darkgreen'] <- 'Signf. Upregulated'
    names(keyvals)[keyvals == 'black'] <- 'Non-Significant'
    names(keyvals)[keyvals == 'purple'] <- 'Signif. Downregulated'
    Volcano_dge <- EnhancedVolcano(resLFC,
                    lab = rownames(resLFC),
                    labSize = 0.0,
                    x = 'log2FoldChange',
                    y = 'pvalue', xlim = c(-10,10),
                    pCutoff = 0.05,
                    cutoffLineCol = "red",
                    FCcutoff = 2.0,
                    colCustom = keyvals)
    
    return(Volcano_dge)
  })
  
  res05 <- eventReactive(input$run_dge,{
    ddx <- ddx()    
    
    res05 <- results(ddx, alpha=0.05)
    
    res05$DGEreg <- ifelse(res05$log2FoldChange > 2 & res05$padj < 0.05, "UP",
                             ifelse(res05$log2FoldChange < -2 & res05$padj < 0.05, "DN", "NC"))
    
    res05 <- cbind(gene_symbol = rownames(res05), res05)
    
    res05 <- as.data.frame(res05)
    
    shinyalert("Success", "DGE Analysis has been completed successfully.", type = "success")
    
    return(res05)
  })
  
  res05_DGEreg <- eventReactive(input$run_dge,{
    res05 <- res05()
    
    return(res05$DGEreg)
  })
  
  genes_up_05 <- eventReactive(input$run_dge,{
    res05 <- res05()    
    
    #genes_up_05 <- as.data.frame(res05[which(res05$log2FoldChange > 2 & res05$padj < .05),])
    genes_up_05 <- res05[ which(res05$DGEreg=='UP'),]
    
    
    return(genes_up_05)
  })
  
  genes_down_05 <- eventReactive(input$run_dge,{
    res05 <- res05()
    
    #genes_down_05 <- as.data.frame(res05[which(res05$log2FoldChange < -2 & res05$padj < .05),])
    genes_down_05 <- res05[ which(res05$DGEreg=='DN'),]
    
    return(genes_down_05)
  })
  
  genes_nc_05 <- eventReactive(input$run_dge,{
    res05 <- res05()
    
    #genes_nc_05 <- as.data.frame(res05[which(res05$log2FoldChange > -2 & res05$log2FoldChange < 2),])
    genes_nc_05 <- res05[ which(res05$DGEreg=='NC'),]
    
    return(genes_nc_05)
  })
  
  output$data_table_13 <- renderTable({
    table(res05_DGEreg())
  })
  
  output$data_table_14 <- renderTable({
    genes_up_05 <- genes_up_05()
    
    if (input$search_term_dge != "") {
      genes_up_05 %>%
        filter_all(any_vars(str_detect(., regex(input$search_term_dge, ignore_case = TRUE))))
    } else {
      genes_up_05
    }
  })
  
  output$data_table_15 <- renderTable({
    genes_down_05 <- genes_down_05()
    
    if (input$search_term2_dge != "") {
      genes_down_05 %>%
        filter_all(any_vars(str_detect(., regex(input$search_term2_dge, ignore_case = TRUE))))
    } else {
      genes_down_05
    }
  })
  
  output$data_table_16 <- renderTable({
    genes_nc_05 <- genes_nc_05()
    
    if (input$search_term3_dge != "") {
      genes_nc_05 %>%
        filter_all(any_vars(str_detect(., regex(input$search_term3_dge, ignore_case = TRUE))))
    } else {
      genes_nc_05
    }
  })
  
  output$plot7 <- renderPlot({
    pca_plot()
  }, width = 1200, height = 750)
  
  output$plot8 <- renderPlot({
    Volcano_dge()
  }, width = 1200, height = 750)
  
  output$plot9 <- renderPlot({
    htmap()
  }, width = 1200, height = 750)
  
  output$download_plot7 <- downloadHandler(
    filename = function() {
      paste("plot_pca", switch(input$file_format, "PNG" = ".png", "PDF" = ".pdf"), sep = "")
    },
    content = function(file) {
      if (input$file_format == "PNG") {
        ggsave(file, plot = pca_plot(), width = 5000, height = 5000, units = "px", bg = "white", dpi = 600)
      } else if (input$file_format == "PDF") {
        ggsave(file, plot = pca_plot(), width = 7, height = 7, device = "pdf", bg = "white")
      }
    }
  )
  
  output$download_plot8 <- downloadHandler(
    filename = function() {
      paste("plot_volcano", switch(input$file_format, "PNG" = ".png", "PDF" = ".pdf"), sep = "")
    },
    content = function(file) {
      if (input$file_format == "PNG") {
        ggsave(file, plot = Volcano_dge(), width = 5000, height = 5000, units = "px", bg = "white", dpi = 600)
      } else if (input$file_format == "PDF") {
        ggsave(file, plot = Volcano_dge(), width = 7, height = 7, device = "pdf", bg = "white")
      }
    }
  )
  
  output$download_plot9 <- downloadHandler(
    filename = function() {
      paste("plot_heatmap", switch(input$file_format, "PNG" = ".png", "PDF" = ".pdf"), sep = "")
    },
    content = function(file) {
      if (input$file_format == "PNG") {
        ggsave(file, plot = htmap(), width = 5000, height = 5000, units = "px", bg = "white", dpi = 600)
      } else if (input$file_format == "PDF") {
        ggsave(file, plot = htmap(), width = 7, height = 7, device = "pdf", bg = "white")
      }
    }
  )
  
  output$download_data_dge <- downloadHandler(
    filename = function() {
      paste("DGE_Genes_UP", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(genes_up_05(), file, row.names = FALSE)
    }
  )
  
  output$download_data2_dge <- downloadHandler(
    filename = function() {
      paste("DGE_Genes_DN", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(genes_down_05(), file, row.names = FALSE)
    }
  )
  
  output$download_data3_dge <- downloadHandler(
    filename = function() {
      paste("DGE_Genes_NC", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(genes_nc_05(), file, row.names = FALSE)
    }
  )
  
  
  ##### GO TERMS #####
  
  de <- eventReactive(input$run_go,{
    
    rm(list = ls(all.names = TRUE), envir = .GlobalEnv)
    
    file <- input$txt_file_go
    if (is.null(file)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    
    gene_list <- read.csv(file$datapath, header = TRUE)
    #gene_list <- gene_list[, 2]
    gene_list <- gene_list$gene_symbol
    
    return(gene_list)
  })
  
  
  ego_BP <- eventReactive(input$run_go,{
    de <- de()
    #pvalue_cutoff <- as.numeric(input$output_type2_go)
    shinyalert("Processing", "GO Analysis has started. Please wait...", type = "info",
               closeOnEsc = FALSE, closeOnClickOutside = FALSE, showConfirmButton = TRUE)
    
    if (input$database_go == "Human") {
      ego <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", pvalueCutoff = Inf)
    } else if (input$database_go == "Mouse") {
      ego <- enrichGO(gene = de, OrgDb = "org.Mm.eg.db", ont = "BP", keyType = "SYMBOL", pvalueCutoff = Inf)
    }
    
    ego@result$GeneRatio_num <- sapply(strsplit(as.character(ego@result$GeneRatio), "/"), function(x) as.numeric(x[1]))
    ego@result$GeneRatio_den <- sapply(strsplit(as.character(ego@result$GeneRatio), "/"), function(x) as.numeric(x[2]))
    ego@result$BgRatio_num <- sapply(strsplit(as.character(ego@result$BgRatio), "/"), function(x) as.numeric(x[1]))
    ego@result$BgRatio_den <- sapply(strsplit(as.character(ego@result$BgRatio), "/"), function(x) as.numeric(x[2]))
    ego@result$fe <- ego@result$GeneRatio_num/ego@result$GeneRatio_den / (ego@result$BgRatio_num/ego@result$BgRatio_den)
    
    shinyalert("Success", "GO Analysis has been completed successfully.", type = "success")
    
    return(ego)
  })
  
  ego_MF <- eventReactive(input$run_go,{
    de <- de()
    #pvalue_cutoff <- as.numeric(input$output_type2_go)
    
    if (input$database_go == "Human") {
      ego <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db", ont = "MF", keyType = "SYMBOL", pvalueCutoff = Inf)
    } else if (input$database_go == "Mouse") {
      ego <- enrichGO(gene = de, OrgDb = "org.Mm.eg.db", ont = "MF", keyType = "SYMBOL", pvalueCutoff = Inf)
    }
    
    ego@result$GeneRatio_num <- sapply(strsplit(as.character(ego@result$GeneRatio), "/"), function(x) as.numeric(x[1]))
    ego@result$GeneRatio_den <- sapply(strsplit(as.character(ego@result$GeneRatio), "/"), function(x) as.numeric(x[2]))
    ego@result$BgRatio_num <- sapply(strsplit(as.character(ego@result$BgRatio), "/"), function(x) as.numeric(x[1]))
    ego@result$BgRatio_den <- sapply(strsplit(as.character(ego@result$BgRatio), "/"), function(x) as.numeric(x[2]))
    ego@result$fe <- ego@result$GeneRatio_num/ego@result$GeneRatio_den / (ego@result$BgRatio_num/ego@result$BgRatio_den)
    
    return(ego)
  })
  
  GO_plot_BP <- eventReactive(input$run_go,{
    ego <- ego_BP()
    
    p <- dotplot(ego, x = "fe", color = "p.adjust", showCategory=20) + scale_color_gradient(low = "red", high = "tan") + labs(size="Count", colour="P.adjust") + xlab("Fold Enrichment")
    p <- p + ggtitle("GO Terms BP")
    p <- p + theme(plot.title = element_text(size = 20, hjust = 0.5),
                   axis.text.x = element_text(size = 15),
                   axis.title.x = element_text(size = 17),
                   axis.text.y = element_text(size = 12),
                   legend.text = element_text(size = 12),
                   legend.title = element_text(size = 17))
    
    return(p)
  })
  
  GO_plot_MF <- eventReactive(input$run_go,{
    ego <- ego_MF()
    
    p <- dotplot(ego, x = "fe", color = "p.adjust", showCategory=20) + scale_color_gradient(low = "#DE2142", high = "tan") + labs(size="Count", colour="P.adjust") + xlab("Fold Enrichment")
    p <- p + ggtitle("GO Terms MF")
    p <- p + theme(plot.title = element_text(size = 20, hjust = 0.5),
                   axis.text.x = element_text(size = 15),
                   axis.title.x = element_text(size = 17),
                   axis.text.y = element_text(size = 12),
                   legend.text = element_text(size = 12),
                   legend.title = element_text(size = 17))
    
    return(p)
  })
  
  output$plot_go <- renderPlot({
    GO_plot_BP()
  }, width = 1200, height = 750)
  
  output$plot2_go <- renderPlot({
    GO_plot_MF()
  }, width = 1200, height = 750)
  
  output$download_plot_go <- downloadHandler(
    filename = function() {
      paste("plot_go_BP", switch(input$file_format, "PNG" = ".png", "PDF" = ".pdf"), sep = "")
    },
    content = function(file) {
      if (input$file_format == "PNG") {
        ggsave(file, plot = GO_plot_BP(), width = 5000, height = 5000, units = "px", bg = "white", dpi = 600)
      } else if (input$file_format == "PDF") {
        ggsave(file, plot = GO_plot_BP(), width = 7, height = 7, device = "pdf", bg = "white")
      }
    }
  )
  
  output$download_plot2_go <- downloadHandler(
    filename = function() {
      paste("plot_go_MF", switch(input$file_format, "PNG" = ".png", "PDF" = ".pdf"), sep = "")
    },
    content = function(file) {
      if (input$file_format == "PNG") {
        ggsave(file, plot = GO_plot_MF(), width = 5000, height = 5000, units = "px", bg = "white", dpi = 600)
      } else if (input$file_format == "PDF") {
        ggsave(file, plot = GO_plot_MF(), width = 7, height = 7, device = "pdf", bg = "white")
      }
    }
  )
  
  
  ##### VENN #####
  
  input_files <- eventReactive(input$run_venn,{
    
    rm(list = ls(all.names = TRUE), envir = .GlobalEnv)
    
    files <- input$txt_file_venn
    
    if (is.null(files)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    } else if (nrow(files) == 1) {
      shinyalert("Error", "Please input more than one file for intersection", type = "error")
      return(NULL)
    } else if (nrow(files) == 2) {
      #df1 <- read.csv(files[4][1, ], header = TRUE)
      #df2 <- read.csv(files[4][2, ], header = TRUE)
      
      df1 <- read.csv(files$datapath[1], header = TRUE)
      df2 <- read.csv(files$datapath[2], header = TRUE)
      
      x <- list(
        dataset1 = df1$gene_symbol,
        dataset2 = df2$gene_symbol
      )
      
      selected_df1 <- df1[, 1, drop = FALSE]
      selected_df2 <- df2[, 1, drop = FALSE]
      
      common_genes <- merge(selected_df1, selected_df2, by.x = 1, by.y = 1)
      common_genes <- unique(common_genes)
      
      plot <- ggvenn(x, fill_color = c("#999999", "#999555"), show_percentage = FALSE, text_size = 7)
      
      return(list(plot = plot, common_genes = common_genes))
    } else if (nrow(files) == 3) {
      #df1 <- read.csv(files[4][1, ], header = TRUE)
      #df2 <- read.csv(files[4][2, ], header = TRUE)
      #df3 <- read.csv(files[4][3, ], header = TRUE)
      
      df1 <- read.csv(files$datapath[1], header = TRUE)
      df2 <- read.csv(files$datapath[2], header = TRUE)
      df3 <- read.csv(files$datapath[3], header = TRUE)
      
      x <- list(
        dataset1 = df1$gene_symbol,
        dataset2 = df2$gene_symbol,
        dataset3 = df3$gene_symbol
      )
      
      selected_df1 <- df1[, 1, drop = FALSE]
      selected_df2 <- df2[, 1, drop = FALSE]
      selected_df3 <- df3[, 1, drop = FALSE]
      
      common_genes <- merge(selected_df1, selected_df2, by.x = 1, by.y = 1)
      common_genes <- merge(common_genes, selected_df3, by.x = 1, by.y = 1)
      common_genes <- unique(common_genes)
      
      plot <- ggvenn(x, fill_color = c("#999999", "#E69F00", "#56B4E9"), show_percentage = FALSE, text_size = 7)
      
      return(list(plot = plot, common_genes = common_genes))
    } else if (nrow(files) == 4) {
      #df1 <- read.csv(files[4][1, ], header = TRUE)
      #df2 <- read.csv(files[4][2, ], header = TRUE)
      #df3 <- read.csv(files[4][3, ], header = TRUE)
      #df4 <- read.csv(files[4][4, ], header = TRUE)
      
      df1 <- read.csv(files$datapath[1], header = TRUE)
      df2 <- read.csv(files$datapath[2], header = TRUE)
      df3 <- read.csv(files$datapath[3], header = TRUE)
      df4 <- read.csv(files$datapath[4], header = TRUE) 
      
      x <- list(
        dataset1 = df1$gene_symbol,
        dataset2 = df2$gene_symbol,
        dataset3 = df3$gene_symbol,
        dataset4 = df4$gene_symbol
      )
      
      selected_df1 <- df1[, 1, drop = FALSE]
      selected_df2 <- df2[, 1, drop = FALSE]
      selected_df3 <- df3[, 1, drop = FALSE]
      selected_df4 <- df4[, 1, drop = FALSE]
      
      common_genes <- merge(selected_df1, selected_df2, by.x = 1, by.y = 1)
      common_genes <- merge(common_genes, selected_df3, by.x = 1, by.y = 1)
      common_genes <- merge(common_genes, selected_df4, by.x = 1, by.y = 1)
      common_genes <- unique(common_genes)
      
      plot <- ggvenn(x, fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"), show_percentage = FALSE, text_size = 7)
      
      return(list(plot = plot, common_genes = common_genes))
    } else if (nrow(files) == 5) {
      # Helper function to display Venn diagram
      display_venn <- function(x, ...){
        grid.newpage()
        venn_object <- venn.diagram(x, filename = NULL, ...)
        grid.draw(venn_object)
      }
      
      #df1 <- read.csv(files[4][1, ], header = TRUE)
      #df2 <- read.csv(files[4][2, ], header = TRUE)
      #df3 <- read.csv(files[4][3, ], header = TRUE)
      #df4 <- read.csv(files[4][4, ], header = TRUE)
      #df5 <- read.csv(files[4][5, ], header = TRUE)
      
      df1 <- read.csv(files$datapath[1], header = TRUE)
      df2 <- read.csv(files$datapath[2], header = TRUE)
      df3 <- read.csv(files$datapath[3], header = TRUE)
      df4 <- read.csv(files$datapath[4], header = TRUE)
      df5 <- read.csv(files$datapath[5], header = TRUE)
      
      x <- list(
        dataset1 = df1$gene_symbol,
        dataset2 = df2$gene_symbol,
        dataset3 = df3$gene_symbol,
        dataset4 = df4$gene_symbol,
        dataset5 = df5$gene_symbol
      )
      
      selected_df1 <- df1[, 1, drop = FALSE]
      selected_df2 <- df2[, 1, drop = FALSE]
      selected_df3 <- df3[, 1, drop = FALSE]
      selected_df4 <- df4[, 1, drop = FALSE]
      selected_df5 <- df5[, 1, drop = FALSE]
      
      common_genes <- merge(selected_df1, selected_df2, by.x = 1, by.y = 1)
      common_genes <- merge(common_genes, selected_df3, by.x = 1, by.y = 1)
      common_genes <- merge(common_genes, selected_df4, by.x = 1, by.y = 1)
      common_genes <- merge(common_genes, selected_df5, by.x = 1, by.y = 1)
      common_genes <- unique(common_genes)
      
      display_venn(x, fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#999555"), alpha = 0.7, margin = 0.02, cex = 1.3)
      
      return(list(plot = plot, common_genes = common_genes))
    }
  })
  
  output$plot_venn <- renderPlot({
    req(input_files()$plot)
  }, width = 1200, height = 750)
  
  output$common_genes <- renderTable({
    common_genes <- input_files()$common_genes
    
    if (input$search_term_common_genes != "") {
      common_genes %>%
        filter_all(any_vars(str_detect(., regex(input$search_term_common_genes, ignore_case = TRUE))))
    } else {
      common_genes
    }
  })
  
  output$download_plot_venn <- downloadHandler(
    filename = function() {
      paste("plot_venn", switch(input$file_format, "PNG" = ".png", "PDF" = ".pdf"), sep = "")
    },
    content = function(file) {
      if (input$file_format == "PNG") {
        ggsave(file, plot = input_files()$plot, width = 5000, height = 5000, units = "px", bg = "white", dpi = 600)
      } else if (input$file_format == "PDF") {
        ggsave(file, plot = input_files()$plot, width = 7, height = 7, device = "pdf", bg = "white")
      }
    }
  )
  
  output$download_common_genes <- downloadHandler(
    filename = function() {
      paste("common_genes", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(input_files()$common_genes, file, row.names = FALSE)
    }
  )
  
  
  ##### APA CORR #####
  plot_corr_apa <- eventReactive(input$run_corr,{
    files1 <- input$corr_apa
    files2 <- input$corr_dge
    if (is.null(files1) || is.null(files2)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    
    shinyalert("Processing", "APA & DGE correlation analysis has started. Please wait...", type = "info",
               closeOnEsc = FALSE, closeOnClickOutside = FALSE, showConfirmButton = TRUE)
    
    df_APA_len <- read.csv(files1$datapath[1])
    df_APA_short <- read.csv(files1$datapath[2])
    
    df_DGE_up <- read.csv(files2$datapath[1])
    df_DGE_nc <- read.csv(files2$datapath[2])
    df_DGE_dn <- read.csv(files2$datapath[3])
    
    merged_APA <- rbind(df_APA_len, df_APA_short)
    merged_DGE <- rbind(df_DGE_up, df_DGE_nc, df_DGE_dn)
    
    merged_all <- merge(merged_APA, merged_DGE, by.x = "gene_symbol", by.y = "gene_symbol")
    merged_all <- merged_all[, c("gene_symbol", "RED", "log2FoldChange", "APAreg")]
    
    df <- data.frame(
      RED = merged_all$RED,
      log2foldchange = merged_all$log2FoldChange,
      APAreg = merged_all$APAreg
    )
    
    correlation_coefficient <- cor.test(df$log2foldchange, df$RED, method = "pearson")
    original_number <- correlation_coefficient$p.value
    
    plot <- ggplot(df, aes(x = RED, y = log2foldchange, color = APAreg)) +
      geom_point() +
      geom_smooth(method = "lm", se = TRUE, color = "black") +
      labs(
        title = paste("3'UTR-APA vs DGE\n\n", "r:", round(correlation_coefficient$estimate, 4), "         ", "p:", format(original_number, scientific = TRUE, digits = 4)),
        x = "RED",
        y = "log2FC DGE"
      ) +
      scale_color_manual(values = c("UP" = "red", "DN" = "blue"))
    
    shinyalert("Success", "APA & DGE correlation analysis has been completed successfully.", type = "success")
    
    return(plot)
  })
  
  output$plot_corr <- renderPlot({
    req(plot_corr_apa())
  }, width = 1200, height = 750)
  
  output$download_plot_corr <- downloadHandler(
    filename = function() {
      paste("plot_corr_3'UTR-APA_DGE", switch(input$file_format, "PNG" = ".png", "PDF" = ".pdf"), sep = "")
    },
    content = function(file) {
      if (input$file_format == "PNG") {
        ggsave(file, plot = plot_corr_apa(), width = 5000, height = 5000, units = "px", bg = "white", dpi = 600)
      } else if (input$file_format == "PDF") {
        ggsave(file, plot = plot_corr_apa(), width = 7, height = 7, device = "pdf", bg = "white")
      }
    }
  )
  
  #output$download_plot_corr <- downloadHandler(
  #  filename = function() {
  #    paste("plot_corr_3'UTR-APA_DGE", ".png", sep = "")
  #  },
  #  content = function(file) {
  #    ggsave(file, plot = plot_corr_apa(), width = 5000, height = 5000, units = c("px"), bg = "white", dpi = 600)
  #  }
  #
  
  
  ##### IPA CORR #####
  plot_corr_ipa <- eventReactive(input$run_corr2,{
    files1 <- input$corr_ipa
    files2 <- input$corr_dge2
    if (is.null(files1) || is.null(files2)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    
    shinyalert("Processing", "IPA & DGE correlation analysis has started. Please wait...", type = "info",
               closeOnEsc = FALSE, closeOnClickOutside = FALSE, showConfirmButton = TRUE)
    
    df_IPA_up <- read.csv(files1$datapath[1])
    df_IPA_up <- df_IPA_up %>% distinct(gene_symbol, .keep_all = TRUE)
    df_IPA_dn <- read.csv(files1$datapath[2])
    df_IPA_dn <- df_IPA_dn %>% distinct(gene_symbol, .keep_all = TRUE)
    
    df_DGE_up <- read.csv(files2$datapath[1])
    df_DGE_nc <- read.csv(files2$datapath[2])
    df_DGE_dn <- read.csv(files2$datapath[3])
    
    merged_IPA <- rbind(df_IPA_up, df_IPA_dn)
    merged_DGE <- rbind(df_DGE_up, df_DGE_nc, df_DGE_dn)
    
    merged_all <- merge(merged_IPA, merged_DGE, by.x = "gene_symbol", by.y = "gene_symbol")
    merged_all <- merged_all[, c("gene_symbol", "RED", "log2FoldChange", "APAreg")]
    
    df <- data.frame(
      RED = merged_all$RED,
      log2foldchange = merged_all$log2FoldChange,
      APAreg = merged_all$APAreg
    )
    
    correlation_coefficient <- cor.test(df$log2foldchange, df$RED, method = "pearson")
    original_number <- correlation_coefficient$p.value
    
    plot <- ggplot(df, aes(x = RED, y = log2foldchange, color = APAreg)) +
      geom_point() +
      geom_smooth(method = "lm", se = TRUE, color = "black") +
      labs(
        title = paste("IPA vs DGE\n\n", "r:", round(correlation_coefficient$estimate, 4), "         ", "p:", format(original_number, scientific = TRUE, digits = 4)),
        x = "RED",
        y = "log2FC DGE"
      ) +
      scale_color_manual(values = c("UP" = "red", "DN" = "blue"))
    
    shinyalert("Success", "IPA & DGE correlation analysis has been completed successfully.", type = "success")
    
    return(plot)
  })
  
  output$plot_corr2 <- renderPlot({
    req(plot_corr_ipa())
  }, width = 1200, height = 750)
  
  output$download_plot_corr2 <- downloadHandler(
    filename = function() {
      paste("plot_corr_IPA_DGE", switch(input$file_format, "PNG" = ".png", "PDF" = ".pdf"), sep = "")
    },
    content = function(file) {
      if (input$file_format == "PNG") {
        ggsave(file, plot = plot_corr_ipa(), width = 5000, height = 5000, units = "px", bg = "white", dpi = 600)
      } else if (input$file_format == "PDF") {
        ggsave(file, plot = plot_corr_ipa(), width = 7, height = 7, device = "pdf", bg = "white")
      }
    }
  )
}

# Run the app
shinyApp(ui, server)
