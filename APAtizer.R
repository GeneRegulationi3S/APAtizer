#https://bioconductor.org/packages/devel/bioc/vignettes/APAlyzer/inst/doc/APAlyzer.html#analysis-of-apa-in-3utrs
#Apalyzer 
#Installing Apalyzer

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("APAlyzer")
BiocManager::install("Rsamtools")
library("TCGAbiolinks")
library("survminer")
library("survival")
library("SummarizedExperiment")
library("tidyverse")
library("DESeq2")
library("APAlyzer")
library("base")
library("data.table")
library("Rsamtools")
library("tidyverse")
library("stats")
library("shinyalert")
library("shiny")
library("purrr")
library("data.table")
library("base")
library("tidyverse")
library("dplyr")
library("splitstackshape")
library("shinythemes")
library("readr")
library("ggplot2")
library("repmis")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("EnhancedVolcano")
library("clusterProfiler")
library("enrichplot")
library("readr")
library("org.Hs.eg.db")
library("VennDiagram")
library("ggvenn")

options(shiny.maxRequestSize=10000000*1024^2)

ui <- fluidPage(theme = shinytheme("darkly"),
                navbarPage(
                  "APAtizer",
                  tabPanel("DAPARS",
                           titlePanel("Use DaPars2 to analyse 3'UTR-APA from RNA-Seq data"),
                           sidebarLayout(
                             sidebarPanel(
                               fileInput("txt_files", "Select Multiple .txt Files", multiple = TRUE),
                               fileInput("txt_file2", "Select Single .txt File"),
                               actionButton("run2", "DaPars Analysis")
                             ),
                             mainPanel(
                               tabsetPanel(
                                 tabPanel("Len genes",
                                          br(),
                                          textInput(inputId = "search_term", label = "Search for a gene"),
                                          downloadButton("download_datax", "Download Len Gene Data"),
                                          tableOutput("data_table_x")
                                 ),
                                 tabPanel("Short Genes",
                                          br(),
                                          textInput(inputId = "search_term2", label = "Search for a gene"),
                                          downloadButton("download_datay", "Download Short Gene Data"),
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
                               textInput("path2", "Directory path:"),
                               fileInput("txt_file3", "Select Single .txt File"),
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
                                            downloadButton(outputId = "download_plot1", label = "Download APA Volcano plot (top 40)"),
                                            plotOutput(outputId = "plot1")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type5 == 'APA Volcano plot'",
                                            br(),
                                            downloadButton(outputId = "download_plot2", label = "Download APA Volcano plot"),
                                            plotOutput(outputId = "plot2")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type5 == 'APA Box'",
                                            br(),
                                            downloadButton(outputId = "download_plot3", label = "Download APA Box"),
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
                               textInput("path", "Directory path:"),
                               fileInput("txt_file", "Select Single .txt File"),
                               actionButton("run", "IPA Analysis"),
                               br(),
                               br(),
                               br(),
                               selectInput(inputId = "output_type", label = "Select Output Type",
                                           choices = c("NvsT_IPA_events_UP", "NvsT_IPA_events_DN","NvsT_IPA_events_NC" )),
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
                                            downloadButton(outputId = "download_plot4", label = "Download IPA Volcano plot (top 40)"),
                                            plotOutput(outputId = "plot4")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type3 == 'IPA Volcano plot'",
                                            br(),
                                            downloadButton(outputId = "download_plot5", label = "Download IPA Volcano plot"),
                                            plotOutput(outputId = "plot5")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type3 == 'IPA Box'",
                                            br(),
                                            downloadButton(outputId = "download_plot6", label = "Download IPA Box"),
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
                               textInput("path_dge", "Directory path:"),
                               fileInput("txt_file_dge", "Select Single .txt File"),
                               actionButton("run_dge", "DGE Analysis"),
                               br(),
                               br(),
                               br(),
                               selectInput(inputId = "output_type_dge", label = "Select Output Type",
                                           choices = c("DGE_Genes_UP", "DGE_Genes_DN", "DGE_Genes_NC")),
                               br(),
                               selectInput(inputId = "output_type2_dge", label = "Select Output Type",
                                           choices = c("PCA Plot", "DGE Volcano Plot", "DGE Heatmap"))
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
                                            downloadButton(outputId = "download_plot7", label = "Download PCA Plot"),
                                            plotOutput(outputId = "plot7")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type2_dge == 'DGE Volcano Plot'",
                                            br(),
                                            downloadButton(outputId = "download_plot8", label = "Download DGE Volcano Plot"),
                                            plotOutput(outputId = "plot8")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type2_dge == 'DGE Heatmap'",
                                            br(),
                                            downloadButton(outputId = "download_plot9", label = "Download DGE Heatmap"),
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
                                            downloadButton(outputId = "download_plot_go", label = "Download GO Plot BP"),
                                            plotOutput(outputId = "plot_go")
                                          ),
                                          conditionalPanel(
                                            condition = "input.output_type_go == 'Molecular Function (MF)'",
                                            br(),
                                            downloadButton(outputId = "download_plot2_go", label = "Download GO Plot MF"),
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
                                            downloadButton(outputId = "download_plot_venn", label = "Download Venn Diagram"),
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
                  tabPanel("SURVIVAL ANALYSIS",
                           titlePanel("Perform Survival analysis on DGE genes"),
                           sidebarLayout(
                             sidebarPanel(
                               fileInput("surv_sample_sheet", "Select sample sheet"),
                               fileInput("surv_clinical_data", "Select clinical data"),
                               textInput(inputId = "path_surv", "HTSeq files directory path:"),
                               textInput(inputId = "surv_gene", label = "Gene of interest"),
                               br(),
                               actionButton("run_surv", "Survival Analysis")
                             ),
                             mainPanel(
                               tabsetPanel(
                                 tabPanel("Plot",
                                          br(),
                                          downloadButton(outputId = "download_plot_surv", label = "Download Survival Analysis Plot"),
                                          plotOutput(outputId = "plot_surv"))
                               )
                             )
                           )
                  )
          )
)
                                            
                                            

server <- function(input, output,session) {
  
  
  ##### IPA #####
  
  df_pacientes <- eventReactive(input$run,{
    file <- input$txt_file
    if (is.null(file)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    df_pacientes <- read.table(file$datapath, sep = "\t", header = TRUE)
    df_pacientes <- dplyr::select(df_pacientes, File.Name,Case.ID,Sample.Type)
    df_pacientes$File.Name<- str_replace(df_pacientes$File.Name, "rna_seq.genomic.gdc_realn.bam", "rna_seq.genomic.gdc_realn.trim.bam")
    
    case_order <- unique(df_pacientes$Case.ID)
    df_pacientes <- df_pacientes %>% arrange(df_pacientes$Sample.Type, match(df_pacientes$Case.ID, case_order))
    
    #normal <- c("Solid Tissue Normal")
    df_pacientes$category <- ifelse(df_pacientes$Sample.Type %in% c("Solid Tissue Normal"), "Normal", "Tumor")
    df_pacientes$category = paste(df_pacientes$Case.ID, df_pacientes$category, sep="_")
    
    return(df_pacientes)
    
    
  })
  
  
  NvsT_IPA <- eventReactive(input$run, {
    datapath <- input$path
    if (is.null(datapath)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    
    datapath <- str_replace_all(datapath, "\\\\", "/")
    flsall <- dir(datapath,".bam")
    flsall <- paste0(datapath, '/', df_pacientes()$File.Name)
    names(flsall) <- df_pacientes()$category
    #Genomic reference
    library("repmis")
    URL="https://github.com/RJWANGbioinfo/PAS_reference_RData/blob/master/"
    file="hg38_REF.RData"
    source_data(paste0(URL,file,"?raw=True"))
    
    refUTRraw=refUTRraw_hg38
    dfIPAraw=dfIPA_hg38
    dfLEraw=dfLE_hg38
    PASREF=REF4PAS(refUTRraw,dfIPAraw,dfLEraw)
    UTRdbraw=PASREF$UTRdbraw
    dfIPA=PASREF$dfIPA
    dfLE=PASREF$dfLE   
    dfIPA=dfIPA
    dfLE=dfLE
    IPA_OUTraw=PASEXP_IPA(dfIPA, dfLE, flsall, Strandtype="forward", SeqType="ThreeMostPairEnd")
    
    x=nrow(df_pacientes())/2
    sampleTable2 = data.frame(samplename = c(names(flsall)),
                              condition = c(rep("KD",x),rep("NT",x)))
    
    NvsT_IPA=APAdiff(sampleTable2,
                     IPA_OUTraw, 
                     conKET='NT',
                     trtKEY='KD',
                     PAS='IPA',
                     CUTreads=5,
                     p_adjust_methods="fdr")
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
    
    return(NvsT_IPA_genes_UP)
  })
  
  # NvsT_IPA_DN
  NvsT_IPA_genes_DN <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    NvsT_IPA_events_DN <- NvsT_IPA_events_DN()
    NvsT_IPA_genes_DN <- distinct(NvsT_IPA_events_DN,select=c(gene_symbol))
    
    return(NvsT_IPA_genes_DN)
  })
  
  # NvsT_IPA_NC
  NvsT_IPA_genes_NC <- eventReactive(input$run,{
    NvsT_IPA <- NvsT_IPA()
    NvsT_IPA_events_NC <- NvsT_IPA_events_NC()
    NvsT_IPA_genes_NC <- distinct(NvsT_IPA_events_NC,select=c(gene_symbol))
    
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
  })
  output$plot5 <- renderPlot({
    f()
  })
  output$plot6 <- renderPlot({
    g()
  })
  
  
  
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
      paste("plot4", ".png")
    },
    content = function(file) {
      ggsave(file, e(), width = 6000, height = 4000,units = c("px"),dpi = 300)
    }
  )
  output$download_plot5 <- downloadHandler(
    filename = function() {
      paste("plot5", ".png")
    },
    content = function(file) {
      ggsave(file, f(),width = 800, height = 700,units = c("px"),dpi = 300)
    }
  )
  output$download_plot6 <- downloadHandler(
    filename = function() {
      paste("plot6", ".png")
    },
    content = function(file) {
      ggsave(file, g(), dpi = 300)
      
    }
  )
  
  
  ##### DAPARS #####
  
  df_pacientes2 <- eventReactive(input$run2,{
    file <- input$txt_file2
    if (is.null(file)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    df_pacientes2 <- read.table(file$datapath, sep = "\t", header = TRUE)
    df_pacientes2 <- dplyr::select(df_pacientes2, File.Name,Case.ID,Sample.Type)
    df_pacientes2$File.Name <- paste("WIG/", df_pacientes2$File.Name, sep = "")
    case_order <- unique(df_pacientes2$Case.ID)
    df_pacientes2 <- df_pacientes2 %>% arrange(df_pacientes2$Sample.Type, match(df_pacientes2$Case.ID, case_order))
    normal <- c("Solid Tissue Normal")
    df_pacientes2$category <- ifelse(df_pacientes2$Sample.Type %in% normal, "Normal", "Tumor")
    df_pacientes2$category = paste(df_pacientes2$Case.ID, df_pacientes2$category, sep="_")
    df_pacientes2$File.Name<- str_replace(df_pacientes2$File.Name, ".bam", "_PDUI")
    return(df_pacientes2)
    
    
  })
  
  # DPDUI
  dpdui <- eventReactive(input$run2,{
    files <- input$txt_files
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
    
    colnames(df)[2:num_cols]<-df_pacientes2()$category
    
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
    res <- df[, grepl("_Tumor", colnames(df))] - df[, grepl("_Normal", colnames(df))]
    
    colnames(res) <- paste(colnames(df[, grepl("_Tumor", colnames(df))]),
                           colnames(df[, grepl("_Normal", colnames(df))]), sep = "-")
    
    df <-cbind(df, res)
    
    num_cols5 = ncol(df)
    
    dpdui <- df[,c(1,num_cols4:num_cols5)]
    
    return(dpdui)
  })
  
  # SHORT GENES
  short_genes <- eventReactive(input$run2,{
    dpdui <- dpdui()
    dpdui$mean <- rowMeans(dpdui[,2:length(dpdui())],na.rm=TRUE) 
    dpdui$Len_Bladder <- dpdui[,c("mean")] >= 0.2
    dpdui$Short_Bladder <- dpdui[,c("mean")] <= -0.2
    gene_len_Bladder <- dpdui[which(dpdui$Len_Bladder == 1), ]                             
    gene_short_Bladder <- dpdui[which(dpdui$Short_Bladder == 1), ]  
    
    
    return(gene_short_Bladder)
  })
  
  # LEN GENES
  len_genes <- eventReactive(input$run2,{
    dpdui <- dpdui()
    dpdui$mean <- rowMeans(dpdui[,2:length(dpdui())],na.rm=TRUE) 
    dpdui$Len_Bladder <- dpdui[,c("mean")] >= 0.2
    dpdui$Short_Bladder <- dpdui[,c("mean")] <= -0.2
    gene_len_Bladder <- dpdui[which(dpdui$Len_Bladder == 1), ]                             
    gene_short_Bladder <- dpdui[which(dpdui$Short_Bladder == 1), ]  
    
    
    return(gene_len_Bladder)
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
    file <- input$txt_file3
    if (is.null(file)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    df_pacientes3 <- read.table(file$datapath, sep = "\t", header = TRUE)
    df_pacientes3 <- dplyr::select(df_pacientes3, File.Name,Case.ID,Sample.Type)
    df_pacientes3$File.Name<- str_replace(df_pacientes3$File.Name, "rna_seq.genomic.gdc_realn.bam", "rna_seq.genomic.gdc_realn.trim.bam")
    
    case_order <- unique(df_pacientes3$Case.ID)
    df_pacientes3 <- df_pacientes3 %>% arrange(df_pacientes3$Sample.Type, match(df_pacientes3$Case.ID, case_order))
    
    #normal <- c("Solid Tissue Normal")
    df_pacientes3$category <- ifelse(df_pacientes3$Sample.Type %in% c("Solid Tissue Normal"), "Normal", "Tumor")
    df_pacientes3$category = paste(df_pacientes3$Case.ID, df_pacientes3$category, sep="_")
    
    return(df_pacientes3)
    
    
  })
  
  NvsT_APA <- eventReactive(input$run3,{
    datapath <- input$path2
    if (is.null(datapath)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    
    datapath <- str_replace_all(datapath, "\\\\", "/")
    flsall <- dir(datapath,".bam")
    flsall <- paste0(datapath, df_pacientes3()$File.Name)
    names(flsall) <- df_pacientes3()$category
    #Genomic reference
    URL="https://github.com/RJWANGbioinfo/PAS_reference_RData/blob/master/"
    file="hg38_REF.RData"
    source_data(paste0(URL,file,"?raw=True"))
    
    refUTRraw=refUTRraw_hg38
    dfIPAraw=dfIPA_hg38
    dfLEraw=dfLE_hg38
    PASREF=REF4PAS(refUTRraw,dfIPAraw,dfLEraw)
    UTRdbraw=PASREF$UTRdbraw
    dfIPA=PASREF$dfIPA
    dfLE=PASREF$dfLE   
    
    #Analysis of APA in 3’UTRs
    refUTRraw=refUTRraw
    UTRdbraw=REF3UTR(refUTRraw)
    DFUTRraw=PASEXP_3UTR(UTRdbraw, flsall, Strandtype="forward")
    
    x=nrow(df_pacientes3())/2
    sampleTable1 = data.frame(samplename = c(names(flsall)),
                              condition = c(rep("KD",x),rep("NT",x)))
    NvsT_APA=APAdiff(sampleTable1,DFUTRraw, 
                     conKET='NT',
                     trtKEY='KD',
                     PAS='3UTR',
                     CUTreads=5,
                     p_adjust_methods="fdr")
    
    NvsT_APA <- as.data.frame(NvsT_APA)
    
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
      paste("plot1", ".png")
    },
    content = function(file) {
      ggsave(file, a(), width = 6000, height = 4000,units = c("px"),dpi = 300)
    }
  )
  output$download_plot2 <- downloadHandler(
    filename = function() {
      paste("plot2", ".png")
    },
    content = function(file) {
      ggsave(file, b(),width = 800, height = 700,units = c("px"),dpi = 300)
    }
  )
  output$download_plot3 <- downloadHandler(
    filename = function() {
      paste("plot2", ".png")
    },
    content = function(file) {
      ggsave(file, d(), dpi = 300)
      
    }
  )
  
  
  ##### DGE #####
  
  df_pacientes_dge <- eventReactive(input$run_dge,{
    
    rm(list = ls(all.names = TRUE), envir = .GlobalEnv)
    
    file <- input$txt_file_dge
    if (is.null(file)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    df_pacientes_dge <- read.table(file$datapath, sep = "\t", header = TRUE)
    df_pacientes_dge <- dplyr::select(df_pacientes_dge, File.Name, Case.ID, Sample.Type)
    df_pacientes_dge$File.Name<- str_replace(df_pacientes_dge$File.Name, "rna_seq.genomic.gdc_realn.bam", "rna_seq.genomic.gdc_realn.trim.htseq.txt")
    
    df_pacientes_dge$Sample.Type <- str_replace(df_pacientes_dge$Sample.Type, "Primary Tumor", "PrimaryTumor")
    df_pacientes_dge$Sample.Type <- str_replace(df_pacientes_dge$Sample.Type, "Solid Tissue Normal", "NormalTissue")
    #normal <- c("NormalTissue")
    df_pacientes_dge$category <- ifelse(df_pacientes_dge$Sample.Type %in% c("NormalTissue"), "NormalTissue", "PrimaryTumor")
    df_pacientes_dge$category2 = paste(df_pacientes_dge$Case.ID, df_pacientes_dge$category, sep="_")
    df_pacientes_dge$category3 = paste(df_pacientes_dge$File.Name, df_pacientes_dge$category, sep="_")
    
    return(df_pacientes_dge)
  })
  
  dds <- eventReactive(input$run_dge,{
    dir <- input$path_dge
    if (is.null(dir)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    dir <- str_replace_all(dir, "\\\\", "/")
    setwd(dir)
    getwd()
    sampleFiles=grep('.htseq.txt', list.files(dir), value=TRUE)
    sampleNames=names(sampleFiles) <- df_pacientes_dge()$category3
    sampleCondition=gsub("[a-zA-Z0-9]*-[a-zA-Z0-9]*-[a-zA-Z0-9]*-[a-zA-Z0-9]*-[a-zA-Z0-9]*.rna_seq.genomic.gdc_realn.trim.htseq.txt_",'',sampleNames)
    #sampleCondition=gsub("\\d[A-Za-z].sorted.htseq.txt_",'',sampleNames)
    
    sampleTable=data.frame(sampleName = df_pacientes_dge()$File.Name, fileName = sampleFiles, condition = sampleCondition)
    reorder_idx <- match(sampleTable$sampleName, sampleTable$fileName) 
    sampleTable$fileName <- sampleTable$fileName[reorder_idx]
    sampleTable <- sampleTable[order(sampleTable$condition), ]
    
    dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                      directory = dir,
                                      design = ~ condition)
    
    
    keep <- rowSums(counts(dds)) >= 5
    dds2 <- dds[keep, ]
    dds2$condition <- factor(dds2$condition, levels = c("NormalTissue","PrimaryTumor"))
    
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
  
  res <- eventReactive(input$run_dge,{
    ddx <- ddx()
    
    res <- results(ddx, contrast=c("condition","NormalTissue","PrimaryTumor"))
    res <- as.data.frame(res)
    
    return(res)
  })
  
  normalized_counts <- eventReactive(input$run_dge,{
    ddx <- ddx()
    
    normalized_counts <- counts(ddx, normalized=TRUE)
    
    return(normalized_counts)
  })
  
  resLFC <- eventReactive(input$run_dge,{
    ddx <- ddx()
    
    resLFC <- lfcShrink(ddx, coef="condition_PrimaryTumor_vs_NormalTissue", type="apeglm")
    
    return(resLFC)
  })
  
  
  vsd <- eventReactive(input$run_dge,{
    ddx <- ddx()
    
    vst(ddx, blind=FALSE)
    
    return(vst)
  })
  
  htmap <- eventReactive(input$run_dge,{
    res <- res()
    normalized_counts <- normalized_counts()
    
    signi <- subset(res, (padj <= 0.05))
    allSig <- merge(normalized_counts, signi, by = 0)
    sigCounts <- allSig[, 2:(ncol(allSig) - 6)]
    row.names(sigCounts) <- allSig$Row.names
    
    htmap <- pheatmap(log2(sigCounts+1), scale = "row", cluster_rows=TRUE, show_rownames=FALSE, show_colnames=FALSE,
             cluster_cols=TRUE, treeheight_row = 0, treeheight_col = 50,  display_numbers=FALSE,
             color = colorRampPalette(c("blue", "white", "red"))(100), cellwidth = 12, cellheight = 0.07)
    
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
      paste("plot_pca", ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = pca_plot(), width = 5000, height = 5000, units = c("px"), bg = "white", dpi = 600)
    }
  )
  
  output$download_plot8 <- downloadHandler(
    filename = function() {
      paste("plot_volcano", ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = Volcano_dge(), width = 5000, height = 5000, units = c("px"), bg = "white", dpi = 600)
    }
  )
  
  output$download_plot9 <- downloadHandler(
    filename = function() {
      paste("plot_htmap", ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = htmap(), width = 5000, height = 5000, units = c("px"), bg = "white", dpi = 600)
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
    gene_list <- gene_list[, 2]
    
    return(gene_list)
  })
  
  
  ego_BP <- eventReactive(input$run_go,{
    de <- de()
    #pvalue_cutoff <- as.numeric(input$output_type2_go)
    
    ego <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", pvalueCutoff = Inf)
    
    ego@result$GeneRatio_num <- sapply(strsplit(as.character(ego@result$GeneRatio), "/"), function(x) as.numeric(x[1]))
    ego@result$GeneRatio_den <- sapply(strsplit(as.character(ego@result$GeneRatio), "/"), function(x) as.numeric(x[2]))
    ego@result$BgRatio_num <- sapply(strsplit(as.character(ego@result$BgRatio), "/"), function(x) as.numeric(x[1]))
    ego@result$BgRatio_den <- sapply(strsplit(as.character(ego@result$BgRatio), "/"), function(x) as.numeric(x[2]))
    ego@result$fe <- ego@result$GeneRatio_num/ego@result$GeneRatio_den / (ego@result$BgRatio_num/ego@result$BgRatio_den)
    
    return(ego)
  })
  
  ego_MF <- eventReactive(input$run_go,{
    de <- de()
    #pvalue_cutoff <- as.numeric(input$output_type2_go)
    
    ego <- enrichGO(gene = de, OrgDb = "org.Hs.eg.db", ont = "MF", keyType = "SYMBOL", pvalueCutoff = Inf)
    
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
      paste("plot_go_BP", ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = GO_plot_BP(), width = 5000, height = 5000, units = c("px"), bg = "white", dpi = 600)
    }
  )
  
  output$download_plot2_go <- downloadHandler(
    filename = function() {
      paste("plot_go_MF", ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = GO_plot_MF(), width = 5000, height = 5000, units = c("px"), bg = "white", dpi = 600)
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
      shinyalert("Warning", "Please input more than one file for intersection", type = "warning")
      return(NULL)
    } else if (nrow(files) == 2) {
      df1 <- read.csv(files[4][1, ], header = TRUE)
      df2 <- read.csv(files[4][2, ], header = TRUE)
      
      x <- list(
        dataset1 = df1[, 2],
        dataset2 = df2[, 2]
      )
      
      selected_df1 <- df1[, 2, drop = FALSE]
      selected_df2 <- df2[, 2, drop = FALSE]
      
      common_genes <- merge(selected_df1, selected_df2, by.x = 1, by.y = 1)
      common_genes <- unique(common_genes)
      
      plot <- ggvenn(x, fill_color = c("#999999", "#999555"), show_percentage = FALSE, text_size = 7)
      
      return(list(plot = plot, common_genes = common_genes))
    } else if (nrow(files) == 3) {
      df1 <- read.csv(files[4][1, ], header = TRUE)
      df2 <- read.csv(files[4][2, ], header = TRUE)
      df3 <- read.csv(files[4][3, ], header = TRUE)
      
      x <- list(
        dataset1 = df1[, 2],
        dataset2 = df2[, 2],
        dataset3 = df3[, 2]
      )
      
      selected_df1 <- df1[, 2, drop = FALSE]
      selected_df2 <- df2[, 2, drop = FALSE]
      selected_df3 <- df3[, 2, drop = FALSE]
      
      common_genes <- merge(selected_df1, selected_df2, by.x = 1, by.y = 1)
      common_genes <- merge(common_genes, selected_df3, by.x = 1, by.y = 1)
      common_genes <- unique(common_genes)
      
      plot <- ggvenn(x, fill_color = c("#999999", "#E69F00", "#56B4E9"), show_percentage = FALSE, text_size = 7)
      
      return(list(plot = plot, common_genes = common_genes))
    } else if (nrow(files) == 4) {
      df1 <- read.csv(files[4][1, ], header = TRUE)
      df2 <- read.csv(files[4][2, ], header = TRUE)
      df3 <- read.csv(files[4][3, ], header = TRUE)
      df4 <- read.csv(files[4][4, ], header = TRUE)
      
      x <- list(
        dataset1 = df1[, 2],
        dataset2 = df2[, 2],
        dataset3 = df3[, 2],
        dataset4 = df4[, 2]
      )
      
      selected_df1 <- df1[, 2, drop = FALSE]
      selected_df2 <- df2[, 2, drop = FALSE]
      selected_df3 <- df3[, 2, drop = FALSE]
      selected_df4 <- df4[, 2, drop = FALSE]
      
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
      
      df1 <- read.csv(files[4][1, ], header = TRUE)
      df2 <- read.csv(files[4][2, ], header = TRUE)
      df3 <- read.csv(files[4][3, ], header = TRUE)
      df4 <- read.csv(files[4][4, ], header = TRUE)
      df5 <- read.csv(files[4][5, ], header = TRUE)
      
      x <- list(
        dataset1 = df1[, 2],
        dataset2 = df2[, 2],
        dataset3 = df3[, 2],
        dataset4 = df4[, 2],
        dataset5 = df5[, 2]
      )
      
      selected_df1 <- df1[, 2, drop = FALSE]
      selected_df2 <- df2[, 2, drop = FALSE]
      selected_df3 <- df3[, 2, drop = FALSE]
      selected_df4 <- df4[, 2, drop = FALSE]
      selected_df5 <- df5[, 2, drop = FALSE]
      
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
      paste("plot_venn", ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = input_files()$plot, width = 5000, height = 5000, units = c("px"), bg = "white", dpi = 600)
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
  
  
  ##### SRA DOWNLOADER #####
  #SRA_project and SRA_id
  
  observeEvent(input$run_sra,{
    setwd("~/Desktop/")
    
    sra_project <- input$SRA_project
    sra_sample <- input$SRA_id
    
    if (is.null(sra_project) | is.null(sra_sample)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    
    sra_project_script <- sprintf('
    #!/bin/bash

    SRA_PROJECT="%s"

    # Search for the SRA Project ID and retrieve the corresponding SRA IDs
    sra_list=$(esearch -db sra -query $SRA_PROJECT | efetch -format runinfo | grep SRR | cut -d "," -f 1 | tr "\\n" " ")

    # Loop through the list and run fastq-dump for each SRA ID
    for sra in $sra_list; do
        echo "Downloading $sra ..."
        fastq-dump --split-3 --gzip $sra
        echo "$sra downloaded ..."
        echo ""
    done

    echo "All SRA files downloaded!"
    ', sra_project)
    
    sra_sample_script <- sprintf('
    #!/bin/bash

    SRA_SAMPLE="%s"

    echo "Downloading $SRA_SAMPLE"
    fastq-dump --split-3 --gzip $SRA_SAMPLE
    echo "$SRA_SAMPLE downloaded"
    ', sra_sample)
    
    # Execute the sample script
    system(sra_sample_script, intern = FALSE, wait = FALSE)
    
  })
  
  ##### SURVIVAL ANALYSIS #####
  
  test_df <- eventReactive(input$run_surv,{
    
    rm(list = ls(all.names = TRUE), envir = .GlobalEnv)
    
    file <- input$surv_sample_sheet
    if (is.null(file)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    df_pacientes <- read.table(file$datapath, sep = "\t", header = TRUE)
    df_pacientes <- dplyr::select(df_pacientes, File.Name, Case.ID, Sample.Type)
    df_pacientes$File.Name<- str_replace(df_pacientes$File.Name, "rna_seq.genomic.gdc_realn.bam", "rna_seq.genomic.gdc_realn.trim.htseq.txt")
    
    df_pacientes$Sample.Type <- str_replace(df_pacientes$Sample.Type, "Primary Tumor", "PrimaryTumor")
    df_pacientes$Sample.Type <- str_replace(df_pacientes$Sample.Type, "Solid Tissue Normal", "NormalTissue")
    #normal <- c("NormalTissue")
    df_pacientes$category <- ifelse(df_pacientes$Sample.Type %in% c("NormalTissue"), "NormalTissue", "PrimaryTumor")
    df_pacientes$category2 <- paste(df_pacientes$Case.ID, df_pacientes$category, sep="_")
    df_pacientes$category3 <- paste(df_pacientes$File.Name, df_pacientes$category, sep="_")
    
    file2 <- input$surv_clinical_data
    if (is.null(file2)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    clinical <- read.csv(file2$datapath, sep = "\t", header = TRUE)
    
    clinical$days_to_last_follow_up <- as.numeric(clinical$days_to_last_follow_up)
    clinical$days_to_death <- as.numeric(clinical$days_to_death)
    
    clinical$deceased <- ifelse(clinical$vital_status == "Alive", FALSE, TRUE)
    clinical$overall_survival <- ifelse(clinical$vital_status == "Alive",
                                        clinical$days_to_last_follow_up,
                                        clinical$days_to_death)
    
    dir <- input$path_surv
    if (is.null(dir)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    setwd(dir)
    getwd()
    
    sampleFiles=grep('.htseq.txt',list.files(dir), value=TRUE)
    sampleFiles
    sampleNames=names(sampleFiles)<-df_pacientes$category3
    sampleNames
    sampleCondition=gsub("[a-zA-Z0-9]*-[a-zA-Z0-9]*-[a-zA-Z0-9]*-[a-zA-Z0-9]*-[a-zA-Z0-9]*.rna_seq.genomic.gdc_realn.trim.htseq.txt_",'',sampleNames)
    #sampleCondition=gsub("\\d[A-Za-z].sorted.htseq.txt_",'',sampleNames)
    sampleCondition
    sampleTable=data.frame(sampleName=df_pacientes$File.Name, fileName=sampleFiles, condition=sampleCondition)
    sampleTable
    
    reorder_idx <- match(sampleTable$sampleName,sampleTable$fileName) 
    sampleTable$fileName <- sampleTable$fileName[reorder_idx]
    
    sampleTable <- sampleTable[order(sampleTable$condition),]
    sampleTable
    # Import sanmples data and counts into DESeqDataSet object
    # "dds" stands for Deseq-Data-Set
    dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                      directory = dir,
                                      design = ~ 1)
    
    keep <- rowSums(counts(dds)) >= 5
    dds2 <- dds[keep,]
    
    vsd <- vst(dds2, blind = FALSE)
    hn_matrix_vst <- assay(vsd)
    
    test <- hn_matrix_vst %>%
      as.data.frame() %>%
      rownames_to_column(var = "gene_id") %>%
      gather(key = "File.Name", value = "counts", -gene_id)
    
    merged_df <- merge(test, df_pacientes, by = 'File.Name', all.x = TRUE)
    
    # Replace the 'filename' column with the corresponding 'submitter_id'
    test['File.Name'] <- merged_df['Case.ID']
    test['Condition'] <- merged_df['Sample.Type']
    
    goi <- input$surv_gene
    if (is.null(goi)) {
      shinyalert("Error", "Please input all required files before running the analysis.", type = "error")
      return(NULL)
    }
    
    test_final <- test %>%
      filter(gene_id == goi) %>%
      filter(Condition == "PrimaryTumor")

    # Get median value
    median_value <- median(test_final$counts)
    
    # Denote which cases have higher or lower expression than the median
    test_final$strata <- ifelse(test_final$counts >= median_value, "HIGH", "LOW")
    
    # Add clinical data to test dataset
    test_final <- merge(test_final, clinical, by.x = "File.Name", by.y = "case_submitter_id")
    #test_final <- test_final[!duplicated(test_final$counts), ]
    test_final <- test_final[!duplicated(test_final$File.Name), ]
    
    assign("test_final", test_final, envir = .GlobalEnv)
    
    return(test_final)
  })
  
  plot_final <- eventReactive(input$run_surv,{
    test_final <- req(test_df())
    
    fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = test_final)
    plot <- ggsurvplot(fit, pval = TRUE, conf.int = FALSE, risk.table = FALSE)
    
    return(plot$plot)
  })
  
  # Use test_df() in the renderPlot
  output$plot_surv <- renderPlot({
    req(plot_final())
  }, width = 1200, height = 750)
  
  output$download_plot_surv <- downloadHandler(
    filename = function() {
      paste("plot_surv_", input$surv_gene, ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot_final(), dpi = 300)
      
    }
  )
}

# Run the app
shinyApp(ui, server)
