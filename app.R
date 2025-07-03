# packages ----------------------------------------------------------------
library(shiny)
library(shinyjs)
library(shinyalert)
library(shinycssloaders)
library(shinyWidgets)
library(bslib)
library(dplyr)
library(stringr)
library(plotly)
library(Hmisc)
library(RColorBrewer)
library(fpc)
library(dbscan)
library(MatrixModels)
library(GWASTools) ; data(centromeres.hg38)

source("R/corrComputation.r") # correlation functions

# ui ----------------------------------------------------------------------
ui = fluidPage(useShinyjs(),
               includeCSS("www/styles.css"), # style
               div(class = "header-bar",
                   HTML("TCM<em>viz</em> &nbsp; &nbsp; &nbsp; Mapping Gene-Gene Transcriptome Correlations in Linear (1D) Chromosomes and Hi-C-Derived (3D) Genomic Views")
               ),
               
               # Home page         
               page_navbar(nav_panel("Home",
                                     div(class="title",
                                         h1(HTML("Welcome to TCM<em>viz</em>")),
                                         h2("TCMviz is a web server designed to map Gene-Gene Transcriptome Correlations in Linear (1D) Chromosomes and Hi-C-Derived (3D) Genomic Views.")
                                     ),
                                     hr(),
                                     div(class="section",
                                         div(class = "section-title", "Overview"),
                                         p("Transcriptome Correlation Maps (TCM) are a new way of visualizing transcriptome data, going beyond simple expression heatmaps by incorporating spatial genome architecture and correlation strength."),
                                         p("TCMviz aims to investigate the complex relationships between the spatial organization of chromosomes and gene expression, with the goal of uncovering epigenetic deregulations associated with tumour progression.")
                                     ),
                                     div(class="section",
                                         div(class = "section-title", "Workflow"),
                                         #p("Insert image"),
                                         imageOutput("image") 
                                     ),
                                     div(class="section",
                                         div(class = "section-title", "Explore functions"),
                                         p(strong("Module 1:")),
                                         p("1D-TCM displays transcriptome correlations as barplots arranged by chromosome, offering a linear, genome-wide view that highlights regions of coordinated gene activity and co-expression."),
                                         p(strong("Module 2:")),
                                         p("3D-TCM extends this approach by integrating clustering and three-dimensional plotting to detect and characterize transcriptional domains, providing a more comprehensive understanding of spatial gene co-regulation across complex biological contexts.")
                                     )
               ),
               
               # 1D-TCM (intrachromosomal study)                           
               nav_panel("1D-TCM",
                         # h3("1D-TCM is ..."),
                         # hr(),
                         fluidRow(column(width = 3,
                                         div(class="sidebar-card",
                                             actionBttn("run1D",
                                                        label = "Run example",
                                                        icon = NULL,
                                                        style = "jelly",
                                                        class="custom-jelly"
                                             )
                                         )
                         ),
                         column(width = 9,
                                div(class="card",
                                    HTML("<p><b>1D-TCM</b> displays transcriptome correlations as barplots arranged by chromosome, offering a linear, genome-wide view that highlights regions of coordinated gene activity and co-expression.</p>")
                                )
                         )
                         ),
                         fluidRow(column(width = 3,
                                         div(class="sidebar-card",
                                             div(class = "section-title", fluidRow("1.Data upload", actionButton("showHelp1D", label = "?", class = "help-btn")
                                             )),
                                             fileInput("coord1D",
                                                       p('1D coordinates'),
                                                       accept = ".txt"
                                             ),
                                             fileInput("expre1D",
                                                       p('Expression matrix'),
                                                       accept = ".txt"
                                             ),
                                             div(class = "section-title", "2.Computation parameters"),
                                             numericInput("neighboursIntra",
                                                          p('Number of Close Neighbours: '),
                                                          10,
                                                          min = 1,
                                                          max = 100
                                             ),
                                             br(),
                                             actionBttn("computeIntra",
                                                        label = "Start",
                                                        icon = NULL,
                                                        style = "jelly",
                                                        class="custom-jelly"
                                             ),
                                             hr(),
                                             p("Correlated genes (.csv)"),
                                             downloadButton("downloadData1D", 
                                                            p("Download")
                                             )
                                         )
                         ),
                         column(width=9,
                                div(class="main",
                                    fluidRow(column(width=6,
                                                    uiOutput("selectChromIntra")
                                    ),
                                    column(width = 6)
                                    ),
                                    withSpinner(plotlyOutput("map1D"))
                                )
                         )
                         )
               ),
               
               # 3D-TCM (interchromosomal study)
               nav_panel("3D-TCM",
                         # h3("3D-TCM is ..."),
                         # hr(),
                         fluidRow(column(width = 3,
                                         div(class="sidebar-card",
                                             actionBttn("run3D",
                                                        label = "Run example",
                                                        icon = NULL,
                                                        style = "jelly",
                                                        class="custom-jelly"
                                             )
                                         )
                         ),
                         column(width = 9,
                                div(class="card",
                                    HTML("<p><b>3D-TCM</b> integrates clustering and three-dimensional plotting to detect and characterize transcriptional domains, providing a more comprehensive understanding of spatial gene co-regulation across complex biological contexts.</p>")
                                )
                         )
                         ),
                         fluidRow(column(width = 3,
                                         div(class="sidebar-card",
                                             div(class = "section-title", fluidRow("1.Data upload", actionButton("showHelp3D", label = "?", class = "help-btn")
                                             )),
                                             fileInput("coord", 
                                                       p('3D coordinates'),
                                                       accept = ".txt"
                                             ),
                                             fileInput("expre",
                                                       p('Expression matrix'),
                                                       accept = ".txt"
                                             ),
                                             div(class = "section-title", "2.Computation parameters"),
                                             numericInput("neighbours",
                                                          p('Number of Close Neighbours:'),
                                                          20,
                                                          min = 1,
                                                          max = 1000
                                             ),
                                             actionBttn("compute3D",
                                                        label = "Start",
                                                        icon = NULL, 
                                                        style = "jelly",
                                                        class="custom-jelly"
                                             )
                                         )
                         ),
                         column(width = 9,
                                div(class="main",
                                    fluidRow(column(width = 6,
                                                    uiOutput("selectChromInter")
                                    ),
                                    column(width = 6,
                                           uiOutput("sliderCorr")
                                    )
                                    ),
                                    withSpinner(plotlyOutput("plot_corr_color"))
                                )
                         )
                         ),
                         fluidRow(column(width = 3,
                                         div(class="sidebar-card",
                                             div(class = "section-title", "3.Clustering"),
                                             numericInput("minpoints",
                                                          p('MinPoints for DBScan: '),
                                                          11,
                                                          min = 1,
                                                          max = 100
                                             ),
                                             sliderInput("neighbour_factor",
                                                         p('Neighbour Coefficient: '),
                                                         min = 0,
                                                         max = 1,
                                                         step = 0.05,
                                                         value = 0.8
                                             ),
                                             actionBttn("computeCluster",
                                                        label = "Start",
                                                        icon = NULL,
                                                        style = "jelly",
                                                        class="custom-jelly"
                                             ),
                                             hr(),
                                             selectInput("dataset",
                                                         p("Clustering Results (.csv)"),
                                                         choices = c("clusters", "top_clusters")
                                             ),
                                             downloadButton("downloadData", 
                                                            p("Download")
                                             )
                                         )
                         ),
                         column(width = 9,
                                div(class="main",
                                    fluidRow(column(width=6,
                                                    uiOutput("selectChromClust")
                                    ),
                                    column(width = 6,
                                           uiOutput("switchHide"),
                                           uiOutput("switchView"),
                                    )
                                    ),
                                    withSpinner(plotlyOutput("plot"))
                                )
                         )
                         )
               ),
               
               # Tutorial
               nav_panel("Tutorial",
                         div(class="section",
                             div(class = "section-title", "Dataset example"),
                             p("bladder cancer from
                           NCBI GEO", tags$a(href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148079", "Series GSE148079"))
                            ),
                         div(class="section",
                             div(class = "section-title", "Import data"),
                             p("data format")
                         ),
                         div(class="section",
                             div(class = "section-title", "Data processing"),
                             p("flamingoR pour coord
                        trancriptome correlation map
                           clustering hdbscan")
                         ),
                         div(class="section",
                             div(class = "section-title", "Analysis"),
                             p("")
                         ),
                         div(class="section",
                             div(class = "section-title", "Interpretation"),
                             p("")
                         ),
               ),
               nav_spacer(),
               
               # Contact information
               # nav_panel("Contact",
               #           h4("Contact info"),
               #           h4("GitHub"),
               #           h4("Logo")
               # ),
               footer = tags$footer("How to cite:")
               )
)


# server ------------------------------------------------------------------

server = function(input, output, session) {

  # INITIALIZATION
  
  # parameters
  options(shiny.maxRequestSize = 30*1024^2)
  
  hideSpinner("plot_corr_color") 
  hideSpinner("plot")
  hideSpinner("map1D")
  
  # variables
  genePos = NULL
  genePosExp = NULL
  genePosExp_Intra = NULL
  genePosExp_Cluster = NULL
  geneExprSing = NULL
  clusters = NULL
  top_clusters = NULL
  merged = NULL
  gene_corr_csv = NULL
  
  # reactive value
  geneStartExp = reactiveVal(NULL)
  genePosExp_reactive = reactiveVal(NULL)
  genePosExp_Cluster = reactiveVal()
  
  # workflow home
  output$image = renderImage({
    list(src = "www/workflow.png", height = "100%")
    }, 
    deleteFile = FALSE) 
  
  # modalDialog 'Data input'
  observeEvent(input$showHelp1D, {
    showModal(modalDialog(title = "Data format",
                          "1D coordinate file",
                          HTML("<table border='1'>
                          <thead>
                          <tr><th>gene_id</th><th>chromosome</th><th>geneStart</th></tr>
                          </thead>
                          <tbody>
                          <tr><td>ENSG00000000003.10</td><td>chrX</td><td>99883667</td></tr>
                          <tr><td>ENSG00000000005.5</td><td>chrX</td><td>99839799</td></tr>
                          <tr><td>ENSG0000000419.8</td><td>chr20</td><td>49551404</td></tr>
                          <tr><td>ENSG0000000457.9</td><td>chr1</td><td>169818772</td></tr>
                          <tr><td>ENSG0000000460.12</td><td>chr1</td><td>169631245</td></tr>
                          <tr><td>ENSG0000000938.8</td><td>chr1</td><td>196621008</td></tr>
                          <tr><td>ENSG000001036.9</td><td>chr6</td><td>143819548</td></tr>
                          <tr><td>ENSG000001084.6</td><td>chr6</td><td>53362139</td></tr>
                          <tr><td>ENSG000001167.10</td><td>chr6</td><td>41040684</td></tr>
                          </tbody>
                               </table>"
                               ),
                          br(),
                          "Expression matrix file",
                          HTML("<table border='1'>
                          <thead>
                          <tr><th>gene_id</th><th>Tumor_T3</th><th>Tumor_T4</th><th>Tumor_T5</th><th>HT1376_1</th><th>HT1376_2</th><th>RT4_1</th><th>RT4_2</th><th>SCABER_1</th></tr>
                          </thead>
                          <tbody>
                          <tr><td>ENSG00000000003.10</td><td>21.28</td><td>3.25</td><td>6.81</td><td>16.62</td><td>12.23</td><td>57.10</td><td>111.94</td><td>22.86</td></tr>
                          <tr><td>ENSG00000000005.5</td><td>0.00</td><td>0.00</td><td>0.00</td><td>0.00</td><td>0.00</td><td>0.00</td><td>0.02</td><td>0.00</td></tr>
                          <tr><td>ENSG0000000419.8</td><td>15.06</td><td>11.42</td><td>13.44</td><td>173.26</td><td>140.47</td><td>103.60</td><td>62.17</td><td>102.00</td></tr>
                          <tr><td>ENSG0000000457.9</td><td>3.28</td><td>3.16</td><td>2.79</td><td>1.00</td><td>0.73</td><td>6.44</td><td>9.32</td><td>6.48</td></tr>
                          <tr><td>ENSG0000000460.12</td><td>2.79</td><td>3.13</td><td>4.68</td><td>22.42</td><td>13.05</td><td>6.00</td><td>8.92</td><td>19.08</td></tr>
                          <tr><td>ENSG0000000938.8</td><td>4.05</td><td>4.34</td><td>5.59</td><td>0.37</td><td>0.23</td><td>1.96</td><td>2.06</td><td>0.22</td></tr>
                          <tr><td>ENSG0000000971.11</td><td>77.40</td><td>13.08</td><td>25.96</td><td>0.21</td><td>0.18</td><td>158.63</td><td>275.09</td><td>24.72</td></tr>
                          <tr><td>ENSG0000001036.9</td><td>15.48</td><td>7.65</td><td>9.40</td><td>49.81</td><td>40.41</td><td>44.15</td><td>44.09</td><td>64.52</td></tr>
                          <tr><td>ENSG0000001084.6</td><td>27.54</td><td>13.05</td><td>4.79</td><td>5.01</td><td>3.92</td><td>157.08</td><td>129.38</td><td>30.44</td></tr>
                          <tr><td>ENSG0000001167.10</td><td>7.61</td><td>6.36</td><td>7.02</td><td>1.92</td><td>1.69</td><td>35.85</td><td>49.88</td><td>22.31</td></tr>
                          </tbody>
                               </table>"
                               ),
                          size = "xl",
                          easyClose = TRUE
                          )
              )
  })
  
  observeEvent(input$showHelp3D, {
    showModal(modalDialog(title = "Data format",
                          "3D coordinate file",
                          HTML("<table border='1'>
                          <thead>
                          <tr><th>gene_id</th><th>symbol</th><th>chromosome</th><th>x</th><th>y</th><th>z</th></tr>
                          </thead>
                          <tbody>
                          <tr><td>ENSG00000000003.10</td><td>TSPAN6</td><td>X</td><td>0.051577508</td><td>0.1147316438</td><td>0.0363846874</td></tr>
                          <tr><td>ENSG00000000005.5</td><td>TNMD</td><td>X</td><td>0.053772894</td><td>0.1213233554</td><td>0.0386941595</td></tr>
                          <tr><td>ENSG0000000419.8</td><td>DPM1</td><td>20</td><td>-0.149840092</td><td>-0.0345255021</td><td>-0.0269464080</td></tr>
                          <tr><td>ENSG0000000460.12</td><td>C1orf112</td><td>1</td><td>0.055503513</td><td>0.3094019242</td><td>-0.0575361874</td></tr>
                          <tr><td>ENSG0000000938.8</td><td>FGR</td><td>1</td><td>-0.183064122</td><td>-0.1190236858</td><td>0.1244854546</td></tr>
                          <tr><td>ENSG0000000971.11</td><td>CFH</td><td>1</td><td>0.155298115</td><td>0.1593384616</td><td>-0.0115371963</td></tr>
                          <tr><td>ENSG0000001036.9</td><td>FUCA2</td><td>6</td><td>-0.089060681</td><td>-0.1469810759</td><td>0.0454745036</td></tr>
                          <tr><td>ENSG0000001084.6</td><td>GCLC</td><td>6</td><td>0.021145584</td><td>0.0742497715</td><td>-0.0741587753</td></tr>
                          <tr><td>ENSG0000001167.10</td><td>NFYA</td><td>6</td><td>0.039597886</td><td>0.1242894061</td><td>-0.1097856644</td></tr>
                          </tbody>
                               </table>"
                               ),
                          br(),
                          "Expression matrix file",
                          HTML("<table border='1'>
                          <thead>
                          <tr><th>gene_id</th><th>Tumor_T3</th><th>Tumor_T4</th><th>Tumor_T5</th><th>HT1376_1</th><th>HT1376_2</th><th>RT4_1</th><th>RT4_2</th><th>SCABER_1</th></tr>
                          </thead>
                          <tbody>
                          <tr><td>ENSG00000000003.10</td><td>21.28</td><td>3.25</td><td>6.81</td><td>16.62</td><td>12.23</td><td>57.10</td><td>111.94</td><td>22.86</td></tr>
                          <tr><td>ENSG00000000005.5</td><td>0.00</td><td>0.00</td><td>0.00</td><td>0.00</td><td>0.00</td><td>0.00</td><td>0.02</td><td>0.00</td></tr>
                          <tr><td>ENSG0000000419.8</td><td>15.06</td><td>11.42</td><td>13.44</td><td>173.26</td><td>140.47</td><td>103.60</td><td>62.17</td><td>102.00</td></tr>
                          <tr><td>ENSG0000000457.9</td><td>3.28</td><td>3.16</td><td>2.79</td><td>1.00</td><td>0.73</td><td>6.44</td><td>9.32</td><td>6.48</td></tr>
                          <tr><td>ENSG0000000460.12</td><td>2.79</td><td>3.13</td><td>4.68</td><td>22.42</td><td>13.05</td><td>6.00</td><td>8.92</td><td>19.08</td></tr>
                          <tr><td>ENSG0000000938.8</td><td>4.05</td><td>4.34</td><td>5.59</td><td>0.37</td><td>0.23</td><td>1.96</td><td>2.06</td><td>0.22</td></tr>
                          <tr><td>ENSG0000000971.11</td><td>77.40</td><td>13.08</td><td>25.96</td><td>0.21</td><td>0.18</td><td>158.63</td><td>275.09</td><td>24.72</td></tr>
                          <tr><td>ENSG0000001036.9</td><td>15.48</td><td>7.65</td><td>9.40</td><td>49.81</td><td>40.41</td><td>44.15</td><td>44.09</td><td>64.52</td></tr>
                          <tr><td>ENSG0000001084.6</td><td>27.54</td><td>13.05</td><td>4.79</td><td>5.01</td><td>3.92</td><td>157.08</td><td>129.38</td><td>30.44</td></tr>
                          <tr><td>ENSG0000001167.10</td><td>7.61</td><td>6.36</td><td>7.02</td><td>1.92</td><td>1.69</td><td>35.85</td><td>49.88</td><td>22.31</td></tr>
                          </tbody>
                               </table>"
                               ),
                          size = "xl",
                          easyClose = TRUE
                          )
              )
  })
  
  # 1D-TCM
  observeEvent(input$computeIntra, {
    if(is.null(input$coord1D) || is.null(input$expre1D)){
      return(NULL)
    }
    req(input$coord1D, input$expre1D)
    
    showSpinner("map1D")
    
    # creation df from uploaded files
    geneExpr1D = read.table(input$expre1D$datapath, header = TRUE) # expression level
    genePos1D = read.table(input$coord1D$datapath, header = TRUE) # 1D coordinate
    
    # indexing the gene id
    row.names(genePos1D) = genePos1D[,1]
    genePos1D = genePos1D[,-1]
    
    # correlation 1D 
    geneStartExpr = transcriptomeMap1D(genePos1D, geneExpr1D, n = input$neighboursIntra)
    merged <<- merge(geneStartExpr, genePos1D, by = "row.names")
  
    # update reactive value
    geneStartExp(merged)
  })
  
  # dynamic selectInput
  output$selectChromIntra = renderUI({
    req(geneStartExp())
    selectInput(
      "choosenchrom",
      label = "Chromosome: ",
      choices = unique(geneStartExp()$chromosome.x),
      selected = unique(geneStartExp()$chromosome.x)[1]
    )
  })
  
  output$map1D = renderPlotly({
    req(geneStartExp, input$choosenchrom)
    
    # filter according to chosen chromosome
    geneStartExpChr = subset(geneStartExp(), chromosome.x == input$choosenchrom)
    
    # GWASTools :: centromeres 
    centro = centromeres.hg38
    centro$chrom = paste0("chr", centro$chrom)
    colnames(centro) = c("chromosome", "centromereStart", "centromereEnd")
    # MAD
    MyList = correlatedGenes(input$choosenchrom, geneStartExp())
    MAD = MyList$MAD
    
    gene_corr_csv <<- merged %>%
      filter(Row.names %in% MyList$CorrelatedGenes) %>%
      transmute(gene = Row.names, chr = chromosome.x, gene_start = geneStart, TC_score=corr1D) %>%
      arrange(gene_start)
    
    # plotly graph
    plot_ly(
      x = geneStartExpChr$geneStart,
      y = geneStartExpChr$corr1D,
      type = 'bar',
      name = "Gene"
    ) %>%
      add_segments(
        x = 0,
        xend = max(geneStartExpChr$geneStart),
        y = MAD,
        yend = MAD,
        name = "MAD value",
        line = list(color = "#6c25be")
      ) %>%
      add_segments(
        x = centro[centro$chromosome == input$choosenchrom, "centromereStart"],
        xend = centro[centro$chromosome == input$choosenchrom, "centromereStart"],
        y = 0,
        yend = max(geneStartExpChr$corr1D),
        name = "Centromere Start",
        line = list(color = "red")
      ) %>%
      add_segments(
        x = centro[centro$chromosome == input$choosenchrom, "centromereEnd"],
        xend = centro[centro$chromosome == input$choosenchrom, "centromereEnd"],
        y = 0,
        yend = max(geneStartExpChr$corr1D),
        name = "Centromere End",
        line = list(color = "red")
      ) %>%
      layout(
        title = "1D Transcriptome Correlation Map Visualisation",
        yaxis = list(title = "Transcription Correlation Score"),
        xaxis = list(title = "Gene Range")
      )
  })
  
  # 3D-TCM
  observeEvent(input$compute3D,{
    if(is.null(input$coord) || is.null(input$expre)){
      return(NULL)
    }
    showSpinner("plot_corr_color") # shows spinner
    
    withProgress(message = "Computation in progress", value = 0, {
      incProgress(0.1, detail = "Reading input files")
      geneExpr = read.table(input$expre$datapath[1], header = TRUE)
      genePos <<- read.table(input$coord$datapath[1])
      
      incProgress(0.4, detail = "Calculating correlation")
      genePosExpr = corrComputation(genePos, geneExpr, n = input$neighbours)
      
      incProgress(0.4, detail = "Merging data")
      genePosExp <<- merge.data.frame(genePosExpr, genePos, by = "row.names") %>%
        dplyr::mutate(chromosome = paste0("chr", as.character(chromosome)))
      write.csv2(genePosExp, "test_mat.csv", row.names = FALSE)
      
      incProgress(0.1, detail = "Finalizing")
      genePosExp_reactive(genePosExp)
      
    })
    
  })
  
  output$selectChromInter = renderUI({
    genePosExp = genePosExp_reactive()
    req(genePosExp)
    
    selectInput(
      "choosenchrom_inter",
      label = "Chromosome:",
      choices = c("All", unique(genePosExp$chromosome)),
      selected = "All"
    )
  })
  
  output$sliderCorr = renderUI({
    genePosExp = genePosExp_reactive()
    req(genePosExp)
    
    max_corr = ceiling(max(genePosExp$corr3D, na.rm = TRUE))  # get max
    
    sliderInput(
      "cor_value",
      "Minimum Correlation Value:",
      min = 0,
      max = max_corr,
      value = 0,
      step = 0.01
    )
  })
  
  output$plot_corr_color = renderPlotly({
    genePosExp = genePosExp_reactive()
    req(genePosExp, input$cor_value)
    
    genePosExpDraw = genePosExp %>%
      filter(corr3D >= input$cor_value)
    
    if (input$choosenchrom_inter != "All") {
      genePosExpDraw = genePosExpDraw %>%
        dplyr::filter(chromosome == input$choosenchrom_inter)
    }
    
    plot_ly(
      genePosExpDraw,
      x = ~x,
      y = ~y,
      z = ~z,
      marker = list(
        color = ~corr3D,
        symbol = "circle",
        size = 3,
        colorscale = "Jet",
        showscale = TRUE
      ),
      text = ~paste(
        "Chromosome ", chromosome,
        "<br>Symbol ", symbol,
        "<br>Correlation Score ", corr3D
      )
    )
  })
  
  observeEvent(input$computeCluster, {
    req(input$coord, input$expre)
    showSpinner("plot")
    
    # Read input data
    geneExpr = read.table(input$expre$datapath[1], header = TRUE)
    genePos = read.table(input$coord$datapath[1])
    
    # Prepare data for clustering
    genePosExp_clu = genePosExp  
    genePosExp_clu$cluster = 0
    
    # Define significance threshold
    MAD_value = median(genePosExp_clu$corr3D) + 2 * mad(genePosExp_clu$corr3D)
    
    genePosExp_rest = genePosExp_clu %>%
      filter(corr3D < input$neighbour_factor * MAD_value)
    
    genePosExp_clu = genePosExp_clu %>%
      filter(corr3D >= input$neighbour_factor * MAD_value)
    
    rownames(genePosExp_clu) = genePosExp_clu$Row.names
    rownames(genePosExp_rest) = genePosExp_rest$Row.names
    genePosExp_clu = genePosExp_clu[, -1]
    genePosExp_rest = genePosExp_rest[, -1]
    
    # Keep only overlapping genes
    overlapGenes = Reduce(intersect, list(rownames(geneExpr), rownames(genePos), rownames(genePosExp_clu)))
    genePos = genePos[overlapGenes, ]
    geneExpr = geneExpr[overlapGenes, ]
    
    # Clustering
    df <- as.matrix(genePosExp_clu[, 4:6])
    db <- hdbscan(df, minPts = input$minpoints)
    genePosExp_clu$cluster <- db$cluster
    
    # Filter clusters
    new_d = genePosExp_clu %>%
      count(cluster, name = "gene_in_cluster") %>%
      filter(gene_in_cluster >= input$minpoints)
    
    genePosExp_clu = filter(genePosExp_clu,
                            cluster %in% new_d$cluster
    ) %>%
      arrange(desc(cluster))
    
    # Saving in a csv file the clustering
    cluster_to_file = genePosExp_clu %>%
      filter(cluster != 0) %>%
      mutate(k_neighbours = input$neighbours,
             filtrage_facteur_mad = input$neighbour_factor,
             HDBSCAN_minpoints = input$minpoints,
             #DBSCAN_distance = input$distance_group
      ) %>%
      relocate(k_neighbours, filtrage_facteur_mad, HDBSCAN_minpoints, .before = corr3D)
    
    clusters <<- cluster_to_file
    
    data_vars = cluster_to_file %>%
      group_by(cluster) %>%
      summarise(mean_corr = mean(corr3D),
                num_pop = n(),
                n_chrom = n_distinct(chromosome)
      )
    data_fin = data_vars %>%
      group_by(cluster) %>%
      summarise(score = mean_corr * (num_pop/(max(data_vars$num_pop))) * (n_chrom/(max(data_vars$n_chrom)))
      ) %>%
      arrange(desc(score))
    
    data_top_five = data_fin[1:5, ]
    new_my_data = filter(cluster_to_file,
                         cluster %in% data_top_five$cluster
    )
    
    top_clusters <<- new_my_data
    
    # Color assignment
    genePosExp_clu = rbind(genePosExp_clu, genePosExp_rest) %>%
      mutate(opacity = 1.0, size = 50, color = "grey")
    
    all_colors = unlist(mapply(brewer.pal, brewer.pal.info$maxcolors, rownames(brewer.pal.info)))
    color_palette = sample(all_colors, length(unique(genePosExp_clu$cluster)))
    
    color_map = list()
    for (cluster_id in unique(genePosExp_clu$cluster)) {
      if (cluster_id != 0) {
        color_map[[as.character(cluster_id)]] <- color_palette[[length(color_map) + 1]]
      }
    }
    
    genePosExp_clu$color = ifelse(
      genePosExp_clu$cluster == 0,
      "grey",
      unlist(color_map[as.character(genePosExp_clu$cluster)])
    )
    
    # Store in reactiveVal
    genePosExp_Cluster(genePosExp_clu)
    
    
    output$switchHide = renderUI({
      switchInput(
        inputId = "checkBoxHideGenes",
        label = "Hide non-clustered Genes", 
        #labelWidth = "100%",
        value = TRUE
      )
    })
    
    output$switchView = renderUI({
      switchInput(
        inputId = "checkBoxViewLevels",
        label = "View by correlation level", 
        #labelWidth = "100%"
      )
    })
  })
  
  # Dynamic chromosome selector
  output$selectChromClust = renderUI({
    df = genePosExp_Cluster()
    req(df)
    
    selectInput(
      "choosenchrom_cluster",
      "Chromosome:",
      choices = c("All", unique(df$chromosome)),
      selected = "All"
    )
  })
  
  # Download csv
  data <- reactive({
    get(input$dataset) # chosen file
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(input$dataset, ".csv")
    },
    content = function(file) {
      write.csv2(data(), file)
    }
  )
  
  data2 <- reactive({
    gene_corr_csv
  })
  
  output$downloadData1D <- downloadHandler(
    filename = function() {
      paste0("correlatedGenes_", input$choosenchrom, ".csv")
    },
    content = function(file) {
      write.csv2(data2(), file)
    }
  )
  
  # 3D cluster plot
  output$plot = renderPlotly({
    df = genePosExp_Cluster()
    req(df)
    
    showSpinner("plot")
    
    
    if (input$checkBoxHideGenes) {
      df = df %>% filter(cluster != 0)
    }
    
    if (input$choosenchrom_cluster != "All") {
      chrom = str_remove(input$choosenchrom_cluster, "Chr")
      df = df %>% filter(chromosome == chrom)
    }
    
    if (input$checkBoxViewLevels) {
      plot = plot_ly(
        df,
        x = ~x, y = ~y, z = ~z,
        size = ~size,
        opacity = ~opacity,
        marker = list(
          color = ~corr3D,
          symbol = "circle",
          size = ~size,
          colorscale = "Jet",
          showscale = TRUE,
          opacity = ~opacity
        ),
        text = ~paste(
          "Chromosome ", chromosome,
          "<br>Symbol ", symbol,
          "<br>Correlation Score ", corr3D,
          "<br>Cluster ", cluster
        )
      )
    } else {
      plot = plot_ly(
        df,
        x = ~x, y = ~y, z = ~z,
        size = ~size,
        opacity = ~opacity,
        marker = list(
          color = ~color,
          symbol = "circle",
          line = list(width = 0)
        ),
        text = ~paste(
          "Cluster ", cluster,
          "<br>corr ", corr3D,
          "<br>Chromosome ", chromosome,
          "<br>Symbol ", symbol
        )
      )
    }
    
    return(plot)
  })
  
  # Run example Button
  observeEvent(input$run1D, {
    print("example10k")
    
    showSpinner("map1D")
    
    # creation df from uploaded files
    geneExpr1D = read.table("hic_bladder_expression_10k.txt", header = TRUE) # expression level
    genePos1D = read.table("hic_bladder_1D_coords_10k.txt", header = TRUE) # 1D coordinate
    
    # indexing the gene id
    row.names(genePos1D) = genePos1D[,1]
    genePos1D = genePos1D[,-1]
    
    # correlation 1D 
    # geneStartExpr = transcriptomeMap1D(genePos1D, geneExpr1D, n = input$neighboursIntra)
    # merged = merge(geneStartExpr, genePos1D, by = "row.names")
    
    merged <<- read.csv2("test_barplot.csv")
    
    # update reactive value
    geneStartExp(merged)
    
    
    
    output$map1D = renderPlotly({
      req(geneStartExp(), input$choosenchrom)
      
      # filter according to chosen chromosome
      geneStartExpChr = subset(geneStartExp(), chromosome.x == input$choosenchrom)
      
      # GWASTools :: centromeres 
      centro = centromeres.hg38
      centro$chrom = paste0("chr", centro$chrom)
      colnames(centro) = c("chromosome", "centromereStart", "centromereEnd")
      
      # MAD
      MyList = correlatedGenes(input$choosenchrom, geneStartExp())
      MAD = MyList$MAD
      
      print(head(merged))
      gene_corr_csv <<- merged %>%
        filter(Row.names %in% MyList$CorrelatedGenes) %>%
        transmute(gene = Row.names, chr = chromosome.x, gene_start = geneStart, TC_score=corr1D) %>%
        arrange(gene_start)
      
      data2 <- reactive({
        gene_corr_csv
      })
      
      output$downloadData1D <- downloadHandler(
        filename = function() {
          paste0("correlatedGenes_", input$choosenchrom, ".csv")
        },
        content = function(file) {
          write.csv2(data2(), file)
        }
      )
      # plotly graph
      plot_ly(
        x = geneStartExpChr$geneStart,
        y = geneStartExpChr$corr1D,
        type = 'bar',
        name = "Gene"
      ) %>%
        add_segments(
          x = 0,
          xend = max(geneStartExpChr$geneStart),
          y = MAD,
          yend = MAD,
          name = "MAD value",
          line = list(color = "#6c25be")
        ) %>%
        add_segments(
          x = centro[centro$chromosome == input$choosenchrom, "centromereStart"],
          xend = centro[centro$chromosome == input$choosenchrom, "centromereStart"],
          y = 0,
          yend = max(geneStartExpChr$corr1D),
          name = "Centromere Start",
          line = list(color = "red")
        ) %>%
        add_segments(
          x = centro[centro$chromosome == input$choosenchrom, "centromereEnd"],
          xend = centro[centro$chromosome == input$choosenchrom, "centromereEnd"],
          y = 0,
          yend = max(geneStartExpChr$corr1D),
          name = "Centromere End",
          line = list(color = "red")
        ) %>%
        layout(
          title = "1D Transcriptome Correlation Map Visualisation",
          yaxis = list(title = "Transcription Correlation Score"),
          xaxis = list(title = "Gene Range")
        )
    })
  })
  
  
  # Run example button
  observeEvent(input$run3D, {
    
    print("example10k")
    
    showSpinner("plot_corr_color") # shows spinner
    showSpinner("plot")
    
    withProgress(message = "Computation in progress", value = 0, {
      incProgress(0.1, detail = "Reading input files")
      # geneExpr = read.table("hic_bladder_expression_10k.txt", header = TRUE)
      # genePos <<- read.table("scabber_3D_coords_10k.txt")
      
      incProgress(0.4, detail = "Calculating correlation")
      # genePosExpr = corrComputation(genePos, geneExpr, n = input$neighbours)
      
      incProgress(0.4, detail = "Merging data")
      # genePosExp <<- merge.data.frame(genePosExpr, genePos, by = "row.names") %>%
      #   dplyr::mutate(chromosome = paste0("chr", as.character(chromosome)))
      
      genePosExp <<- read.csv2('test_mat.csv', )
      print(head(genePosExp))
      incProgress(0.1, detail = "Finalizing")
      genePosExp_reactive(genePosExp)
      
    })
    
    output$selectChromInter = renderUI({
      genePosExp = genePosExp_reactive()
      req(genePosExp)
      
      selectInput(
        "choosenchrom_inter",
        label = "Chromosome:",
        choices = c("All", unique(genePosExp$chromosome)),
        selected = "All"
      )
    })
    
    output$sliderCorr = renderUI({
      genePosExp = genePosExp_reactive()
      req(genePosExp)
      
      max_corr = ceiling(max(genePosExp$corr3D, na.rm = TRUE))  # safely get max
      
      sliderInput(
        "cor_value",
        "Minimum Correlation Value:",
        min = 0,
        max = max_corr,
        value = 0,
        step = 0.01
      )
    })
    
    output$plot_corr_color = renderPlotly({
      genePosExp = genePosExp_reactive()
      req(genePosExp, input$cor_value)
      
      genePosExpDraw = genePosExp %>%
        filter(corr3D >= input$cor_value)
      
      if (input$choosenchrom_inter != "All") {
        genePosExpDraw = genePosExpDraw %>%
          dplyr::filter(chromosome == input$choosenchrom_inter)
      }
      
      plot_ly(
        genePosExpDraw,
        x = ~x,
        y = ~y,
        z = ~z,
        marker = list(
          color = ~corr3D,
          symbol = "circle",
          size = 3,
          colorscale = "Jet",
          showscale = TRUE
        ),
        text = ~paste(
          "Chromosome ", chromosome,
          "<br>Symbol ", symbol,
          "<br>Correlation Score ", corr3D
        )
      )
    })
    
    
    showSpinner("plot")
    
    # read input data
    geneExpr = read.table("hic_bladder_expression_10k.txt", header = TRUE)
    genePos = read.table("scabber_3D_coords_10k.txt")
    
    # prepare data for clustering
    genePosExp_clu = genePosExp  
    genePosExp_clu$cluster = 0
    
    # define significance threshold
    MAD_value = median(genePosExp_clu$corr3D) + 2 * mad(genePosExp_clu$corr3D)
    
    genePosExp_rest = genePosExp_clu %>%
      filter(corr3D < input$neighbour_factor * MAD_value)
    
    genePosExp_clu = genePosExp_clu %>%
      filter(corr3D >= input$neighbour_factor * MAD_value)
    
    rownames(genePosExp_clu) = genePosExp_clu$Row.names
    rownames(genePosExp_rest) = genePosExp_rest$Row.names
    genePosExp_clu = genePosExp_clu[, -1]
    genePosExp_rest = genePosExp_rest[, -1]
    
    # keep only overlapping genes
    # overlapGenes = Reduce(intersect, list(rownames(geneExpr), rownames(genePos), rownames(genePosExp_clu)))
    # genePos = genePos[overlapGenes, ]
    # geneExpr = geneExpr[overlapGenes, ]
    
    # clustering
    df = as.matrix(genePosExp_clu[, 4:6])
    db = hdbscan(df, minPts = input$minpoints)
    genePosExp_clu$cluster = db$cluster
    
    # filter clusters
    new_d <- genePosExp_clu %>%
      count(cluster, name = "gene_in_cluster") %>%
      filter(gene_in_cluster >= input$minpoints) 
    
    genePosExp_clu <- filter(genePosExp_clu, cluster %in% new_d$cluster) %>%
      arrange(desc(cluster))
    
    # saving in a csv file the clustering
    cluster_to_file = genePosExp_clu %>%
      filter(cluster != 0) %>% 		  
      mutate(k_neighbours = input$neighbours,
             filtrage_facteur_mad = input$neighbour_factor,
             HDBSCAN_minpoints = input$minpoints,
             #DBSCAN_distance = input$distance_group
      ) %>%
      relocate(k_neighbours, filtrage_facteur_mad, HDBSCAN_minpoints, .before = corr3D)
    
    clusters <<- cluster_to_file
    
    data_vars = cluster_to_file %>%
      group_by(cluster) %>%
      summarise(mean_corr = mean(corr3D), 
                num_pop = n(), 
                n_chrom = n_distinct(chromosome)
      )
    data_fin = data_vars %>%
      group_by(cluster) %>%
      summarise(score = mean_corr * (num_pop/(max(data_vars$num_pop))) * (n_chrom/(max(data_vars$n_chrom)))
      ) %>% 
      arrange(desc(score))
    
    data_top_five = data_fin[1:5, ] 
    new_my_data = filter(cluster_to_file,
                         cluster %in% data_top_five$cluster
    )
    
    top_clusters <<- new_my_data
    
    # color assignment
    genePosExp_clu <- rbind(genePosExp_clu, genePosExp_rest) %>%
      mutate(opacity = 1.0, size = 50, color = "grey")
    
    all_colors <- unlist(mapply(brewer.pal, brewer.pal.info$maxcolors, rownames(brewer.pal.info)))
    color_palette <- sample(all_colors, length(unique(genePosExp_clu$cluster)))
    
    color_map <- list()
    for (cluster_id in unique(genePosExp_clu$cluster)) {
      if (cluster_id != 0) {
        color_map[[as.character(cluster_id)]] <- color_palette[[length(color_map) + 1]]
      }
    }
    
    genePosExp_clu$color <- ifelse(
      genePosExp_clu$cluster == 0,
      "grey",
      unlist(color_map[as.character(genePosExp_clu$cluster)])
    )
    
    # store in reactiveVal
    genePosExp_Cluster(genePosExp_clu)
    
    
    output$switchHide = renderUI({
      switchInput(
        inputId = "checkBoxHideGenes",
        label = "Hide non-clustered Genes",
        labelWidth = "100%",
        value = TRUE
      )
    })
    
    output$switchView = renderUI({
      switchInput(
        inputId = "checkBoxViewLevels",
        label = "View by correlation level",
        labelWidth = "100%"
      )
    })
    
    output$selectChromClust = renderUI({
      df <- genePosExp_Cluster()
      req(df)
      
      selectInput(
        "choosenchrom_cluster",
        "Chromosome:",
        choices = c("All", unique(df$chromosome)),
        selected = "All"
      )
    })
    
    
    # 3D cluster plot
    output$plot = renderPlotly({
      df <- genePosExp_Cluster()
      req(df)
      
      # Download csv
      data <- reactive({
        get(input$dataset) # chosen file
      })
      
      output$downloadData <- downloadHandler(
        filename = function() {
          paste0(input$dataset, ".csv")
        },
        content = function(file) {
          write.csv2(data(), file)
        }
      )  
      
      if (input$checkBoxHideGenes) {
        df = df %>% filter(cluster != 0)
      }
      
      if (input$choosenchrom_cluster != "All") {
        chrom = str_remove(input$choosenchrom_cluster, "Chr")
        df = df %>% filter(chromosome == chrom)
      }
      
      if (input$checkBoxViewLevels) {
        plot <- plot_ly(
          df,
          x = ~x, y = ~y, z = ~z,
          size = ~size,
          opacity = ~opacity,
          marker = list(
            color = ~corr3D,
            symbol = "circle",
            size = ~size,
            colorscale = "Jet",
            showscale = TRUE,
            opacity = ~opacity
          ),
          text = ~paste(
            "Chromosome ", chromosome,
            "<br>Symbol ", symbol,
            "<br>Correlation Score ", corr3D,
            "<br>Cluster ", cluster
          )
        )
      } else {
        plot = plot_ly(
          df,
          x = ~x, y = ~y, z = ~z,
          size = ~size,
          opacity = ~opacity,
          marker = list(
            color = ~color,
            symbol = "circle",
            line = list(width = 0)
          ),
          text = ~paste(
            "Cluster ", cluster,
            "<br>corr ", corr3D,
            "<br>Chromosome ", chromosome,
            "<br>Symbol ", symbol
          )
        )
      }
      
      return(plot)
    })
  })
  
}


# shinyApp ----------------------------------------------------------------

shinyApp(ui, server)