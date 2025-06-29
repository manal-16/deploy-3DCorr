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
               includeCSS("www/styles.css"),
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
                                         p("Insert image")
                                         ),
                                     div(class="section",
                                         div(class = "section-title", "Explore functions"),
                                         p("Module 1: 1D-TCM"),
                                         p("Module 2: 3D-TCM")
                                         )
                                     ),
               
               # 1D-TCM (intrachromosomal study)                           
               nav_panel("1D-TCM",
                         h3("1D-TCM is ..."),
                         hr(),
                         fluidRow(column(width = 3,
                                         div(class="sidebar-card",
                                             actionBttn("run1D",
                                                        label = "Run example",
                                                        icon = NULL,
                                                        style = "jelly",
                                                        class="custom-jelly"
                                                        )
                                             )
                                         )
                                  ),
                         fluidRow(column(width = 3,
                                  div(class="sidebar-card",
                                      div(class = "section-title", fluidRow("1.Data upload", actionButton("showHelp", label = "?", class = "help-btn")
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
                                             h3("3D-TCM is ...")
                                             )
                                         )
                                  ),
                         fluidRow(column(width = 3,
                                         div(class="sidebar-card",
                                             div(class = "section-title", fluidRow("1.Data upload", actionButton("showHelp", label = "?", class = "help-btn")
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
                                                          15,
                                                          min = 1,
                                                          max = 100
                                                          ),
                                             numericInput("distance_group",
                                                          p('Distance of Circle: '),
                                                          0.02,
                                                          min = 0.001,
                                                          max = 1,
                                                          step = 0.001
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
                         h4("Dataset example"),
                         h4("Import data"),
                         h4("Data processing"),
                         h4("Analysis"),
                         h4("Interpretation")
               ),
               nav_spacer(),
               
               # Contact information
               nav_panel("Contact",
                         h4("Contact info"),
                         h4("GitHub"),
                         h4("Logo")
               ),
               footer = tags$footer("How to cite:")
               )
)


# server ------------------------------------------------------------------

server = function(input, output, session) {
  
  # INITIALIZATION
  
  # parameters
  options(shiny.maxRequestSize = 30*1024^2)
  shinyjs::disable("checkBoxViewLevels")
  
  hideSpinner("plot_corr_color") # hidden until "draw3D" button pressed
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
  merged=NULL
  gene_corr_csv=NULL
  
  # reactive value
  geneStartExp = reactiveVal(NULL)
  genePosExp_reactive = reactiveVal(NULL)
  genePosExp_Cluster = reactiveVal()
  
  observeEvent(input$showHelp, {
           showModal(modalDialog(
                title = "Data format",
                "...",
                size = "l",
                easyClose = TRUE
          ))
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
    
    #write.csv2(merged, "test_barplot.csv", row.names = FALSE)
    
    
    
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
      name = "Correlation Value"
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
        yaxis = list(title = "Correlation Value"),
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
    db <- fpc::dbscan(df, eps = input$distance_group, MinPts = input$minpoints)
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
             DBSCAN_minpoints = input$minpoints,
             DBSCAN_distance = input$distance_group) %>%
      relocate(k_neighbours, filtrage_facteur_mad, DBSCAN_minpoints, DBSCAN_distance, .before = corr3D)
    
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
    
    top_clusters <<- my_new_data
    
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
        name = "Correlation Value"
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
          yaxis = list(title = "Correlation Value"),
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
    db = fpc::dbscan(df, eps = input$distance_group, MinPts = input$minpoints)
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
             DBSCAN_minpoints = input$minpoints,
             DBSCAN_distance = input$distance_group) %>%
      relocate(k_neighbours, filtrage_facteur_mad, DBSCAN_minpoints, DBSCAN_distance, .before = corr3D)
    
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