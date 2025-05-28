library(shiny)
library(bslib)
library(GWASTools)
data(centromeres.hg38)
library(shinyjs)
#library(shinythemes)
library(shinyalert)
library(shinycssloaders)
library(dplyr)
library(stringr)
library(plotly)
source("R/corrComputation.r") # correlation functions

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .header-bar {
      font-family: 'Times New Roman';
        background-color: #D9E6F2;
        padding: 15px 30px;
        font-size: 20px;
        font-weight: bold;
        color: #1F3B73;
      }
      .navbar {
        background-color: #74A9DC !important;
        font-family: 'Impact'!important;
      }
      .navbar-default .navbar-nav > li > a {
        color: black !important;
        
      }
      .navbar-default .navbar-nav > .active > a {
        background-color: white !important;
        color: black !important;
        font-weight: bold;
      }
      h1 {
        text-align: center;
        font-size: 28px;
        margin-top: 30px;
        font-weight: bold;
      }
      h2 {
      text-align: center;
        font-size: 20px;
        font-weight: normal;
      }
      .section-title {
        font-weight: bold;
        margin-top: 10px;
        font-size: 18px;
      }
      .sidebar-card {
        background-color: #fff;
        padding: 20px;
        border: 1px solid #ddd;
        border-radius: 10px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.05);
        height: 100%;
      }
    "))
  ),
  
  #header above navbar
  div(class = "header-bar",
      HTML("TCM<em>viz</em> &nbsp; &nbsp; &nbsp; Mapping Gene-Gene Transcriptome Correlations in Linear (1D) Chromosomes and Hi-C-Derived (3D) Genomic Views")
  ),
  
  # Navbar
  page_navbar(
    #title = "",#HTML("TCM<em>viz</em>"),  # Title with no space, italic 'viz'
    
    nav_panel("Home",
             h1(HTML("Welcome to TCM<em>viz</em>")),
             h2("TCMviz is a web server designed to map Gene-Gene Transcriptome Correlations in Linear (1D) Chromosomes and Hi-C-Derived (3D) Genomic Views."),
             hr(),
             div(style = "padding: 20px;",
                 div(class = "section-title", "Overview"),
                 p(HTML("Transcriptome Correlation Maps (TCM) are a new way of visualizing transcriptome data, going beyond simple expression heatmaps by incorporating spatial genome architecture and correlation strength.")),
                 p(HTML("TCMviz aims to investigate the complex relationships between the spatial organization of chromosomes and gene expression, with the goal of uncovering epigenetic deregulations associated with tumour progression.")),
             ),
             div(style = "padding: 20px;",
                 div(class = "section-title", "Explore functions"),
                 p(HTML("Module 1: 1D-TCM")),
                 p(HTML("Module 2: 3D-TCM")),
             ),
    ),
    
    nav_panel("1D-TCM",
              h4("1D-TCM is ..."),
              hr(),
              fluidRow(
              column(width = 4,
                     div(class="sidebar-card",
                actionButton("example", "Run example"),
                p("OR"),
                p("Data upload"),
                fileInput("coord1D", "Upload 1D Coordinates file as .txt",accept = ".txt"),
                fileInput("expre1D", "Upload Expression file as .txt", accept = ".txt"),
                numericInput("neighboursIntra", "Number of Close Neighbours:", 10, min = 1, max = 100),
                actionButton("computeIntra","Start", icon = icon("redo"))
                     )
              ),
             column(width=8,
                    fluidRow(column(width=6,
                      selectInput("choosenchrom",
                                "Chromosome:",
                                c("chr1","chr2","chr3","chr4","chr5",
                                  "chr6","chr7","chr8","chr9","chr10","chr11",
                                  "chr12","chr13","chr14","chr15","chr16","chr17",
                                  "chr18","chr19","chr20","chr21","chr22","chrX","chrY")
                                )
                    ),
                    column(width = 6,
                    sliderInput("cor_value_intra", "Minimum Correlation Value:", min = 0, max = 20, value = 0),)),
                    withSpinner(plotlyOutput("map1D")),
              )
                
              
              )
             ),
    nav_panel("3D-TCM",
              h4("3D-TCM is ..."),
              hr(),
              fluidRow(
                column(width = 4,
                       div(class="sidebar-card",
                           actionButton("example2", "Run example"),
                           hr(),
                           h4("Correlation"),
                           p("Data upload"),
                           fileInput("coord", "Upload Coordinates file as .txt",accept = ".txt"),
                           fileInput("expre", "Upload Expression Matrix file as .txt",accept = ".txt"),
                           numericInput("neighbours", "Number of Close Neighbours:", 20, min = 1, max = 1000),
                           actionButton("compute3D", "Start", icon = icon("redo")),
                           hr(),
                           div(class="clust",
                               h4("Clustering"),
                               numericInput("minpoints", "MinPoints for DBScan:", 15, min = 1, max = 100),
                               numericInput("distance_group", "Distance of Circle:", 0.02, min = 0.001, max = 1, step = 0.001),
                               sliderInput("neighbour_factor", "Neighbour Coefficient:", min = 0, max = 1, step  =0.05, value = 0.8),
                               numericInput("neighbours2", "Number of Close Neighbours:", 20, min = 1, max = 1000),
                               actionButton("computeCluster", "Start", icon = icon("redo")),
                               checkboxInput("checkBoxHideGenes", "Hide non-clustered Genes", value = FALSE),
                               checkboxInput("checkBoxViewLevels", "View by correlation level", value = FALSE),
                               )
                       ),
                      
                ),
                column(width=8,
                       fluidRow(column(width=6,
                                       selectInput("choosenchrom_inter",
                                                   "Chromosome:",
                                                   c("All","Chr1","Chr2","Chr3","Chr4","Chr5",
                                                     "Chr6","Chr7","Chr8","Chr9","Chr10","Chr11",
                                                     "Chr12","Chr13","Chr14","Chr15","Chr16","Chr17",
                                                     "Chr18","Chr19","Chr20","Chr21","Chr22","ChrX","ChrY"
                                                   )
                                       ),
                       ),
                       column(width = 6,
                              sliderInput("cor_value", "Minimum Correlation Value:", min = 0, max = 100, value = 0),)),
                       # sliderInput("x_coord_cut_intra", "X coordinates Cut:", min = -1.0, max = 1.0, value = c(-1.0, 1.0), step = 0.1),
                       # sliderInput("y_coord_cut_intra", "Y coordinates Cut:", min = -1.0, max = 1.0, value = c(-1.0, 1.0),step = 0.1),
                       # sliderInput("z_coord_cut_intra", "Z coordinates Cut:", min = -1.0, max = 1.0, value = c(-1.0, 1.0),step = 0.1),
                       withSpinner(plotlyOutput("plot_corr_color")),
                       selectInput("choosenchrom_cluster",
                                   "Chromosome:",
                                   c("All","Chr1","Chr2","Chr3","Chr4","Chr5",
                                     "Chr6","Chr7","Chr8","Chr9","Chr10","Chr11",
                                     "Chr12","Chr13","Chr14","Chr15","Chr16","Chr17",
                                     "Chr18","Chr19","Chr20","Chr21","Chr22","ChrX","ChrY"
                                   )
                       ),
                       withSpinner(plotlyOutput("plot"))
                )
                
                
              )
              ),
    nav_panel("Tutorial",
              h4("Dataset example"),
              h4("Import data"),
              h4("Data processing"),
              h4("Analysis"),
              h4("Interpretation")),
    nav_spacer(),
    nav_panel("Contact",
              h4("Contact info"),
              h4("GitHub"),
              h4("Logo")
              ),
    footer = "How to cite:"
  )
)

server <- function(input, output, session) {
  output$map1D <- renderPlotly({
    if(is.null(input$coord1D) || is.null(input$expre1D)){
      return(NULL)
    }
    # GWASTools :: centromeres 
    centro = centromeres.hg38
    centro$chrom = paste0("chr",
                          centro$chrom
    )
    colnames(centro) = c("chromosome",
                         "centromereStart",
                         "centromereEnd"
    )
    
    
    # creation df from uploaded files
    geneExpr1D = read.table(input$expre1D$datapath[1], header = TRUE) # expression level
    genePos1D = read.table(input$coord1D$datapath[1], header = TRUE) # 1D coordinate
    
    # indexing the gene id
    row.names(genePos1D) = genePos1D[,1]
    genePos1D = genePos1D[,-1]
    
    # correlation 1D
    geneStartExpr = transcriptomeMap1D(genePos1D, geneExpr1D, n=input$neighboursIntra)
    geneStartExp = merge.data.frame(geneStartExpr, genePos1D, by = "row.names")
    
    
    
    MyList <- correlatedGenes(input$choosenchrom, geneStartExp)
    MAD <- MyList$MAD
    print(MAD)
    
    geneStartExpChr <-subset(geneStartExp,
                             geneStartExp$chromosome.x == input$choosenchrom)
    print(head(geneStartExpChr))
    
    # plotly graph		
    plot <- plot_ly(x=geneStartExpChr$geneStart,
                    y=geneStartExpChr$corr1D,
                    width = 900,
                    height = 350,
                    type = 'bar',
                    name = "Correlation Value"
    ) %>%
      add_segments(x = 0,
                   xend = max(geneStartExpChr$geneStart),
                   y = MAD,
                   yend = MAD, 
                   name = "MAD value", 
                   line=list(color="#6c25be")
      ) %>%
      # centromere's location
      add_segments(x = centro[which(centro$chromosome== input$choosenchrom),"centromereStart"],
                   xend= centro[which(centro$chromosome== input$choosenchrom),"centromereStart"],
                   y = 0,
                   yend = max(geneStartExpChr$corr1D),
                   line=list(color="Red"),
                   name = "Centromere Start"
      ) %>%
      add_segments(x = centro[which(centro$chromosome== input$choosenchrom),"centromereEnd"],
                   xend= centro[which(centro$chromosome== input$choosenchrom),"centromereEnd"],
                   y = 0, yend = max(geneStartExpChr$corr1D), 
                   line=list(color="Red"),
                   name = "Centromere End"
      ) %>%
      layout(title = "Intrachromosomal Study 1D Visualisation",
             yaxis= list(title="Correlation Value"),
             xaxis=list(title="Gene Range")
      )
    
    return(plot)
    
  })
  dataframe = observeEvent(input$compute3D,{
    if(is.null(input$coord) || is.null(input$expre)){
      return(NULL)
    }
    withProgress(message = "Computation in progress", value = 0, {
      incProgress(0.1, detail = "Reading input files")
      geneExpr = read.table(input$expre$datapath[1], header = TRUE)
      genePos <<- read.table(input$coord$datapath[1])
      
      incProgress(0.4, detail = "Calculating correlation")
      genePosExpr = corrComputation(genePos, geneExpr, n = input$neighbours)
      
      incProgress(0.4, detail = "Merging data")
      genePosExp <<- merge.data.frame(genePosExpr, genePos, by = "row.names")
      
      incProgress(0.1, detail = "Finalizing")
    })
    
  })
  
  draw = observeEvent(input$compute3D,{ 
    output$plot_corr_color <- renderPlotly({
      if(is.null(genePosExp)){
        return(NULL)
      }
      
      showSpinner("plot_corr_color") # shows spinner
      
      # filter according to the user choices
      genePosExpDraw <- genePosExp %>%
        filter(corr3D >= input$cor_value)
      #          x >= input$x_coord_cut[1],
      #          x <= input$x_coord_cut[2],
      #          y >= input$y_coord_cut[1],
      #          y <= input$y_coord_cut[2],
      #          z >= input$z_coord_cut[1],
      #          z <= input$z_coord_cut[2]
      #   )
      # print(head(genePosExpDraw))
      
      if(input$choosenchrom_inter != "All"){
        chrom = str_remove(input$choosenchrom_inter, "Chr") 
        genePosExpDraw <- genePosExpDraw %>%
          filter(chromosome == chrom) 
      }
      
      # plotly graph
      plot = plot_ly(genePosExpDraw, 
                     x = genePosExpDraw$x, 
                     y = genePosExpDraw$y, 
                     z = genePosExpDraw$z,
                     # color scale 
                     marker = list(color = genePosExpDraw$corr3D, 
                                   symbol = "circle",
                                   size = 3,
                                   colorscale = "Jet",
                                   showscale = TRUE
                     ),
                     # label
                     text = ~paste("Chromosome ",
                                   chromosome,
                                   "<br>Symbol",
                                   symbol,
                                   "<br>Correlation Score",
                                   corr3D
                     )
      )
      gc() # for the report on memory usage
      return(plot)
    }
    )
    
  })
  
  ataframe_cluster = observeEvent(input$computeCluster,{
    if(is.null(input$coord) || is.null(input$expre)){
      return(NULL)
    }
    geneExpr = read.table(input$expre$datapath[1], header = TRUE) 
    genePos = read.table(input$coord$datapath[1])
    
    # calculation correlation between expression level and position
    genePosExpr = corrComputation(genePos, geneExpr, n = input$neighbours2)
    
    # merging by gene id, "Row.names" column created consequently 
    genePosExp <<- merge.data.frame(genePosExpr, genePos, by = "row.names")
    
    # initialization of a cluster column in the merged df
    genePosExp_clu = genePosExp %>%
      mutate(cluster = 0)
    
    # definition of a significance level
    MAD_value = median(genePosExp_clu$corr3D) + 2 * (mad(genePosExp_clu$corr3D))
    genePosExp_rest = genePosExp_clu %>%
      filter(corr3D < (input$neighbour_factor * MAD_value)) # below the threshold
    genePosExp_clu = genePosExp_clu %>%
      filter(corr3D >= (input$neighbour_factor * MAD_value)) # above the threshold
    
    # deleting the "Row.names" column by indexing it 
    rownames(genePosExp_clu) = genePosExp_clu$Row.names 
    rownames(genePosExp_rest) = genePosExp_rest$Row.names
    genePosExp_clu = genePosExp_clu[,-1] 
    genePosExp_rest = genePosExp_rest[,-1]
    
    
    # keeping only overlapping genes
    overlapGenes <- Reduce(intersect, 
                           list(rownames(geneExpr),
                                rownames(genePos),
                                rownames(genePosExp_clu)
                           )
    ) 
    genePos <- genePos[overlapGenes,] 
    geneExpr <- geneExpr[overlapGenes,] 
    
    # density based clustering (dbscan)
    df <- as.matrix(genePosExp_clu[,4:6]) 
    print(head(df)) # "x", y", "z" columns
    db <- fpc::dbscan(df,
                      eps = input$distance_group,
                      MinPts = input$minpoints
    ) 
    genePosExp_clu$cluster = c(db$cluster) # assigning a gene to a cluster
    
    print("finished")
    
    # filtering the clusters
    new_d = genePosExp_clu %>%
      count(cluster, name = "gene_in_cluster") %>%
      filter(gene_in_cluster >= input$minpoints) 
    
    genePosExp_clu = filter(genePosExp_clu, 
                            cluster %in% new_d$cluster)
    
    genePosExp_clu <- genePosExp_clu %>%
      arrange(desc(cluster))
    # coloration of cluster
    genePosExp_clu = rbind(genePosExp_clu, genePosExp_rest) %>%
      mutate(opacity = 1.0,
             size = 50,
             color = "grey"
      )
    all_colors = unlist(mapply(brewer.pal,
                               brewer.pal.info$maxcolors,
                               rownames(brewer.pal.info)
    )
    )
    color_palette = sample(all_colors,
                           length(unique(genePosExp_clu$cluster))
    )
    print(head(genePosExp_clu))
    
    color_chooser = 0
    for(gene in rownames(genePosExp_clu)){
      if(genePosExp_clu[gene,]$cluster != 0){ # cluster 0 = grey
        temp = filter(genePosExp_clu,
                      cluster == genePosExp_clu[gene,]$cluster
        ) # focus on one cluster
        if(temp[1,]$color == "grey"){ # no color yet for the cluster
          color_chooser = color_chooser + 1 
          genePosExp_clu[gene,]$color = color_palette[color_chooser]
        }else{
          genePosExp_clu[gene,]$color = temp[1,]$color # already a color
        }
      }
    }
    
    print("pass")
    genePosExp_Cluster <<- genePosExp_clu 
  })
  draw_cluster = observeEvent(input$computeCluster,{
    output$plot <- renderPlotly({
      if(is.null(genePosExp_Cluster)){
        return(NULL)
      }
      showSpinner("plot")
      print(input$checkBoxHideGenes)
      
      # filtering according to the user choices
      if(input$checkBoxHideGenes){
        genePosExp_Cluster = genePosExp_Cluster %>%
          filter(cluster != 0)
      }
      if(input$choosenchrom_cluster !="All"){
        chrom = str_remove(input$choosenchrom_cluster, "Chr")
        print(chrom)
        genePosExp_Cluster<- genePosExp_Cluster %>%
          filter(chromosome == chrom)
      }
      
      # plotly graph
      if(input$checkBoxViewLevels){ # correlation + cluster
        plot = plot_ly(genePosExp_Cluster,
                       x = genePosExp_Cluster$x,
                       y = genePosExp_Cluster$y,
                       z = genePosExp_Cluster$z,
                       size = genePosExp_Cluster$size,
                       opacity = genePosExp_Cluster$opacity,
                       marker = list(color = genePosExp_Cluster$corr3D,
                                     symbol = "circle",
                                     size = genePosExp_Cluster$size,
                                     colorscale = "Jet",
                                     showscale = TRUE,
                                     opacity = genePosExp_Cluster$opacity
                       ),
                       # label
                       text = ~paste("Chromosome ",
                                     chromosome, 
                                     "<br>Symbol", 
                                     symbol, 
                                     "<br>Correlation Score", 
                                     corr3D, 
                                     "<br> Cluster ", 
                                     cluster
                       )
        )
      }else{
        plot = plot_ly(genePosExp_Cluster, # cluster
                       x = genePosExp_Cluster$x, 
                       y = genePosExp_Cluster$y, 
                       z = genePosExp_Cluster$z, 
                       size = genePosExp_Cluster$size,
                       opacity = genePosExp_Cluster$opacity,
                       marker = list(color = genePosExp_Cluster$color,
                                     symbol = "circle",
                                     line = list(width = 0)
                       ),
                       text = ~paste("Cluster ",
                                     cluster,
                                     "<br> corr",
                                     corr3D,
                                     "<br> Chromosome ",
                                     chromosome,
                                     "<br> symbol",
                                     symbol
                       )
        )
      }
      
      
      return(plot)
      
    })
  })
  
}

shinyApp(ui, server)
