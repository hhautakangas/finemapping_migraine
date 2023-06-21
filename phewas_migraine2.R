# Load packages ----
library(shiny)
library(data.table)
library(RColorBrewer)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(DT)
library(ggplot2)
library(GetoptLong)

# load data
# for 3 color heatmap
hm <- fread("data/temp_for_phewas_candidate_variants_FGR10_heatmap_3color.txt", data.table=F)
matrix = as.matrix(hm[,2:202])
row.names(matrix) <- hm[,1]

dt <- fread("data/phewas_FGR10_meta_betas_flipped_for_shiny_plot.txt", data.table=F)
colnames(dt)[ncol(dt)] <- "category_FG"
colnames(dt)[3:4] <- c("effect allele", "other allele")
# change columns for factors, for table filtering
for( ii in c(1:5,9:11)){
  dt[,ii] <- as.factor( dt[,ii])
}


# Source helpers ----


# User interface ----
ui <- fluidPage(
  titlePanel("PheWAS: Migraine candidate variants in FinnGen R10"),
  h4("Migraine lead and/or top configuration variants for migraine risk loci."),
  
  sidebarLayout(
    #
    sidebarPanel(
      fluidRow(  
        selectInput("vars","Select locus",
                    choices = unique(dt[,2]) ),
        selectInput("cat","Select category",
                    choices = unique(dt[,"category_FG"] ))),
    ),
    
    # mainPanel(plotOutput("heatmap", click = "plot_click") )#,
    #       verbatimTextOutput("info"))
    mainPanel(
      fluidRow(
        InteractiveComplexHeatmapOutput(width1 = 700, height1 = 450, output_ui_float = TRUE,
                                        output_ui = htmlOutput("go_info"))),
      h3("Selected locus"),
      fluidRow(
        plotOutput("tile1"),
        h3("Selected category"),
        plotOutput("tile2")),
      plotOutput("heatmap"),
      fluidRow(
        div(DT::dataTableOutput("mytable"),style = "font-size: 75%; width: 75%"))# <- this line is changed
    ),
    
    
  )
)

click_action = function(df, output) {
  output[["go_info"]] = renderUI({
    if(!is.null(df)) {
      rsid = rownames(matrix)[df$row_index]
      phenoc = colnames(matrix)[df$column_index]
      
      oe = try(beta_f <- dt[ dt[,1] == rsid & dt[,"phenocode"] == phenoc, "beta_FGR10"] , silent = TRUE)
      if(inherits(oe, "try-error")) {
        beta_f = ""
      }
      oe = try(beta_m <- dt[ dt[,1] == rsid & dt[,"phenocode"] == phenoc, "beta_meta"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        beta_m = ""
      }
      oe = try(effect_allele <- dt[ dt[,1] == rsid & dt[,"phenocode"] == phenoc, "effect_allele"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        effect_allele = ""
      }
      oe = try(phenotype <- dt[ dt[,1] == rsid & dt[,"phenocode"] == phenoc, "phenotype"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        phenotype = ""
      }
      oe = try(catg <- dt[ dt[,1] == rsid & dt[,"phenocode"] == phenoc, "category_FG"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        catg = ""
      }
      oe = try(gene <- dt[ dt[,1] == rsid & dt[,"phenocode"] == phenoc, "nearest_genes"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        gene = ""
      }
      
      
      HTML(qq(
        "<div style='padding:5px 10px;border:1px solid black; width:600px; background-color:white;'>
        <h5>FinnGen R10 </h5>
        <p> <b>SNP</b>: <a href='https://r10.finngen.fi/@{rsid}' target='_blank'>@{rsid}</a></p>
        <p> <b>Phenocode</b>: <a href='https://r10.finngen.fi/pheno/@{phenoc}' target='_blank'>@{phenoc}</a> </p>
        <p> <b>Phenotype</b>: @{phenotype} </p>
      <p> <b>Category</b>: @{catg} </p>
      <p> <b>Nearest genes</b>: @{gene} </p>
      <p> <b>Effect allele</b>: @{effect_allele} </p>
      <p><b>Log odds-ratio (FG)</b>: @{beta_f}</p>
      <p><b>Log odds-ration (migraine)</b>: @{beta_m}</p>
    </div>"
      ))
    }
  })
}

# Server logic
server <- function(input, output, session){
  observe({       # <- this line is changed
    cluster_rows = FALSE
    cluster_columns = FALSE
    row_order = NULL
    column_order = NULL
    row_km = NULL
    column_km = NULL
    matrix = as.matrix(hm[,2:202])
    row.names(matrix) <- hm[,1]
    ht = Heatmap(matrix,  name = "migraine",
                 cluster_rows = cluster_rows, cluster_columns = cluster_columns,
                 row_order = row_order, column_order = column_order, 
                 row_km = row_km, column_km = column_km, show_row_names = FALSE, show_column_names = FALSE)
    
    makeInteractiveComplexHeatmap(input, output, session, ht,
                                  click_action = click_action)
  })
  
  output$tile1 <- renderPlot({
    hm2 <- dt[ dt[,2] ==  input$vars,]
    ggplot(hm2,aes(x=phenocode, y=rsid_candidate, fill=beta_FGR10))+
      geom_tile(color = "white",
                lwd = 1.5,
                linetype = 1) +
      coord_fixed()+
      theme_classic()+
      xlab("")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$tile2 <- renderPlot({
    hm2 <- dt[ dt[,"category_FG"] ==  input$cat,]
    ggplot(hm2,aes(x=phenocode, y=rsid_candidate, fill=beta_FGR10))+
      geom_tile(color = "white",
                lwd = 1.5,
                linetype = 1) +
      coord_fixed()+
      theme_classic()+
      xlab("")+
      theme(axis.text.x = element_text(angle = 45, hjust=1))
    
    
  })
  
  output$mytable = DT::renderDataTable({
    #Display table with select
    DT::datatable(dt, width = 7, rownames = FALSE, filter = 'top',
                  options = list(orderClasses = TRUE,
                                 lengthMenu = c( 10, 25, 50),
                                 pageLength = 10 ,
                                 
                                 drawCallback= JS(
                                   'function(settings) {
                                     Shiny.bindAll(this.api().table().node());}')
                  ),selection='none',escape=F)
    
    
  } )
  
  output$heatmap <- renderPlot({
    data <- hm
    cols <- colorRampPalette(brewer.pal(8,"RdBu"))(3)
    heatmap(as.matrix(data[,2:202]), Colv = NA, Rowv = NA, labRow =data$rsid, col=cols, cexRow=0.7, cexCol = 0.7, scale="none")
    
  }, height = 800, width = 800)
  
  
  
}

# Run the app
shinyApp(ui = ui, server = server)





##
# Server logic
#server <- function(input, output){

#  output$heatmap <- renderPlot({
#    data <- hm
#    cols <- colorRampPalette(brewer.pal(8,"RdBu"))(3)
#    heatmap(as.matrix(data[,2:202]), Colv = NA, Rowv = NA, labRow =data$rsid, col=cols, cexRow=0.7, cexCol = 0.7, scale="none")

#  }, height = 800, width = 800)


#  output$info <- renderPrint({
#    nearPoints(hm, input$plot_click, threshold = 10, maxpoints = 1,
#               addDist = TRUE)
#  })

#}

