library(shiny)
library(plotly)
library(knitr)

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 9MB.
options(shiny.maxRequestSize = 9*1024^2)

#|--------------------------------------------------- UI CODE -------------------------------------|
ui <- fluidPage(
  titlePanel("Analyze genetic marker data with adegenet"),
  hr(),
  tabsetPanel(
    tabPanel("PCA",
             sidebarLayout(
  sidebarPanel(
    fileInput('file1', 'Choose CSV File',
              accept=c('text/csv',
                       'text/comma-separated-values,text/plain',
                       '.csv')),
    tags$hr() , 
    checkboxInput('scale', 'Scaling:', TRUE),
    tags$hr(),
    checkboxInput('center', 'Centering:', TRUE),
    tags$hr(),
    selectInput('z', 'NA method:', choices = c("mean","asis","zero"), selected = "mean"),
    tags$hr(),
    sliderInput('axes', 'Number of axes retained:', min = 3, max = 50,
                value = 3, step = 1, round = 0),
    tags$hr(),
    radioButtons('xax', 'Which PC on X axes:', c("1"="Axis1", "2"="Axis2"), "Axis1"),
    tags$hr(),
    radioButtons('yax', 'Which PC on Y axes:', c("2"="Axis2", "3"="Axis3"), "Axis2"),
    tags$hr(),
    textInput('title', "Plot Title", ''),
    
    width = 2),
  
  mainPanel(
    plotlyOutput("PCA",width = "100%", height = "150%")
    #tableOutput('table')
)
)
),
############## start pcoa panel ###############
tabPanel("PCoA",
         sidebarLayout(
           sidebarPanel(
             # fileInput('file1', 'Choose CSV File',
             #           accept=c('text/csv',
             #                    'text/comma-separated-values,text/plain',
             #                    '.csv')),
             
             checkboxInput('freq', 'Use allele frequency?', TRUE),
             tags$hr() , 
             selectInput('q', 'NA method:', choices = c("mean","asis","zero"), selected = "mean"),
             tags$hr(),
             selectInput('method', 'Distance calculation method:',
                         choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                         selected = "euclidean"),
             tags$hr(),
             sliderInput('axes', 'Number of axes retained:', min = 3, max = 50,
                         value = 3, step = 1, round = 0),
             tags$hr(),
             radioButtons('xax', 'Which PC on X axes:', c("1"="Axis1", "2"="Axis2"), "Axis1"),
             tags$hr(),
             radioButtons('yax', 'Which PC on Y axes:', c("2"="Axis2", "3"="Axis3"), "Axis2"),
             tags$hr(),
             textInput('title', "Plot Title", ''),
             
             width = 2),
           
           mainPanel(
             plotlyOutput("PCoA",width = "100%", height = "150%")
             #tableOutput('table')
           )
         )
),



tabPanel("Help",
         tags$h1("File Formatting Issues", style = "color:darkblue"),
         tags$p(" It is anticipated that most issues that arise while trying to use this application will be related
                        to file formatting. A very specific file format is required. This format has several variations, 
                depending on the type of data included in the file. The figures below illustrate acceptable file formats, 
                as they would appear in a spreadsheet program (such as Excel) before being coverted into one of the accepted
                file formats (.csv, .tsv, or .txt).", style="font-size:20px"))

),

tags$style(type="text/css",
           ".shiny-output-error { visibility: hidden; }",
           ".shiny-output-error:before { visibility: hidden; }"
)

)

          
    
#|----------------------------------------------- SERVER CODE ------------------------------------|

server <- function(input, output, session) {
  upData <- reactive({if(is.null(input$file1)) return(mtcars) 
    inFile <- input$file1
    dat <- read.genalex(inFile$datapath)
    
    #return(dat[1:5,1:5])
    
  })
  

  
  # observe({
  #   updateSelectInput(session,
  #                     "variable",
  #                     choices=names(upData()))
  #})
  
  output$PCA <- renderPlotly({
    dat.genind <- genclone2genind(upData())
    dat.X <- scaleGen(dat.genind, "scale"=input$scale, "center"=input$center, NA.method=input$z)
    dat.pca <- dudi.pca(dat.X,cent=FALSE,scale=FALSE,scannf = F,nf = input$axes)
    
    plot_ly(dat.pca$li, x=dat.pca$li[,input$xax],y=dat.pca$li[,input$yax], text=rownames(dat.pca$li),
                mode = "markers", color = pop(dat.genind), marker = list(size = 11)) %>%
    layout(p, title = input$title, 
               xaxis = list(title = "PC 1"),
              yaxis = list(title = "PC 2"))
    
  })
  
  # Display data table below plot.
  #output$table <- renderUI(upData())
  ####### PCoA #########
  
  output$PCoA <- renderPlotly({
    dat.genind <- genclone2genind(upData())
    dat.X.pco <- tab(dat.genind, "freq"=input$freq, NA.method=input$q)
    dat.pco <- dudi.pco(dist(dat.X.pco, method=input$method), scannf=FALSE, nf=input$axes)
    
    plot_ly(dat.pco$li, x=dat.pco$li[,1],y=dat.pco$li[,2], text=rownames(dat.pco$li),
            mode = "markers", color = pop(dat.genind), marker = list(size = 11)) %>%
      layout(p, title = input$title, 
             xaxis = list(title = "PC 1"),
             yaxis = list(title = "PC 2"))
    
  })
}

shinyApp(ui = ui, server = server)

