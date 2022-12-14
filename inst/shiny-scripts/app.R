library(shiny)

# Define UI for data upload app ----
ui <- fluidPage(

  # App title
  titlePanel("linProt - Linear Regression for Protein Function"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(
      tags$p("Welcome to Shiny App of linProt R package."),
      br(),

      tags$p("linProt rovides a modular pipeline that allows for that seamless training and
    evaluation of linear models for protein function prediction from user provided
    labelled AA sequences. This package provides visualizations of the patterns
    learned by the linear models which can be used to generate novel proteins with
    customizable functional characteristics."),

      tags$div("The included dataset is from",
      tags$a(href = "https://www.nature.com/articles/s41598-018-33984-w#MOESM3",
       "(Karasuyama et al., 2018)"),
      " and containes a number of peptide sequences and thier experimentally
      determined peak absorbancess. These sequences were then aligned using
      the CLUSTAL OMEGA algorithm. You may choose to enter your own dataset
      below. If no dataset is provided, the program will use the one
      provided."),


      # Alignment input panel.
      fileInput("alignment", "Choose your aligned protein sequences.",
                multiple = FALSE,
                accept = c("text/fasta",
                           "text/fst")
                ),
      # Label input panel.
      fileInput("labels", "Enter the function labels for each of your
                protein sequences.",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),

      # Horizontal line ----
      tags$hr(),

      tags$p("Basic Configuration info"),

      numericInput('train_part',
                   'Percent of your dataset should be used for training,
                   remainder will be used for validation.',
                   value=75,
                   min=95,
                   max=1,
                   width='100%'
      ),

      radioButtons("encoding", "Encoding Method",
                   choices = c(onehot = "One-hot",
                               vhse = 'VHSE8'),
                   selected = "One-hot"),

      tags$hr(),

      tags$p("Enter Hyperparameters"),

      numericInput('lr',
                   'Learning Rate',
                   value=0.01,
                   min=0.0001,
                   max=1,
                   width='100%'
                   ),

      numericInput('num_iter',
                   'Number of Training Cycles',
                   value=1000,
                   min=100,
                   max=1000000,
                   width='100%'
      ),

      numericInput('l1',
                   'L1 Regularisation Constant',
                   value=0.01,
                   min=0.0001,
                   max=1,
                   width='100%'
      ),

      numericInput('l2',
                   'L2 Regularisation Constant',
                   value=0.01,
                   min=0.0001,
                   max=1,
                   width='100%'
      ),


      # Horizontal line ----
      tags$hr(),

      tags$p("Enter Visualisation Configuration"),

      numericInput('start',
                   'First amino acid index to be included in heatmap.',
                   value=1,
                   width='100%'
      ),

      numericInput('stop',
                   'Last amino acid index to be included in heatmap.',
                   value=5,
                   width='100%'
      )),

    # Main panel for displaying outputs ----
    mainPanel(

      tabsetPanel(type = "tabs",
                  tabPanel("Training Progress",
                           br(),
                           plotOutput("train_prog")
                           ),

                  tabPanel("Function Heatmap",
                           br(),
                           plotOutput("function_heatmap")
                           )
      )
    )
  )
)

# Define server logic to read selected file ----
server <- function(input, output) {
  print(input$alignment)

  if (is.null(input$alignment) | is.null(input$labels)){

  }

}
# Run the app ----
shinyApp(ui, server)