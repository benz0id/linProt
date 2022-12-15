# library(shiny)
# library(seqinr)
# library(assertthat)

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
                   choices = c(one_hot = "one-hot",
                               vhse = 'VHSE8'),
                   selected = "one-hot"),

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
      ),

      actionButton("run", "Run Analysis")


      ),

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
                           ),

                  tabPanel("Maximum and Minimum Sequences",
                           br(),
                           textOutput("max_seq"),
                           textOutput("min_seq")
                  ),

      )
    )
  )
)

# Define server logic to read selected file ----
server <- function(input, output) {

  observeEvent(input$run, {

    # Check that user has input required files.
    if (is.null(input$alignment) | is.null(input$labels)){
      # No externally provided data. Use local data.
      examples_path <- '../extdata/peptide_alignment.fasta'
      labels_path <- '../extdata/lambda_max.csv'
    } else {
      # Use externally provided data
      examples_path <- input$labels$datapath
      labels_path <- input$alignment$datapath
    }

    # Safely parse external data.
    tryCatch(
      {
        labels <- as.numeric(unlist(read.csv(labels_path, header = FALSE)$V1))
        examples <- read.fasta(examples_path)
      },
      error = function(e) {
        stop(safeError(e))
      }
    )

    # Convert to list of strings.
    examples_list <- list()
    for (i in seq_along(examples)){
      examples_list <- append(examples_list, paste(examples[[i]], sep='',
                                                   collapse = ''))
    }
    examples <- examples_list

    # Extract run parameters.
    num_train = as.integer(length(examples) * input$train_part / 100)

    if (input$encoding == "one-hot"){
      encode_fxn = encode_onehot
    } else {
      encode_fxn = encode_physchem
    }

    reg_hypers <- c(l1 = input$l1, l2 = input$l2)

    # Shuffle, partition, and encode data set.
    shuffled_datasets <- shuffled_partitions(examples, labels, num_train,
                                             encode=encode_fxn)

    # Designate test and train data.
    train_data <- shuffled_datasets$e1
    train_labels <- shuffled_datasets$l1
    valid_data <- shuffled_datasets$e2
    valid_labels <- shuffled_datasets$l2



    # Train a linear model to perform regression.
    model <- linear_train(train_data,
                          train_labels,
                          valid_data,
                          valid_labels,
                          reg = 'elastic',
                          reg_hypers = reg_hypers,
                          num_iter = input$num_iter,
                          rec_loss_every = 10,
                          learning_rate = input$lr)

    # View the expected influence of each residue on the function, (lambda max
    # in this case).

    output$train_prog <- renderPlot({plot_cost_over_rep(model)})

    output$function_heatmap <- renderPlot({residue_effect_heatmap(model,
                                                                  input$start,
                                                                    input$stop)})
    text1 <- paste(c('Maximum Sequence: ',
                    maximal_sequence(model), collapse=''))
    output$max_seq <- renderText({text1})


    text2 <- paste(c('Minimal Sequence: ', maximal_sequence(model, do_min=TRUE),
                    collapse=''))
    output$min_seq <- renderText({text2})
  })
}



# Run the app ----
shinyApp(ui, server)
