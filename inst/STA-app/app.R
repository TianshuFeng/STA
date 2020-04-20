library(STA)
library(visNetwork)
library(igraph)
library(plotly)
library(rhdf5)

options(shiny.maxRequestSize = 300*1024^2)

#### User Interface ####

ui <- fluidPage(
  titlePanel("Semi-supervised topological analysis"),

  sidebarLayout(

    sidebarPanel(
      helpText("Interactive analysis of network from STA."),

      # Input: Select a network file ----
      helpText("The h5 file from function save_network_h5"),
      fileInput(inputId = "network_file",
                label = "Choose a h5 file",
                multiple = FALSE,
                accept = c(".h5")),

      # Input: Select a description file ----
      helpText("Samples in the description csv file should follow the same
               order as data used to generate the network. The
               first row is always treated as header."),
      fileInput(inputId = "descript_file",
                label = "Choose a description csv file",
                multiple = FALSE,
                accept = c(".csv")),

      # Input: Select a categorical variable ----
      conditionalPanel(condition = "input.networkType == 'cate'",
                       # Input: Select a categorical variable

                       selectInput(inputId = "discrete_feature",
                                   label = "Select a categorical feature from the original dataset:",
                                   choices = "None",
                                   selected = "None"),

                       selectInput(inputId = "group_var",
                                   label = "Select a categorical variable as the group index:",
                                   choices = "None",
                                   selected = "None"),

                       # Input: whether or not using a color mixer
                       helpText("Whether to mix the color of samples within nodes."),
                       checkboxInput(inputId = "color_mixer",
                                     label = "Color mixer")
      ),
      # Input: Select a continuous variable ----
      conditionalPanel(condition = "input.networkType == 'conti'",

                       selectInput(inputId = "continuous_feature",
                                   label = "Select a continuous feature from the original dataset:",
                                   choices = "None",
                                   selected = "None"),

                       selectInput(inputId = "continuous_var",
                                   label = "Select a continuous variable:",
                                   choices = "None",
                                   selected = "None")
      ),

      # Input: whether or not to compare nodes ----
      helpText("Whether to compare two consecutively selected nodes."),
      checkboxInput(inputId = "if_node_compare",
                    label = "Node comparison"),

      # input: Selection of brewer palettes ----
      hr(),
      selectInput(inputId = "color_palettes",
                  label = "Select a color palette:",
                  choices = c("Spectral",
                              "Blues","BuGn","BuPu","GnBu","Greens","Greys","Oranges","OrRd","PuBu","PuBuGn","PuRd","Purples","RdPu","Reds","YlGn","YlGnBu YlOrBr","YlOrRd",
                              "BrBG","PiYG","PRGn","PuOr","RdBu","RdGy","RdYlBu","RdYlGn",
                              "Set1", "Set3"),
                  selected = "Spectral"),

      # input: Selection of the degree of depth when highlighting ----
      sliderInput(inputId = "degree_network",
                  label = "Select a degree of depth when highlighting nearest nodes.",
                  min = 0, max = 9,value = 2, step = 1)

    ),

    # Main Panel ----
    mainPanel(
      tabsetPanel(tabPanel("Categorical",
                           value = 'cate'),
                  tabPanel("Continuous",
                           value = 'conti'),
                  id = "networkType"
      ),
      visNetworkOutput("network_proxy"),
      conditionalPanel(condition = "input.networkType == 'cate'",
                       "Nodes comparison:",
                       verbatimTextOutput("chisq_test_res"),
                       "",
                       "Node summary:",
                       fluidRow(style='height:200px',
                         column(5, tableOutput('summary_node_categorical')),
                         column(7, plotlyOutput('pie_chart', width = "70%"))
                       )

      ),
      conditionalPanel(condition = "input.networkType == 'conti'",
                       plotOutput('legend_continuous', height = 150),
                       tableOutput('summary_node_continuous'),
                       verbatimTextOutput("cor_res")
      )
    )
  )
)

# Node click history

node_history <- c()

#### Server ####
server <- function(input, output, session) {

  # Server: read h5 file ----
  h5_mapper <- reactive({
    req(input$network_file)
    name_list <- h5ls(input$network_file$datapath)$name

    if (!"obj_mapper" %in% name_list){
      stop("Invalid h5 file: obj_mapper not found")
    } else if (!"colname_feature" %in% name_list) {
      stop("Invalid h5 file: colname_feature not found")
    } else {
      STA:::load_network_h5(file = input$network_file$datapath)
    }
  })

  # Server: read network file ----
  obj_mapper <- reactive({

    req(input$network_file)

    h5_mapper()$obj_mapper
  })

  # Server: load colnames of features of the original dataset

  colname_feature <- reactive({
    req(input$network_file)

    h5_mapper()$colname_feature
  })

  # Server: Read description file ----
  description <- reactive({
    req(input$descript_file)

    read.csv(file = input$descript_file$datapath,
             header = TRUE)
  })

  # Server: Update the feature selection widgets
  observe({
    if(colname_feature()[1] != "None") {

      discrete_feature_name <- colname_feature()[colname_feature()[,2] == "character",1]
      continuous_feature_name <- colname_feature()[colname_feature()[,2] != "character",1]

      updateSelectInput(session = session,
                        inputId = "discrete_feature",
                        choices = c("None", discrete_feature_name),
                        selected = "None")

      updateSelectInput(session = session,
                        inputId = "continuous_feature",
                        choices = c("None", continuous_feature_name),
                        selected = "None")
    }
  })


  # Server: Update the categorical variable selection widget ----
  observe({
    req(input$descript_file)

    is.fact <- sapply(description(), is.factor)
    updateSelectInput(session = session,
                      inputId = "group_var",
                      choices = c("None", colnames(description())[is.fact]),
                      selected = "None")
  })

  # Server: Update the continuous variable selection widget ----
  observe({
    req(input$descript_file)

    is.fact <- sapply(description(), is.factor)
    updateSelectInput(session = session,
                      inputId = "continuous_var",
                      choices = c("None", colnames(description())[!is.fact]),
                      selected = "None")
  })


  # Server: Update nodes if None variable is selected ----

  observe({
    req(input$network_file)

    if((input$group_var == "None" & input$discrete_feature == "None") |
       (input$continuous_var == "None" & input$continuous_feature == "None")) {

      update_node <- data.frame(id = 1:length(obj_mapper()$points_in_vertex),
                                color = rep(color_map_Spectral(1),
                                            length(obj_mapper()$points_in_vertex))
      )

      visNetworkProxy("network_proxy") %>%
        visUpdateNodes(nodes = update_node)
    }

  })

  # Server: Update the categorical variable from description ----
  #

  groups_ind_descript <- reactive({
    req(input$network_file)

    if(input$group_var != "None" & input$networkType == 'cate') {
      description()[input$group_var][,1]
    }
  })

  # Server: Update the categorical variable from original data ----

  groups_ind_feature <- reactive({
    req(input$network_file)

    if(input$discrete_feature != "None" &
       colname_feature()[1] != "None" &
       input$networkType == 'cate') {

      feature <- h5read(file = input$network_file$datapath,
                        name = paste0("dataset/", input$discrete_feature))
      h5closeAll()

      feature
    }
  })

  groups_ind <- reactive({
    if(input$group_var != "None") {
      groups_ind_descript()
    } else if (input$discrete_feature != "None") {
      groups_ind_feature()
    } else if (input$group_var != "None" & input$discrete_feature != "None") {
      groups_ind_descript()
    }
  })

  # Server: create color code if under categorical label\
  color_code <- reactive({
    if((input$group_var != "None" | input$discrete_feature != "None") &
       input$networkType == 'cate') {

      STA:::auto_set_colorcode(groups = groups_ind(),
                               palette = input$color_palettes)
    }
  })

  # Server: Update nodes with selected categorical variable ----
  observe({
    req(input$network_file)

    if(input$networkType == 'cate' &
       (input$group_var != "None" |
        input$discrete_feature != "None")) {
      if(!input$color_mixer) {

        dom_grp <- c()
        for (i in obj_mapper()$points_in_vertex) {
          dom_grp <-
            c(dom_grp, names(sort(table(groups_ind()[i]), decreasing = T))[1])
        }

        # dom_grp <- as.numeric(as.factor(dom_grp)) - 1
        #
        # if(max(dom_grp) == 0) {
        #   dom_grp <- dom_grp + 1
        # }

        update_node <- data.frame(id = 1:length(dom_grp),
                                  color = color_map(dom_grp,
                                                    color_code = color_code())
        )

        visNetworkProxy("network_proxy") %>%
          visUpdateNodes(nodes = update_node)

      } else {
        sample_color <- color_map(groups_ind(),
                                  color_code = color_code())

        avg_color <- c()
        for (i in obj_mapper()$points_in_vertex) {
          avg_color <- c(avg_color, STA:::color_mixer(sample_color[i], na.rm = TRUE))
        }

        update_node <- data.frame(id = 1:length(avg_color),
                                  color = avg_color)

        visNetworkProxy("network_proxy") %>%
          visUpdateNodes(nodes = update_node)

      }
    }
  })



  # Server: Update the continuous variable from description ----
  #

  conti_var_descript <- reactive({
    req(input$network_file)
    req(input$descript_file)

    if(input$continuous_var != "None" &
       input$continuous_var != "" & input$networkType == 'conti') {
      description()[input$continuous_var][,1]
    }
  })

  # Server: Update the continuous variable from original data ----

  conti_var_feature <- reactive({
    req(input$network_file)

    if(input$continuous_feature != "None" &
       colname_feature()[1] != "None" &
       input$networkType == 'conti') {

      feature <- h5read(file = input$network_file$datapath,
                        name = paste0("dataset/", input$continuous_feature))
      h5closeAll()

      feature
    }
  })

  conti_var <- reactive({
    if(input$continuous_var != "None") {
      conti_var_descript()
    } else if (input$continuous_feature != "None") {
      conti_var_feature()
    } else if (input$continuous_var != "None" & input$continuous_feature != "None") {
      conti_var_descript()
    }
  })

  # Server: Update LEGEND with selected continuous variable ----

  observe({
    req(input$network_file)

    if((input$continuous_var != "None" |
        input$continuous_feature != "None") &
       input$networkType == 'conti') {

      temp_var <- conti_var()

      op <- par(mar=c(1,0,0,0))

      lth <- 50

      output$legend_continuous <- renderPlot({
        plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,1.5),xaxt="n",yaxt="n",bty="n")

        rect(
          xleft = head(seq(1, 2, length.out = lth),-1),
          ybottom = 1,
          xright =  tail(seq(1, 2, length.out = lth),-1),
          ytop = 1.5,
          col=STA:::color_map_Spectral((1:lth)/lth),
          border = NA
        )

        mtext(round(c(min(temp_var, na.rm = TRUE),
                      median(temp_var, na.rm = TRUE),
                      max(temp_var, na.rm = TRUE)), digits = 3),
              side=1,
              at=c(1.05, 1.5, 1.95),
              las=1,cex=1.2)
      })
    }
  })

  # Server: Update nodes with selected continuous variable ----

  observe({
    req(input$network_file)

    if((input$continuous_var != "None" |
        input$continuous_feature != "None") &
       input$networkType == 'conti') {

      avg_value <- c()
      for (i in obj_mapper()$points_in_vertex) {
        avg_value <-
          c(avg_value, mean(conti_var()[i], na.rm = TRUE))
      }

      # Standardize to (0, 1)
      avg_value <- (avg_value - min(avg_value))/(max(avg_value) - min(avg_value))


      update_node <- data.frame(id = 1:length(avg_value),
                                color = STA:::color_map_Spectral(avg_value,
                                                                 name = input$color_palettes))


      visNetworkProxy("network_proxy") %>%
        visUpdateNodes(nodes = update_node)
    }
  })


  # Server: When a node is selected/clicked ----

  current_node_id <- reactive({
    input$current_node_id$nodes[[1]]
  })

  # Update node clicking history, append to the beginning of the list
  observe({
    node_history <<- c(current_node_id(), node_history)

    if( length(node_history) > 10) {
      node_history <<- node_history[1:10]
    }
    print(node_history)
  })

  dist_between_nodes <- reactive({
    ad_igraph <- graph_from_adjacency_matrix(as.matrix(obj_mapper()$adjacency),
                                             mode = "undirected")
    igraph::distances(graph = ad_igraph)
  })

  observe({
    if(!is.null(current_node_id())) {

      if((input$group_var != "None" |
          input$discrete_feature != "None") & input$networkType == 'cate') {

        # If it is under a categorical variable ----
        temp_node_summary <- table(groups_ind()[obj_mapper()$points_in_vertex[[current_node_id()]]])
        # temp_node_summary <- sort(temp_node_summary, decreasing = T)

        output$summary_node_categorical <- renderTable({
          if(!is.null(current_node_id())) {
            temp_node_summary
          }
        })

        # Recent two nodes comparison
        output$chisq_test_res <- NULL
        if(length(node_history) >= 2 & !is.null(current_node_id()) & input$if_node_compare) {
          group_ind_factor <- as.factor(groups_ind())
          var_node1 <- group_ind_factor[obj_mapper()$points_in_vertex[[node_history[1]]]]
          var_node2 <- group_ind_factor[obj_mapper()$points_in_vertex[[node_history[2]]]]
          contingency <- rbind(table(var_node1), table(var_node2))

          contingency <- contingency[, colSums(contingency != 0, na.rm = T) > 0]
          chi_res <- chisq.test(contingency)

          output$chisq_test_res <- renderPrint({
            if(length(node_history) >= 2 & !is.null(current_node_id()) & input$if_node_compare) {
              print(paste("IDs of compared nodes:", node_history[1], node_history[2]))
              chi_res
            } else {
              "Two consecutively selected nodes will be compared with Chi-sq test."
            }
          })
        }

        # Plot the pie chart with plotly
        output$pie_chart <- renderPlotly({
          if(!is.null(current_node_id()) & input$networkType == 'cate'){

            plotly_color <- color_map(names(temp_node_summary),
                                      color_code = color_code())
            names(plotly_color) <- names(temp_node_summary)

            fig <- plotly::plot_ly(data = as.data.frame(temp_node_summary),
                                   labels = ~Var1,
                                   values = ~Freq,
                                   type = 'pie',
                                   marker = list(colors = plotly_color),
                                   textposition = 'inside',
                                   textinfo = 'label+percent',insidetextfont = list(color = '#FFFFFF'),
                                   showlegend = FALSE,
                                   width = 200, height = 200)
            fig <- fig %>% layout(title = 'Pie chart for the selected node',
                                  xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                  yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                  autosize = F,
                                  margin = list(l = 0, r = 0, b = 0, t = 0, pad = 0))
            fig
          } else {
            NULL
          }
        })
      } else if((input$continuous_var != "None" |
                 input$continuous_feature != "None") &
                input$networkType == 'conti') {

        # If it is under a continuous variable ----

        # Calculate the summary of the node
        output$summary_node_continuous <- renderTable({
          if(!is.null(current_node_id())) {
            temp_node_summary <- summary(conti_var()[ obj_mapper()$points_in_vertex[[current_node_id()]] ],
                                         digits = 3)
            t(as.matrix(temp_node_summary))
          }
        }, align = 'c')

        # Calculate the correlation between topological distance and average values
        avg_value <- c()

        for(points in obj_mapper()$points_in_vertex) {
          temp_values <- conti_var()[points]
          avg_value <- c(avg_value, mean(temp_values[is.finite(temp_values)],
                                         na.rm = TRUE))
        }

        same_graph_nodes <- is.finite(dist_between_nodes()[current_node_id(),])
        relative_dist <- dist_between_nodes()[current_node_id(),same_graph_nodes]
        avg_var_value <- avg_value[same_graph_nodes]
        cor_res <- cor.test(relative_dist,
                            avg_var_value)

        # If nodes are compared

        if(length(node_history) >= 2 &
           !is.null(current_node_id()) &
           input$if_node_compare) {

          conti_var1 <- conti_var()[ obj_mapper()$points_in_vertex[[node_history[1]]] ]
          conti_var2 <- conti_var()[ obj_mapper()$points_in_vertex[[node_history[2]]] ]

          ks_test_res <- ks.test(conti_var1, conti_var2)
        }


        output$cor_res <- renderPrint({
          if (length(node_history) >= 2 &
              !is.null(current_node_id()) &
              input$if_node_compare) {
            print(paste("IDs of compared nodes:", node_history[1], node_history[2]))
            print(ks_test_res)
            print("Correlation result:")
            print(cor_res)
          } else if(!is.null(current_node_id())) {
            cor_res
          } else {
            "Please select a node."
          }
        })

      } else if ( is.null(current_node_id() )) {
        output$summary_node <- NULL
      }
    }
  })

  # Server: Generate the network ----
  output$network_proxy <- renderVisNetwork({
    n_obs <- STA:::num_obs_network(obj_mapper())
    simple_visNet(obj_mapper = obj_mapper(),
                  color_filter = F,
                  groups_ind = rep(1, n_obs),
                  save_network = FALSE)%>%
      visEvents(selectNode = "function(nodes) {
        Shiny.onInputChange('current_node_id', nodes);
      ;}", deselectNode = "function(nodes) {
                Shiny.onInputChange('current_node_id', null);
                ;}")
  })

  # Server: Adjust the degree of depth of highlight when a node is selected ----
  observe({
    req(input$network_file)

    visNetworkProxy("network_proxy") %>%
      visOptions(
        highlightNearest = list(
          enabled = TRUE,
          degree = input$degree_network,
          hover = T
        ))

  })

  output$shiny_return <- renderPrint({
    input$current_node_id
  })
}

# Create Shiny app ----
shinyApp(ui, server)
