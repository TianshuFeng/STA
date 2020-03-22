library(STA)
library(visNetwork)

#### User Interface ####

ui <- fluidPage(
  titlePanel("Semi-supervised topological analysis"),

  sidebarLayout(

    sidebarPanel(
      helpText("Interactive analysis of network from STA."),

      # Input: Select a network file ----
      helpText("The network rds file containing a TDAmapper class object"),
      fileInput(inputId = "network_file",
                label = "Choose a network RDS file",
                multiple = FALSE,
                accept = c(".rds")),

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
                       selectInput(inputId = "group_var",
                                   label = "Select a categorical variable as the group index:",
                                   choices = ""),
                       # Input: whether or not using a color mixer
                       checkboxInput(inputId = "color_mixer",
                                     label = "Color mixer")
                       ),
      # Input: Select a continuous variable ----
      conditionalPanel(condition = "input.networkType == 'conti'",
                       selectInput(inputId = "continuous_var",
                                   label = "Select a continuous variable:",
                                   choices = "")
                       ),

      # input: Selection of brewer palettes ----
      hr(),
      selectInput(inputId = "color_palettes",
                  label = "Select a color palette:",
                  choices = c("Spectral",
                              "Blues","BuGn","BuPu","GnBu","Greens","Greys","Oranges","OrRd","PuBu","PuBuGn","PuRd","Purples","RdPu","Reds","YlGn","YlGnBu YlOrBr","YlOrRd",
                              "BrBG","PiYG","PRGn","PuOr","RdBu","RdGy","RdYlBu","RdYlGn",
                              "Set1", "Set3"),
                  selected = "Spectral")

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
      verbatimTextOutput("shiny_return")
    )
  )
)


#### Server ####
server <- function(input, output, session) {

  # Server: read network file ----
  obj_mapper <- reactive({

    req(input$network_file)

    readRDS(input$network_file$datapath)
    })

  # Server: Read description file ----
  description <- reactive({
    req(input$descript_file)

    read.csv(file = input$descript_file$datapath,
             header = TRUE)
  })

  # Server: Update the categorical variable selection widget ----
  observe({
    req(input$descript_file)

    is.fact <- sapply(description(), is.factor)
    updateSelectInput(session = session,
                      inputId = "group_var",
                      choices = colnames(description())[is.fact])
  })

  # Server: Update the continuous variable selection widget ----
  observe({
    req(input$descript_file)

    is.fact <- sapply(description(), is.factor)
    updateSelectInput(session = session,
                      inputId = "continuous_var",
                      choices = colnames(description())[!is.fact])
  })

  # Server: Update nodes with selected categorical variable ----
  observe({
    req(input$network_file)
    req(input$descript_file)

    if(input$group_var != "" & input$networkType == 'cate') {
      groups_ind <- description()[input$group_var]

      if(!input$color_mixer) {
        dom_grp <- c()
        for (i in obj_mapper()$points_in_vertex) {
          dom_grp <-
            c(dom_grp, names(sort(table(groups_ind[i,]), decreasing = T))[1])
        }

        dom_grp <- as.numeric(as.factor(dom_grp)) - 1

        if(max(dom_grp) == 0) {
          dom_grp <- dom_grp + 1
        }

        update_node <- data.frame(id = 1:length(dom_grp),
                                  color = STA:::color_map_Spectral(dom_grp / max(dom_grp),
                                                                   name = input$color_palettes))

        visNetworkProxy("network_proxy") %>%
          visUpdateNodes(nodes = update_node)

      } else {
        sample_color <- as.numeric(as.factor(groups_ind[,1])) - 1

        if(max(sample_color) == 0) {
          sample_color <- sample_color + 1
        }

        sample_color <- STA:::color_map_Spectral(sample_color/max(sample_color),
                                                 name = input$color_palettes)

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

  # Server: Update nodes with selected continuous variable ----

  observe({
    req(input$network_file)
    req(input$descript_file)

    if(input$continuous_var != "" & input$networkType == 'conti') {
      conti_var <- description()[input$continuous_var]

      avg_value <- c()
      for (i in obj_mapper()$points_in_vertex) {
        avg_value <-
          c(avg_value, mean(conti_var[i,], na.rm = TRUE))
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

  # Server: Generate the network ----
  output$network_proxy <- renderVisNetwork({
    n_obs <- STA:::num_obs_network(obj_mapper())
    simple_visNet(obj_mapper = obj_mapper(),
                  color_filter = F,
                  groups_ind = rep(1, n_obs),
                  save_network = FALSE)%>%
      visEvents(selectNode = "function(nodes) {
        Shiny.onInputChange('current_node_id', nodes);
      ;}")
  })

  output$shiny_return <- renderPrint({
    input$current_node_id$nodes[[1]]
  })
}

# Create Shiny app ----
shinyApp(ui, server)
