

auto_set_colorcode = function(groups, palette = "Set1") {
  #' Create color codes for groups.
  #'
  #' \code{auto_set_colorcode} returns a dataframe with group names and color
  #' codes.
  #'
  #' @param groups A vector giving the names of groups. Duplicate group names
  #'   will be disregarded.
  #' @param palette A string giving the name of palette provided in
  #'   \code{\link[RColorBrewer:brewer.pal]{RColorBrewer}}.
  #'
  #' @return A dataframe that contains color codes for samples from different
  #'   groups. The first column consists of names of the groups and the second
  #'   column contains the coresponding color codes.
  #'
  #' @examples
  #' auto_set_colorcode(c('a','b','c'))
  #'
  #' @export


  require(RColorBrewer)

  num_groups = length(unique(groups))
  if (num_groups < 10) {
    color_code = data.frame(
      Abbrev = unique(groups),
      Hex.Colors = brewer.pal(num_groups, palette),
      stringsAsFactors = FALSE
    )
  } else {
    color_code = data.frame(
      Abbrev = unique(groups),
      Hex.Colors = colorRampPalette(brewer.pal(9, palette))(num_groups),
      stringsAsFactors = FALSE
    )
  }
  return(color_code)
}


read_color_code = function(file, header = FALSE, sep = "\t") {
  #' Read color codes from a given file
  #'
  #' \code{read_color_code} reads a file with given group names and color codes
  #' and returns a dataframe with group names and color codes that can be
  #' accepted by visualization functions.
  #'
  #' The given files should contain two columns. The first column should consist
  #' of the names of the labels, and the second column should contain standard
  #' Hex color codes.
  #'
  #' @inheritParams utils::read.table
  #'
  #' @return A dataframe that contains color codes for samples from different
  #'   groups. The first column consists of names of the groups and the second
  #'   column contains the coresponding color codes.
  #'
  #' @export
  #'
  #' @examples
  #' test1 <- data.frame(c('a','b','c'), c("#f97f6c", "#6dcff6", "#fdbc4b"))
  #' tf <- tempfile()
  #' write.table(test1, tf, row.names = FALSE, col.names = FALSE, sep = "\t")
  #'
  #' check_color_code(read_color_code(tf))
  #' unlink(tf)

  color_code = read.table(
    file,
    header = header,
    sep = sep,
    as.is = TRUE,
    comment.char = ""
  )

  if (any(nchar(color_code[, 2]) != 7 |
          substr(color_code[, 2],
                 start = 1, stop = 1) != "#")) {
    stop("Invalid hex color code")
  }

  colnames(color_code) = c("Abbrev", "Hex.Colors")
  color_code[, 2] = toupper(color_code[, 2])
  return(color_code)
}


#' Check the validation of a color code dataframe
#'
#' \code{check_color_code} check the validation of a color code dataframe
#'
#' A valid color code dataframe contains two columns: \code{Abbrev}
#' and \code{Hex.Colors}. The column \code{Abbrev} contains names of
#' groups that should match the names of groups in the data. The column
#' \code{Hex.Colors} contains self defined colors of groups where colors should
#' be defined as Hex color codes. A tool for picking up colors can be found
#' \href{https://htmlcolorcodes.com/}{here}.
#'
#' @param color_code The color code dataframe
#'
#' @export
#'
#' @examples
#' temp = auto_set_colorcode(c('a','b','c'))
#' check_color_code(temp)
#'
check_color_code = function(color_code) {
  if (!is.data.frame(color_code)) {
    warning("color_code should be data frame.")
    return(FALSE)
  }

  if (!(colnames(color_code)[1] == "Abbrev" &
        colnames(color_code)[2] == "Hex.Colors")) {
    warning("Invalid column names in color_code")
    return(FALSE)
  }

  if (!all(nchar(color_code[, 2]) == 7 &
           substr(color_code[, 2],
                  start = 1, stop = 1) == "#")) {
    stop("Invalid hex color code.")
    return(FALSE)
  }

  return((
    colnames(color_code)[1] == "Abbrev" &
      colnames(color_code)[2] == "Hex.Colors"
  ) &
    all(
      nchar(color_code[, 2]) == 7 &
        substr(color_code[, 2],
               start = 1, stop = 1) == "#"
    ))
}


#' Assign color codes to samples
#'
#' \code{color_map} assigns color codes to samples based on the given color code.
#'
#' @param samples_group A vector of group names of samples
#' @param color_code A color code dataframe
#'
#' @return A vector of color codes matching the group names
#' @export
#' @examples
#' temp_code = auto_set_colorcode(c('a','b','c'))
#' temp_groups = c('a','a','b','c')
#' color_map(temp_groups, temp_code)
#'
color_map = function(samples_group, color_code) {
  # Accept a list of acronym and return the corresponding color codes
  samples_group = data.frame(id = 1:length(samples_group),
                             Abbrev = samples_group)
  res = merge(
    samples_group,
    color_code,
    by = "Abbrev",
    all = FALSE,
    sort = FALSE
  )

  return(res$Hex.Colors[order(res$id)])
}


most_freq = function(vec) {
  return(names(table(vec))[as.vector(table(vec)) == max(table(vec))][1])
}


#' Create statistics summary of nodes
#'
#' \code{stat_summery} returns statistics summaries of nodes as well as a vector
#' of Javascript codes describing the summaries.
#'
#' This function summarizes the nodes returned by
#' \code{\link[TDAmapper]{mapper}}. Basic statistics included in this function
#' are the number of samples \code{N}, the percentage of samples from different
#' groups \code{groups_per_node}, the groups dominate the nodes
#' \code{dominant_group}, the percentage of samples from the dominated group
#' \code{dominant_percent}. It also generates pie plots for nodes which are
#' saved in \code{/figures}.
#'
#' If additional information is provided, the function can provide further
#' statistics including median time at diagnosis and median time to death.
#'
#' The users can also define and add customized statistics to this function. The
#' customized statistics should be wraped into a vector of Javascript codes for
#' each of the nodes and passed to \code{add_analysis_js} following the same
#' order of nodes as in \code{obj_mapper}.
#'
#' @param obj_mapper An object of class \code{TDAmapper}.
#' @param groups_ind A vector of group names each of the samples belongs to.
#' @param dat Optional. A dataframe contains additional information of samples.
#'   If \code{add_surv_analysis==TRUE}, \code{dat} must be provided with columns
#'   \code{age_at_diagnosis} and \code{days_to_death} describing the time at
#'   diagnosis and time to death at the time of observation.
#' @param folder The name of the folder to save the generated networks.
#' @param add_surv_analysis A logical object. \code{TRUE} if median time at diagnosis
#'   and median time to death should be included.
#' @param add_analysis_js A vector of Javascript codes summerizing customized
#'   statistics.
#' @param color_code A color code dataframe.
#'
#' @return A list includes the statistics and pie plots in folder
#'   \code{/figures}. Details of elements in the list can be found in Detail.
#'
#' @export
#'
#' @examples
#' # See ?network_visualization
#'
stat_summery = function(obj_mapper,
                        groups_ind,
                        dat = NULL,
                        folder = "",
                        add_surv_analysis = FALSE,
                        add_analysis_js = NULL,
                        color_code) {
  # First get a dataframe where each column is a statistic. Then Turn them into
  # html.
  # Return: list of HTMLs
  # Also return: legend
  require(survival)

  group_by_indi = groups_ind
  n_row = length(groups_ind)

  N = c()  # number of obs
  groups = c()  # Groups in the node dand Percentage of them, string
  dominant_group = c()
  stats_sum = c()
  percent = c()
  k = 0


  if (add_surv_analysis) {
    if (is.null(dat)) {
      stop("No data with survival information provided")
    }

    if (!all(c("age_at_diagnosis", "days_to_death") %in% colnames(dat))) {
      stop("age_at_diagnosis or days_to_death not provided")
    }

    median_diagnosis = NULL
    HR = NULL
    median_death = NULL
  }

  for (vtx in obj_mapper$points_in_vertex) {
    k = k + 1
    temp_tb = sort(table(group_by_indi[vtx]), decreasing = T)
    temp_per = round(temp_tb / sum(temp_tb), 4) * 100
    temp_gp = names(temp_tb)
    groups = c(
      groups,
      paste0(
        temp_gp,
        ": <b>",
        temp_tb,
        " (",
        temp_per,
        "%)",
        "</b>",
        collapse = "<br>"
      )
    ) # Maybe display the first few

    percent = c(percent, temp_per[1])

    savepdf(paste0("pie_", k), folder)
    pie(temp_tb,
        labels = NA,
        col = color_map(temp_gp,
                        color_code = color_code))
    dev.off()

    dominant_group = c(dominant_group, temp_gp[1])

    N = c(N, paste0(
      "N : <b>",
      length(vtx),
      " (",
      round(length(vtx) / n_row, 4) * 100,
      "%)</b><br>"
    ))

    if (!is.null(dat) & add_surv_analysis) {
      med_diag_temp = median(dat[vtx, "age_at_diagnosis"],
                             na.rm = T)
      median_diagnosis = c(median_diagnosis, med_diag_temp)

      median_death = c(median_death,
                       median(as.numeric(dat[vtx, "days_to_death"]),
                              na.rm = T))
      # hr_temp_data = dat[vtx, c("vital_status",
      #                           "days_to_death",
      #                           "days_to_last_followup")]
      # hr_temp_data$Status = 0
      # hr_temp_data$Status[hr_temp_data$vital_status == "Dead"] = 1
      # hr_temp_data$Time = hr_temp_data$days_to_last_followup
      # hr_temp_data$Time[hr_temp_data$vital_status == "Dead"] =
      #   hr_temp_data$days_to_death[hr_temp_data$vital_status == "Dead"]
      # hr_temp_data$Time = as.numeric(hr_temp_data$Time)
      #
      # surv_temp = survdiff(Surv(Time, Status)~1, data = hr_temp_data)
    }
  }

  #stats_sum = paste0("N: ", N, "<br>Acronym:<br>", groups, "<br>")


  if (add_surv_analysis) {
    stats_sum = list(
      N = N,
      groups_per_node = groups,
      dominant_group = dominant_group,
      dominant_percent = percent,
      median_diagnosis = median_diagnosis,
      median_death = median_death
    )

    add_analysis_js = paste0(
      'Median diagnosis age: <b>',
      round(stats_sum$median_diagnosis, 1),
      '</b><br>',
      'Median death (days): <b>',
      round(stats_sum$median_death, 1),
      '</b><br>'
    )
  } else {
    stats_sum = list(
      N = N,
      groups_per_node = groups,
      dominant_group = dominant_group,
      dominant_percent = percent
    )
  }

  java_description = paste0(
    '<div style=\"text-align:center;\">',
    stats_sum$N,
    'Majority: <b>',
    stats_sum$dominant_group,
    ' (',
    stats_sum$dominant_percent,
    '%)</b><br>',
    add_analysis_js,
    # Mean : <b>1</b><br>
    # Variance : <b>0</b>',
    '<hr class = \"rPartvisNetwork\">
    <img style="height: 140px; width: 140px; object-fit: contain" src="',
    paste0("figures/pie_", 1:obj_mapper$num_vertices, ".png"),
    '" /><br>
    <div class =\"showOnMe2\">
    <div style=\"text-align:center;\">
    <U style=\"color:blue;\" class = \"classActivePointer\">Details</U>
    </div>
    <div class=\"showMeRpartTTp2\" style=\"display:none;\">
    ',
    stats_sum$groups_per_node ,
    '
    </script>
    <script type=\"text/javascript\">
    $(document).ready(function(){\n$(\".showOnMe2\").click(function(){\n$(\".showMeRpartTTp2\").toggle();\n$.sparkline_display_visible();\n});\n  });
    </script>
    </div>
    </div>',
    '</div>',
    "<br>"
    )

  stats_sum$java_desp = java_description
  return(stats_sum)
    }


#' Add legends to the generated graph
#'
#' @param stats_sum Results from stat_summery
#' @param color_code Color code dataframe
#'
#' @return Dataframe of legends for the generated graph
#'
legend_node = function(stats_sum = NULL, color_code) {
  #  stats_sum: results from stat_summery
  #group_list = unique(stats_sum$dominant_group)

  legend_res = data.frame(
    label = color_code$Abbrev,
    color = color_code$Hex.Colors,
    # color_map(group_list),
    shape = 'dot',
    size = 22,
    font.size = 16
  )
  return(legend_res)

}


#' Generate graphs from results of Mapper
#'
#' \code{network_visualization} generates an interactive graph from the provided
#' Mapper object.
#'
#' \code{network_visualization} generates an interactive graph based on the
#' provided Mapper object with Javascript tools from \code{visNetwork}. It
#' accepts statistics summary from the \code{\link{stat_summery}} function and
#' display them as tooltips. The tooltips can also be customized by the users by
#' passing Javascript codes with additional summerise of nodes to the argument
#' \code{add_analysis_js}.
#'
#' Nodes are colored with the colors associated with the dominated groups within
#' each of the nodes. The colors of groups can either be defined by users or by
#' function \code{\link{auto_set_colorcode}}. Self defined color codes should
#' follow the format introduced in \code{\link{check_color_code}}, and we
#' recommend reading color code files with \code{\link{read_color_code}}.
#'
#'
#' @param obj_mapper An object of class \code{TDAmapper}.
#' @param groups_ind A vector of group names each of the samples belongs to.
#' @param folder The name of the folder to save the generated networks.
#' @param dat,add_surv_analysis,add_analysis_js Arguments passed to
#'   \code{\link{stat_summery}} for customizing statistics summary of nodes.
#' @param palette A string giving the name of palette provided in
#'   \code{\link[RColorBrewer:brewer.pal]{RColorBrewer}} to automatically assign
#'   colors to nodes.
#' @param legend_ncol Number of columns of legends.
#' @param color_code The dataframe of color codes for groups of samples. If not
#'   provided, the function will automatically assign colors to different
#'   groups.
#'
#' @param color_mix Boolean. If to display the color of nodes as a mixer of the
#'   colors of samples within the nodes, where colors of samples are determined
#'   by their associated groups
#'
#' @return An HTML file and a set of pie plots will be saved under the location
#'   given in \code{folder}. The HTML file contains the interactive graph
#'   generated based on the Mapper object, and the pie plots are for the
#'   summerise of nodes.
#'
#' @export
#'
#' @examples
#' tp_data = chicken_generator(1)
#' ff = filter_coordinate(tp_data[,-1], 2)
#' tp_data_mapper = mapper.kmeans(dat = tp_data[,2:4],
#'                                filter_values = ff,
#'                                num_intervals = 10,
#'                                percent_overlap = 70)
#' network_visualization(tp_data_mapper, groups_ind = tp_data$Group, dat = tp_data[,2:4],
#'                       folder = "Exp_network")
#'
#' # Add additional analysis to nodes
#'
#' add_analysis_js = paste0('Node Index:<b>',
#'                            1:length(tp_data_mapper$points_in_vertex),
#'                            '</b><br>')
#' network_visualization(tp_data_mapper, groups_ind = tp_data$Group, dat = tp_data[,2:4],
#'                       folder = "Exp_network", add_analysis_js = add_analysis_js)
#'
#' @references Feng, T., Davila, J.I., Liu, Y., Lin, S., Huang, S. and Wang, C.,
#' 2019. Semi-supervised Topological Analysis for Elucidating Hidden Structures
#' in High-Dimensional Transcriptome Datasets. _IEEE/ACM transactions on
#' computational biology and bioinformatics._
#'
network_visualization = function(obj_mapper,
                                 groups_ind,
                                 dat = NULL,
                                 folder = "",
                                 add_surv_analysis = FALSE,
                                 add_analysis_js = NULL,
                                 palette = "Set1",
                                 legend_ncol = 2,
                                 color_code = NULL,
                                 color_mix = FALSE) {
  #  obj_mapper: the mapper object from TDAMapper
  #  groups_ind: the vector of group each individual sample belongs to


  require(visNetwork)
  require(RColorBrewer)

  if (class(obj_mapper) != "TDAmapper") {
    stop("Invalid obj_mapper.")
  }

  if (is.null(color_code)) {
    color_code = auto_set_colorcode(groups = groups_ind,
                                    palette = palette)
  } else if (!check_color_code(color_code)) {
    stop("Invalid color code.")
  }

  dir.create(file.path(folder), showWarnings = FALSE)

  obj_mapper = null_remover(obj_mapper)

  MapperNodes <- mapperVertices(obj_mapper, 1)
  MapperLinks <- mapperEdges(obj_mapper)

  members = c()
  for (i in obj_mapper$points_in_vertex) {
    members = c(members, paste0(i, collapse = ", "))
  }

  stats_sum = stat_summery(
    obj_mapper,
    groups_ind,
    dat = dat,
    folder = folder,
    add_surv_analysis = add_surv_analysis,
    add_analysis_js = add_analysis_js,
    color_code = color_code
  )

  if (!color_mix) {
    nodes <-
      data.frame(
        id = 1:nrow(MapperNodes),
        #label = MapperNodes$Nodename,
        value = MapperNodes$Nodesize,
        group = stats_sum$dominant_group,
        color = color_map(stats_sum$dominant_group,
                          color_code = color_code),
        title = stats_sum$java_desp
    )
  } else if (color_mix) {
    sample_color <- color_map(groups_ind, color_code = color_code)

    avg_color <- c()
    for (i in obj_mapper$points_in_vertex) {
      avg_color <- c(avg_color, color_mixer(sample_color[i], na.rm = TRUE))
    }

    nodes <-
      data.frame(
        id = 1:nrow(MapperNodes),
        #label = MapperNodes$Nodename,
        value = MapperNodes$Nodesize,
        group = stats_sum$dominant_group,
        color = color_map(stats_sum$dominant_group,
                          color_code = color_code),
        title = stats_sum$java_desp
      )
  } else {
    stop("Invalid color_mix. Should be Boolean.")
  }

  edges <-
    data.frame(from = MapperLinks$Linksource + 1, to = MapperLinks$Linktarget +
                 1)

  visNetwork(nodes, edges, width = "100%", height = "700px") %>%
    visInteraction(
      tooltipDelay = 500,
      selectConnectedEdges = FALSE,
      tooltipStyle = 'position: fixed;visibility:hidden;padding: 5px;
      white-space: normal;width: 150px;
      font-family: cursive;font-size:12px;font-color:purple;background-color: #E6E6E6;
      border-radius: 15px;'
    ) %>%
    visOptions(
      selectedBy = "group",
      highlightNearest = list(
        enabled = TRUE,
        degree = 2,
        hover = T
      )
    ) %>%
    visLegend(
      addNodes = legend_node(color_code = color_code),
      useGroups = FALSE,
      enabled = TRUE,
      width = 0.1,
      ncol = legend_ncol,
      position = "left"
    )  %>% sparkline::spk_add_deps() %>%
    visSave(file = "network.html", background = "white")
  save_logic = file.rename(from = "network.html", to = file.path(folder, "network.html"))


  if (save_logic) {
    cat(
      "The generated HTML file can be found in:\n",
      file.path(folder, "network.html"),
      "\n"
    )
  } else {
    warning("Cannot save file in the target folder,
            please check the working directory.")
  }
  }
