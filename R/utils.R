
#' Remove NULL nodes in the graph
#'
#' @param obj_mapper The Mapper object
#'
#' @return A Mapper object
#' @keywords internal
#'
null_remover <- function(obj_mapper) {

  null_idx <- which(sapply(obj_mapper$points_in_vertex, is.null), arr.ind = TRUE)
  if (length(null_idx) > 0) {
    obj_mapper$adjacency <- obj_mapper$adjacency[-null_idx, -null_idx]

    obj_mapper$points_in_vertex <- obj_mapper$points_in_vertex[-null_idx]

    obj_mapper$num_vertices <- length(obj_mapper$points_in_vertex)

    obj_mapper$level_of_vertex <- obj_mapper$level_of_vertex[-null_idx]
  }
  return(obj_mapper)
}


#' Spectral color map function
#'
#' \code{color_map_Spectral} maps numeric values between 0 and 1 to hex codes
#'
#' @param x A numeric vector whose entries are between 0 and 1.
#' @param name Name of the color palette. Default is Spectral.
#'
#' @return A vector of hex codes.
#' @export
#'
#' @examples
#' color_map_Spectral((1:5)/5)
#'
color_map_Spectral <- function(x, name = "Spectral") {
  require(RColorBrewer)
  color_temp <- colorRamp(brewer.pal(9, name))(x)
  color_hex <- rgb(color_temp[, 1], color_temp[, 2], color_temp[, 3],
                   maxColorValue = 255)
  return(color_hex)
}



#' Color mixer
#'
#' @param colors vector of Hex colors
#' @param weight vector of weight of each color in calculating the average
#'
#' @return A char of Hex color code.
#' @export
#' @keywords internal
#'
#' @examples
color_mixer <- function(col_vec, weight = NULL, na.rm = TRUE) {

  rgb_mat <- col2rgb(col_vec)
  if (is.null(weight)) {
    weight <- rep(1, length(col_vec))
  } else if (length(weight) != length(col_vec)) {
    stop("The length of weight is not equal")
  }

  if(na.rm) {
    col_vec <- na.omit(col_vec)
  }

  weight <- weight/sum(weight)
  avg_col <- rgb_mat %*% weight

  avg_col_hex <- rgb(red = avg_col[1],
                     green = avg_col[2],
                     blue = avg_col[3],
                     maxColorValue = 255)
  return(avg_col_hex)
}

savepdf <- function(file, folder = "", width = 10, height = 10) {
  dir.create(file.path(folder, "figures"), showWarnings = FALSE)
  fname <- paste(folder, "/figures/", file, ".png", sep = "")
  png(fname)
  par(mgp = c(0, 0, 0), tcl = -0.4, mar = c(0, 0, 0, 0), oma = c(0, 0,
                                                                 0, 0), bg = NA)

}


#' mapperVertices function
#'
#' The input to this function is a TDAmapper class object and the output
#' is a data frame of vertices that can be used as input to the networkD3
#' plot utility.
#'
#' @param m An object of class TDAmapper that is the output of the mapper
#' function.
#'
#' @return A data frame describing the vertices in the graph of the mapper
#' output and the point labels that will be displayed when the mouse
#' hovers over a vertex in the graph.
#'
#' @author Paul Pearson, \email{pearsonp@@hope.edu}
#' @references \url{https://github.com/paultpearson/TDAmapper}
#' @seealso \code{\link{mapperEdges}}
#' @keywords mapperVertices
#'
#' @keywords internal
#'
mapperVertices <- function(m, pt_labels) {

  # Hovering over vertices gives the point labels: convert the list of
  # vectors of point indices to a list of vectors of labels
  labels_in_vertex <- lapply(m$points_in_vertex, FUN = function(v) {
    pt_labels[v]
  })
  nodename <- sapply(sapply(labels_in_vertex, as.character), paste0,
                     collapse = ", ")
  nodename <- paste0("V", 1:m$num_vertices, ": ", nodename)

  # Hovering over vertices gives the point indices: list the points in
  # each vertex nodename <- sapply( sapply(m$points_in_vertex,
  # as.character), paste0, collapse=', ') concatenate the vertex number
  # with the labels for the points in each vertex nodename <- paste0('V',
  # 1:m$num_vertices, ': ', nodename )

  nodegroup <- m$level_of_vertex
  nodesize <- sapply(m$points_in_vertex, length)

  return(data.frame(Nodename = nodename, Nodegroup = nodegroup, Nodesize = nodesize))

}



#' mapperEdges function
#'
#' The input to this function is a TDAmapper class object and the output
#' is a data frame of edges that can be used as input to the networkD3
#' plot utility.
#'
#' @param m An object of class TDAmapper that is the output of the mapper function.
#'
#' @return A data frame describing the edges in the graph of the mapper output.
#'
#' @author Paul Pearson, \email{pearsonp@@hope.edu}
#' @references \url{https://github.com/paultpearson/TDAmapper}
#' @seealso \code{\link{mapperVertices}}
#' @keywords mapperEdges
#'
#' @keywords internal
#'
mapperEdges <- function(m) {
  linksource <- c()
  linktarget <- c()
  linkvalue <- c()
  k <- 1
  for (i in 2:m$num_vertices) {
    for (j in 1:(i - 1)) {
      if (m$adjacency[i, j] == 1) {
        linksource[k] <- i - 1
        linktarget[k] <- j - 1

        # Calculate the percentage of overlaped samples
        node1 <- m$points_in_vertex[[i]]
        node2 <- m$points_in_vertex[[j]]
        inter <- intersect(node1, node2)

        linkvalue[k] <- length(inter)/min(length(node1), length(node2))
        k <- k + 1
      }
    }
  }
  return(data.frame(Linksource = linksource, Linktarget = linktarget,
                    Linkvalue = linkvalue))

}


num_obs_network <- function(obj_mapper) {
  nobs <- sapply(obj_mapper$points_in_vertex, max)
  return(max(nobs))
}


#' Save results into a h5 file
#'
#' @param obj_mapper An object of class \code{TDAmapper}.
#' @param dataset The original dataset used to generate the network.
#' @param file The filename (character) of the file in which the dataset will be
#'   located.
#'
#' @return
#' @export
#'
#' @examples
save_network_h5 <- function(obj_mapper, dataset = NULL, file = "network.h5") {

  h5closeAll()

  if(class(obj_mapper) != "TDAmapper") {
    stop("The obj_mapper must be the result from mapper or
         mapper.sta of STA package.")
  } else {
    h5createFile(file)

    # Processing obj_mapper ----
    obj_mapper_list <- obj_mapper

    # Remove its class
    attr(obj_mapper_list, "class") <- NULL

    # Write to h5 file
    h5write(obj_mapper_list, file = file, name = "obj_mapper")

    # Processing raw data ----
    if(is.null(dataset)) {

      h5write("None", file = file, name = "dataset")
      h5write("None", file = file, name = "colname_feature")

    } else if(!is.data.frame(dataset)) {

      stop("Arg 'dataset' should be a data.frame.")

    } else {

      # change factors to strings
      ii <- sapply(dataset, is.factor)
      dataset[ii] <- lapply(dataset[ii], as.character)

      # Write to h5 file
      h5write(dataset, file = file, name = "dataset", DataFrameAsCompound = FALSE)

      type_col <- sapply(dataset, class)
      colnames_matrix <- data.frame(colname = names(type_col),
                                    type = type_col,
                                    stringsAsFactors = FALSE)
      # Write to h5 file
      h5write(colnames_matrix, file = file, name = "colname_feature")
    }
    h5closeAll()
  }
}


#' Load a h5 file
#'
#' Load a h5 file from save_network_h5.
#'
#' @param file The filename (character) of the file in which the dataset will be
#'   located.
#'
#' @return
#' @export
#'
#' @examples
load_network_h5 <- function(file = "network.h5") {

  require(gtools)
  h5closeAll()

  # read the obj_mapper ----
  obj_mapper <- h5read(file = file, name = "obj_mapper")

  # correct the order of list
  obj_mapper$points_in_level_set <- obj_mapper$points_in_level_set[mixedsort(names(obj_mapper$points_in_level_set))]
  obj_mapper$points_in_vertex <- obj_mapper$points_in_vertex[mixedsort(names(obj_mapper$points_in_vertex))]
  obj_mapper$vertices_in_level_set <- obj_mapper$vertices_in_level_set[mixedsort(names(obj_mapper$vertices_in_level_set))]

  # Read the colnames of features

  colname_feature <- h5read(file = file, name = "colname_feature")

  h5closeAll()

  return(list(obj_mapper = obj_mapper, colname_feature = colname_feature))
}
