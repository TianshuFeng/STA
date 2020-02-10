
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
#'
#' @return A vector of hex codes.
#' @export
#'
#' @examples
#' color_map_Spectral((1:5)/5)
#'
color_map_Spectral <- function(x) {
  require(RColorBrewer)
  color_temp <- colorRamp(brewer.pal(11, "Spectral"))(x)
  color_hex <- rgb(color_temp[, 1], color_temp[, 2], color_temp[, 3],
                   maxColorValue = 255)
  return(color_hex)
}



#' Color mixer
#'
#' @param colors vector of Hex colors
#' @param weight vector of weight of each color in calculating the average
#'
#' @return
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
        linkvalue[k] <- 2
        k <- k + 1
      }
    }
  }
  return(data.frame(Linksource = linksource, Linktarget = linktarget,
                    Linkvalue = linkvalue))

}
