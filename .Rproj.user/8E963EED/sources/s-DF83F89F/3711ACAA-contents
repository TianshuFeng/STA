lsmi_from_lsfi <- function( lsfi, num_intervals ) {
  # inputs:
  # lsfi = an integer in the range 1:prod(v)
  # num_intervals = c(i1,i1,...) a vector of numbers of intervals
  # output:
  # f+1 = a vector of multiindices with length filter_output_dim
  j <- c(1,num_intervals) # put 1 in front to make indexing easier in the product prod(j[1:k])
  f <- c()
  for (k in 1:length(num_intervals)) {
    # use lsfi-1 to shift from 1-based indexing to 0-based indexing
    f[k] <- floor( (lsfi-1) / prod(j[1:k])) %% num_intervals[k]
  }
  #print(f+1)
  # lsmi = f+1 = level set multi index
  return(f+1) # shift from 0-based indexing back to 1-based indexing
}


lsfi_from_lsmi <- function( lsmi, num_intervals ) {
  lsfi <- lsmi[1]
  if (length(num_intervals) > 1) {
    for (i in 2:length(num_intervals)) {
      lsfi <- lsfi + prod(num_intervals[1:(i-1)]) * (lsmi[i]-1)
    }
  }
  return(lsfi)
}


#' Mapper function with multiple cluster methods
#'
#' This function is adopted from \code{mapper} function of \code{TDAmapper} with
#' different clustering methods (mainly k-means).
#'
#' This function is adopted from \code{mapper} function of \code{TDAmapper} by
#' replacing its cluster method with the cluster function
#' \code{\link[NbClust]{NbClust}} from R package \code{NbClust}.
#'
#' The advantage of \code{NbClust} is that it provides users with 8 different
#' cluster methods, 6 different distance measures and 30 indices for determining
#' the number of clusters. This allows users to select the best clustering
#' scheme from the different results obtained by varying all combinations of
#' number of clusters, distance measures, and clustering methods. Details of the
#' distance measures, clustering methods and cluster indices can be found in
#' \code{\link[NbClust]{NbClust}}.
#'
#' @inheritParams TDAmapper::mapper
#' @param dat Matrix or dataset where rows are data points and columns are
#'   predictive variables.
#' @param dist_method The distance measure to be used to compute the
#'   dissimilarity matrix. By default, distance="euclidean". It must be one of
#'   This must be one of: "euclidean", "maximum", "manhattan", "canberra",
#'   "binary", "minkowski" or "NULL". Details can be found in
#'   \code{\link[NbClust]{NbClust}}.
#' @param cluster_method The cluster analysis method to be used. This should be
#'   one of: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty",
#'   "median", "centroid", "kmeans".Details can be found in
#'   \code{\link[NbClust]{NbClust}}.
#' @param cluster_index The index to be calculated. This should be one of :
#'   "kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew",
#'   "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2",
#'   "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain",
#'   "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw",
#'   "all" (all indices except GAP, Gamma, Gplus and Tau), "alllong" (all
#'   indices with Gap, Gamma, Gplus and Tau included). Details can be found in
#'   \code{\link[NbClust]{NbClust}}.
#' @param n_class number of clusters. By default, n_class=0. If n_class>0, this
#'   function will instead call \code{\link[stats]{kmeans}} and pass
#'   \code{n_class} to argument \code{centers} of \code{\link[stats]{kmeans}}.
#' @param ... Further arguments for either \code{\link[NbClust]{NbClust}} or
#'   \code{\link[stats]{kmeans}}.
#'
#' @return An object of class \code{TDAmapper} which is a list of items named
#'   \code{adjacency} (adjacency matrix for the edges), \code{num_vertices}
#'   (integer number of vertices), \code{level_of_vertex} (vector with
#'   \code{level_of_vertex[i]} = index of the level set for vertex i),
#'   \code{points_in_vertex} (list with \code{points_in_vertex[[i]]} = vector of
#'   indices of points in vertex i), \code{points_in_level} (list with
#'   \code{points_in_level[[i]]} = vector of indices of points in level set i,
#'   and \code{vertices_in_level} (list with \code{vertices_in_level[[i]]} =
#'   vector of indices of vertices in level set i.
#'
#' @export
#'
#' @references Malika Charrad, Nadia Ghazzali, Veronique Boiteau, Azam Niknafs
#'   (2014). NbClust: An R Package for Determining the Relevant Number of
#'   Clusters in a Data Set. Journal of Statistical Software, 61(6), 1-36. URL
#'   http://www.jstatsoft.org/v61/i06/.
#'
#' @examples
mapper.kmeans <- function(dat, dist_method = "euclidean", filter_values, num_intervals, percent_overlap,
                          cluster_method = "kmeans", cluster_index = "all",
                          n_class = 0, ...) {

  # ...: further argument for nbclust

  ##### begin documentation ############
  # inputs
  # f : X \subset R^n \to R^k, a filter function on a data set with numpoints observations
  # filter_values = data.frame(y_1, y_2,..., y_k), where each y_i is a vector of length num_points
  # num_intervals = c(i_1, i_2,..., i_k), a vector of number of intervals for each variable y_i
  # percent_overlap = c(p_1, p_2,..., p_k), a vector of percent overlap for adjacent intervals within each variable y_i
  ##### end documentation ###############


  #     #filter_output_dim <- length(filter_values)
  #     if (length(num_intervals) == 1) {
  #         num_points <- length(filter_values)
  #         filter_output_dim <- 1
  #         num_levelsets <- num_intervals
  #
  #         # define some vectors of length k = number of columns = number of variables
  #         filter_min <- min(filter_values)
  #         filter_max <- max(filter_values)
  #         interval_width <- (filter_max - filter_min) / num_intervals
  #
  #         } else {
  # #    filter_values <- as.matrix(filter_values)
  #         num_points <- dim(filter_values)[1] # number of rows = number of observations
  #         filter_output_dim <- dim(filter_values)[2] # number of columns = number of variables = length(num_intervals)
  #         num_levelsets <- prod(num_intervals)
  #
  #         # define some vectors of length k = number of columns = number of variables
  #         filter_min <- as.vector(sapply(filter_values,min))
  #         filter_max <- as.vector(sapply(filter_values,max))
  #         interval_width <- (filter_max - filter_min) / num_intervals
  #
  #    }

  # class(filter_values[,1]) = numeric, which has dim(filter_values[,1]) = NULL,
  # so we coerce filter_values to a data.frame so that its dim is not NULL
  require(NbClust)
  if(n_class == 0){
    cat("No number of cluters specified, use NbClust instead. \n")}

  filter_values <- data.frame(filter_values)
  num_points <- dim(filter_values)[1] # number of rows = number of observations
  filter_output_dim <- dim(filter_values)[2] # number of columns = number of variables = length(num_intervals)
  num_levelsets <- prod(num_intervals)

  # define some vectors of length k = number of columns = number of variables
  filter_min <- as.vector(sapply(filter_values,min))
  filter_max <- as.vector(sapply(filter_values,max))
  interval_width <- (filter_max - filter_min) / num_intervals


  # initialize variables
  vertex_index <- 0
  level_of_vertex <- c()
  points_in_vertex <- list()
  points_in_level_set <- vector( "list", num_levelsets )
  vertices_in_level_set <- vector( "list", num_levelsets )
  # for future development
  # cutree_in_level_set <- vector( "list", num_levelsets )


  #### begin plot the filter function ##############
  #     # Reality check
  #     # Plot the filter values
  #     plot(filter_values[,1], filter_values[,2], type="n")
  #     # cex = font size as a proportion of default
  #     text(filter_values[,1], filter_values[,2], labels=1:num_points, cex=0.5)
  #     # midpoint of overlapping intervals
  #     abline(v = filter_min[1]+interval_width[1]*(0:num_intervals[1]),
  #            h = filter_min[2]+interval_width[2]*(0:num_intervals[2]), col="red")
  #     # left and right interval boundaries
  #     abline(v = filter_min[1]+interval_width[1]*(0:num_intervals[1])
  #            -0.5*interval_width[1]*percent_overlap[1]/100, col = "blue", lty = 3)
  #     abline(v = filter_min[1]+interval_width[1]*(0:num_intervals[1])
  #            +0.5*interval_width[1]*percent_overlap[1]/100,
  #            col = "blue", lty = 3)
  #     # bottom and top interval boundaries
  #     abline(h = filter_min[2]+interval_width[2]*(0:num_intervals[2])
  #            -0.5*interval_width[2]*percent_overlap[2]/100, col = "blue", lty = 3)
  #     abline(h = filter_min[2]+interval_width[2]*(0:num_intervals[2])
  #            +0.5*interval_width[1]*percent_overlap[2]/100,
  #            col = "blue", lty = 3)
  #### end plot the filter function ##########



  # begin loop through all level sets
  for (lsfi in 1:num_levelsets) {

    ################################
    # begin covering

    # level set flat index (lsfi), which is a number, has a corresponding
    # level set multi index (lsmi), which is a vector
    lsmi <- lsmi_from_lsfi( lsfi, num_intervals )

    lsfmin <- filter_min + (lsmi - 1) * interval_width - 0.5 * interval_width * percent_overlap/100
    lsfmax <- lsfmin + interval_width + interval_width * percent_overlap/100

    # begin loop through all the points and assign them to level sets
    for (point_index in 1:num_points) {
      # compare two logical vectors and get a logical vector,
      # then check if all entries are true
      if ( all( lsfmin <= filter_values[point_index,] &
                filter_values[point_index,] <= lsfmax ) ) {
        points_in_level_set[[lsfi]] <- c( points_in_level_set[[lsfi]],
                                          point_index )
      }
    }
    # end loop through all the points and assign them to level sets

    # end covering
    ######################################

    ######################################
    # begin clustering

    points_in_this_level <- points_in_level_set[[lsfi]]
    num_points_in_this_level <- length(points_in_level_set[[lsfi]])

    if (num_points_in_this_level == 0) {
      num_vertices_in_this_level <- 0
    }

    if (num_points_in_this_level <4 | num_points_in_this_level < n_class) {
      #warning('Level set has only one point')
      num_vertices_in_this_level <- 1
      level_internal_indices <- c(1)
      level_external_indices <- points_in_level_set[[lsfi]]
    }

    if (num_points_in_this_level > 3 & num_points_in_this_level > n_class) {
      # heirarchical clustering
      # level_dist_object <- as.dist(
      #   as.matrix(dist_object)[points_in_this_level,points_in_this_level])
      # level_max_dist <- max(level_dist_object)
      # level_hclust   <- hclust( level_dist_object, method="single" )
      # level_heights  <- level_hclust$height
      #
      # # cut the cluster tree
      # # internal indices refers to 1:num_points_in_this_level
      # # external indices refers to the row number of the original data point
      # level_cutoff   <- cluster_cutoff_at_first_empty_bin(level_heights, level_max_dist, num_bins_when_clustering)
      # level_external_indices <- points_in_this_level[level_hclust$order]
      # level_internal_indices <- as.vector(cutree(list(
      #   merge = level_hclust$merge,
      #   height = level_hclust$height,
      #   labels = level_external_indices),
      #   h=level_cutoff))
      # num_vertices_in_this_level <- max(level_internal_indices)

      if(n_class == 0){

        pdf(file = NULL)
        options(warn=-1)
        log = capture.output({
          level_kmeans = try(NbClust(data = dat[points_in_this_level,],
                                     distance = dist_method,
                                     method = cluster_method, index = cluster_index, ...),
                             silent = TRUE)
        })
        dev.off()
        options(warn=0)

        if("try-error" %in% class(level_kmeans)) {
          level_internal_indices <- c(1)
        } else {
          level_internal_indices = as.numeric(level_kmeans$Best.partition)
        }
      } else {
        level_kmeans = kmeans(x = dat[points_in_this_level,],
                              centers = n_class, ...)
        level_internal_indices = level_kmeans$cluster
      }

      level_external_indices = points_in_this_level
      num_vertices_in_this_level <- max(level_internal_indices)
    }




    # end clustering
    ######################################

    ######################################
    # begin vertex construction

    # check admissibility condition
    if (num_vertices_in_this_level > 0) {

      vertices_in_level_set[[lsfi]] <- vertex_index + (1:num_vertices_in_this_level)

      for (j in 1:num_vertices_in_this_level) {

        vertex_index <- vertex_index + 1
        level_of_vertex[vertex_index] <- lsfi
        points_in_vertex[[vertex_index]] <- level_external_indices[level_internal_indices == j]

      }
    }

    # end vertex construction
    ######################################

  } # end loop through all level sets


  ########################################
  #  begin simplicial complex

  # create empty adjacency matrix
  adja <- mat.or.vec(vertex_index, vertex_index)

  # loop through all level sets
  for (lsfi in 1:num_levelsets) {

    # get the level set multi-index from the level set flat index
    lsmi <- lsmi_from_lsfi(lsfi,num_intervals)

    # Find adjacent level sets +1 of each entry in lsmi
    # (within bounds of num_intervals, of course).
    # Need the inverse function lsfi_from_lsmi to do this easily.
    for (k in 1:filter_output_dim) {

      # check admissibility condition is met
      if (lsmi[k] < num_intervals[k]) {
        lsmi_adjacent <- lsmi + diag(filter_output_dim)[,k]
        lsfi_adjacent <- lsfi_from_lsmi(lsmi_adjacent, num_intervals)
      } else { next }

      # check admissibility condition is met
      if (length(vertices_in_level_set[[lsfi]]) < 1 |
          length(vertices_in_level_set[[lsfi_adjacent]]) < 1) { next }

      # construct adjacency matrix
      for (v1 in vertices_in_level_set[[lsfi]]) {
        for (v2 in vertices_in_level_set[[lsfi_adjacent]]) {
          adja[v1,v2] <- (length(intersect(
            points_in_vertex[[v1]],
            points_in_vertex[[v2]])) > 0)
          adja[v2,v1] <- adja[v1,v2]
        }
      }

    }


  }

  #  end simplicial complex
  #######################################

  mapperoutput <- list(adjacency = adja,
                       num_vertices = vertex_index,
                       level_of_vertex = level_of_vertex,
                       points_in_vertex = points_in_vertex,
                       points_in_level_set = points_in_level_set,
                       vertices_in_level_set = vertices_in_level_set
  )

  class(mapperoutput) <- "TDAmapper"

  return(mapperoutput)


} # end mapper function
