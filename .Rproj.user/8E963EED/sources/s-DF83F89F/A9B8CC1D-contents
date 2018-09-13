
#' Shannon index for Mapper
#'
#' \code{mapper_shannon_index} computes Shannon indices for nodes of a given
#' graph from Mapper.
#'
#' how it is computed
#'
#' @param obj_mapper An object of class \code{TDAmapper}.
#' @param group_ind A vector of group names each of the samples belongs to.
#'
#' @return A list object with two objects. The first object is the weighted
#'   average Shannon index for the whole graph. The second object is the matrix
#'   of Shannon indices for each of the nodes in the graph, where the first
#'   column is the vector of indices and the second column is the vector of node
#'   sizes.
#'
#' @export
#'
#' @examples
#' library(TDAmapper)
#' tp_data = chicken_generator(1)
#' tp_data_mapper = mapper.kmeans(dat = tp_data[,2:4],
#'                                filter_values = tp_data$Z,
#'                                num_intervals = 10,
#'                                percent_overlap = 70)
#' tp_shannon = mapper_shannon_index(tp_data_mapper, tp_data$Group)
#' tp_shannon
mapper_shannon_index = function(obj_mapper, group_ind) {

  if(class(obj_mapper) != "TDAmapper") {
    stop("Invalid obj_mapper.")
  }

  require(entropy)

  obj_mapper = null_remover(obj_mapper)

  shannon_indices = NULL
  for(i in obj_mapper$points_in_vertex){
    gp = group_ind[i]
    entropy_temp = entropy.empirical(table(gp), unit="log2")
    shannon_indices = rbind(shannon_indices, c(entropy_temp, length(i)))
  }

  colnames(shannon_indices) = c("Index", "Weight")
  return(list(avg_index = sum(shannon_indices[,1] * shannon_indices[,2])/sum(shannon_indices[,2]),
              index_per_node = shannon_indices))
}



#' Spread measure function
#'
#' \code{spread_measure} measures how data points from the same groups
#' concentrate in the network.
#'
#' How it is computed
#'
#' @inheritParams mapper_shannon_index
#'
#' @return A scalar which is the spread measure of the graph.
#' @export
#'
#' @examples
#' library(TDAmapper)
#' tp_data = chicken_generator(1)
#' tp_data_mapper = mapper.kmeans(dat = tp_data[,2:4],
#'                                filter_values = tp_data$Z,
#'                                num_intervals = 10,
#'                                percent_overlap = 70)
#' tp_spread = spread_measure(tp_data_mapper, tp_data$Group)
#' tp_spread
#'
spread_measure = function(obj_mapper, group_ind) {

  if(class(obj_mapper) != "TDAmapper") {
    stop("Invalid obj_mapper.")
  }

  require(igraph)

  obj_mapper = null_remover(obj_mapper)

  ad_igraph <- graph_from_adjacency_matrix(as.matrix(obj_mapper$adjacency),
                                           mode="undirected")
  comp_mapper = components(graph = ad_igraph)


  groups_vertices = data.frame(matrix(0, nrow = length(unique(group_ind))))
  colnames(groups_vertices) = "group_ind"
  groups_vertices[,1] = unique(group_ind)

  for(i in 1:length(obj_mapper$points_in_vertex) ) {
    temp_ind = obj_mapper$points_in_vertex[[i]]
    temp_tb = as.data.frame(table(group_ind[temp_ind]))
    colnames(temp_tb) = c("group_ind",i)
    groups_vertices = merge(groups_vertices, temp_tb, by = "group_ind", all = TRUE)
  }

  groups_vertices[is.na(groups_vertices)] = 0

  rownames(groups_vertices) = groups_vertices[,1]
  groups_vertices = groups_vertices[,-1]

  mean_dist = c()
  main_network = comp_mapper$membership == which.max(comp_mapper$csize)
  for(i in 1:nrow(groups_vertices)) {
    dist_shortest = distances(graph = ad_igraph,
                              v = which(groups_vertices[i,]>0 & main_network),
                              to = which(groups_vertices[i,]>0 & main_network))

    v_temp = as.matrix(groups_vertices[i,groups_vertices[i,]>0 & main_network,drop = F])
    dist_shortest = sum(diag(v_temp[1,]) %*% dist_shortest %*% diag(v_temp[1,]))/
      sum(t(v_temp)%*%(v_temp))
    mean_dist = c(mean_dist, dist_shortest)
  }

  weight = table(group_ind)[rownames(groups_vertices)]
  return(sum(weight * mean_dist)/sum(weight))
}


#' Find the pareto frontier
#'
#' \code{pareto_opt} finds the pareto frontier of filter functions based on
#' Shannon indices and spread measures with function \code{\link[rPref]{psel}}
#' from R package `rPref`.
#'
#' @param res_filter The data.frame of Shannon indices and spread measures under
#'   different filter functions. The first column must be names of filter
#'   functions, and the second and third columns should be the shannon indices
#'   and spread measures.
#' @param ... Additional arguments for \code{\link[rPref]{psel}}.
#'
#' @return A data.frame of the names of filter functions in the pareto frontier,
#'   as well as their Shannon indices and spread measures.
#'
#' @export
#'
#' @examples
#' filter_names = c("ref","Linf","gaussian","coordinate","dtm","eccen","pca")
#' res_filter = data.frame(Filter = filter_names,
#'                         weighted_shannon = c(1.708074, 1.753691, 1.788572, 1.783137, 1.798342, 1.729896, 1.645932),
#'                         spread_index = c(3.675639, 4.459636, 3.194998, 3.355519, 3.583395, 3.868635, 5.459008))
#' pareto_opt(res_filter)
#'
pareto_opt = function(res_filter, ...) {
  require(rPref)
  # First column must be the names of filters

  if(ncol(res_filter) != 3) {
    stop("There should be 3 columns: name of filters, shannon indices and spread measures.")
  }

  classes = sapply(res_filter, class)

  if(sum(classes == "numeric")<2) {
    stop("Both shannon indices and spread measures should be numeric.")
  }


  c_name = colnames(res_filter)
  pref = low_(as.expression(c_name[2])) * low_(as.expression(c_name[3]))
  res = psel(res_filter, pref, ...)

  return(res)
}
