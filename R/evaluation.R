
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
#' tp_data = chicken_generator(1)
#' tp_data_mapper = mapper.kmeans(dat = tp_data[,2:4],
#'                                filter_values = tp_data$Y,
#'                                num_intervals = 10,
#'                                percent_overlap = 70)
#' tp_shannon = mapper_shannon_index(tp_data_mapper, tp_data$Group)
#' tp_shannon
#'
mapper_shannon_index = function(obj_mapper, group_ind) {

  if(class(obj_mapper) != "TDAmapper") {
    stop("Invalid obj_mapper.")
  }

  require(entropy)

  shannon_indices = NULL
  for(i in obj_mapper$points_in_vertex){
    gp = group_ind[i]
    entropy_temp = entropy.empirical(table(gp), unit="log2")
    shannon_indices = rbind(shannon_indices, c(entropy_temp, length(i)))
  }

  colnames(shannon_indices) = c("Index", "Weight")

  res = list(avg_index = sum(shannon_indices[,1] * shannon_indices[,2])/sum(shannon_indices[,2]),
             index_per_node = shannon_indices)
  class(res) <- "shannon_index"
  return(res)
}


#' @export
print.shannon_index = function(x) {
  cat("The weighted average of shannon indices is", signif(x$avg_index, 4), "\n")
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
#' tp_data = chicken_generator(1)
#' tp_data_mapper = mapper.kmeans(dat = tp_data[,2:4],
#'                                filter_values = tp_data$Y,
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

  # Remove the classes not exist in the main graph
  mean_dist[is.na(mean_dist)] = 0
  weight[is.na(mean_dist)] = 0

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
#' filter_names = c("Coordinate","Eccen","LInf","Ref","DTM","Gaussian","PCA")
#' res_filter = data.frame(Filter = filter_names,
#'                         weighted_shannon = c(1.389, 1.421, 1.453, 1.158, 1.345, 1.399, 1.349),
#'                         spread_measure = c(2.767, 4.101, 2.607, 4.001, 4.119, 3.957, 2.034))
#' res = pareto_opt(res_filter)
#' res
#'
#' # Illustration of Pareto frontier
#' library(ggplot2)
#' res = pareto_opt(res_filter, top = nrow(res_filter))
#' res$.level[res$.level>1] = "Others"
#' res$.level[res$.level==1] = "Frontier"
#' class(res) = "data.frame"
#' gp <- ggplot(res, aes(x = weighted_shannon, y = spread_measure,
#'                       color = factor(.level),
#'                       label = Filter)) +
#'   geom_point(size = 3) +
#'   geom_step(aes(group = .level), data = res[res$.level=="Frontier",],
#'             direction = "vh", size = 2) +
#'   geom_label(aes(fill = factor(.level)), colour = "white", fontface = "bold") +
#'   labs(fill='') + labs(color='') + xlab("Shannon index") +
#'   ylab("Spread measure")
#' gp
#' # The red line(s) is the Pareto frontier, and all the points on the top right
#' # side of these line(s) are not potentially optimal.
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

  class(res) = "Pareto_frontier"
  return(res)
}

#' Print Pareto frontier
#'
#' @param res The Pareto_frontier object
#'
#' @return Summary of the Pareto_frontier object
#' @export
#'
#' @seealso \code{\link{pareto_opt}}
print.Pareto_frontier = function(res) {
  cat("Filter functions in the Pareto frontier:\n\t", paste(res$Filter, collapse = ", "), "\n")
  cat("\nEvaluation results of filter functions in the Pareto frontier:\n")
  class(res) = "data.frame"
  print(res)
}



#' Evaluate filter functions
#'
#' \code{filter_evaluate} evaluates filter function with provided filter
#' vectors.
#'
#' @param ... Filter objects. The classes of the objects should be \code{filter}
#'   and include an attribute \code{filter} which is the name of the
#'   corresponding filter function.
#' @inheritParams mapper.kmeans
#' @param arg_mapper A list for additional arguments for
#'   \code{\link{mapper.kmeans}}.
#' @param group_ind A vector of group names each of the samples belongs to.
#'
#' @return A data.frame of Shannon indices and spread measures under given
#'   filter functions. The first column contents names of filter functions, and
#'   the second and third columns are the shannon indices and spread measures,
#'   respectively.
#' @export
#'
#' @examples
#' tp_data = chicken_generator(1)
#' tp_dist = dist(tp_data[,-1])
#' a = filter_eccen(dist = tp_dist, p = 2)
#' b = filter_coordinate(tp_data[,-1], 2)
#' c = filter_gaussian(dist=tp_dist, sigma=1)
#' filter_evaluate(a,b,c,
#'                 dat = tp_data, group_ind = tp_data$Group,
#'                 num_intervals = 10, percent_overlap = 70)
#'
#' # Add additional arguments (NOT RUN)
#' if(FALSE) {
#'   filter_evaluate(a,b,c,
#'                   dat = tp_data, group_ind = tp_data$Group,
#'                   num_intervals = 10, percent_overlap = 70,
#'                   arg_mapper = list(n_class = 1))
#' }
#'
filter_evaluate = function(..., dat, group_ind, num_intervals, percent_overlap,
                           arg_mapper = list()) {

  filter_list = list(...)
  res = NULL

  arg_mapper$dat = dat
  arg_mapper$group_ind = group_ind
  arg_mapper$num_intervals = num_intervals
  arg_mapper$percent_overlap = percent_overlap

  for(ff in filter_list) {
    if(class(ff) != "filter") {
      stop("filter_evaluate can only accept 'filter' object.")
    }

    tp_filter = attributes(ff)$filter
    class(ff) = "numeric"

    arg_mapper$filter_values = ff
    tp_mapper = do.call(mapper.kmeans, arg_mapper)

    tp_shannon = mapper_shannon_index(obj_mapper = tp_mapper, group_ind = group_ind)
    tp_spread = spread_measure(obj_mapper = tp_mapper, group_ind = group_ind)


    res = rbind(res, c(tp_filter, tp_shannon$avg_index, tp_spread))
  }

  colnames(res) = c("Filter", "weighted_shannon", "spread_index")
  res = as.data.frame(res, stringsAsFactors = FALSE)
  res[,2:3] = lapply(res[,2:3], as.numeric)
  return(res)
}
