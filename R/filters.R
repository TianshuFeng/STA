#' PCA filter
#'
#' \code{filter_pca} returns the first k principle components as filter values.
#'
#' The PCA filter function is defined as \eqn{f(x_i) = x_i^T\phi_{(1:k)}}, where
#' \eqn{\phi_{(1:k)}} is the matrix of \eqn{k} eigenvectors associated with the
#' \eqn{k} largest eigenvalues of \eqn{cov(X)}, \eqn{x_i} is some data point and
#' \eqn{X} is the data matrix.
#'
#' @param dat A numeric dataset matrix, rach row represents a data point and
#'   each column represents a predictive variable.
#' @param k A scaler deciding the number of principle components to be returned
#' @param ... Optional arguments to \code{\link[stats]{cov}}.
#'
#' @return A matrix object of filter values.
#'
#' @export
#'
#' @examples
#' tp_data = chicken_generator(1)
#' filter_pca(dat=tp_data[,-1])
#'
filter_pca = function(dat, k = 1, ...) {
  require(rARPACK)
  cov_temp = cov(dat,...)
  eigen_vector = eigs_sym(cov_temp ,k)$vectors

  res = as.matrix(dat) %*% eigen_vector
  res = as.matrix(res)

  attr(res, "filter") = "PCA"
  class(res) = "filter"
  return(res)
}


#' L-infinity filter
#'
#' \code{filter_Linf} computes the L-infinity filter.
#'
#' The L-infinity filter function is defined as \eqn{f(x_i) = max_{j} d(x_i,x_j)} where
#' \eqn{x_i,x_j} are some data points.
#'
#' @param dist The distance matrix.
#' @param ... Optional arguments to \code{\link[base]{max}}.
#'
#' @return A matrix object of filter values.
#' @export
#'
#' @examples
#' tp_data = chicken_generator(1)
#' tp_dist = dist(tp_data[,-1])
#' filter_Linf(dist=tp_dist)
#'
filter_Linf = function(dist, ...) {
  # dist: distance matrix

  res = apply(as.matrix(dist), 1, max, ...)
  res = as.matrix(res)

  attr(res, "filter") = "Linf"
  class(res) = "filter"
  return(res)
}


#' Reference distance filter
#'
#' \code{filter_ref} computes the reference distance filter.
#'
#' Reference distance filter function is defined as \eqn{f(x_i;y,g) =
#' 1-\textrm{median}_{y_j=g} d(x_i,x_j)}, where
#'
#' @param dist The distance matrix.
#' @param groups_ind A vector of group names each of the samples belongs to.
#' @param ref A string object specifying the name of the reference group.
#' @param ... Additional arguments.
#'
#' @return A matrix object of filter values.
#' @export
#'
#' @examples
#' tp_data <- chicken_generator(1)
#' tp_dist <- dist(tp_data[,-1])
#' filter_ref(dist=tp_dist, groups_ind=tp_data$Group, ref = "Shank")
#'
filter_ref = function(dist, groups_ind,
                      ref, ...) {
  #  Find the average distance between samples under center level and other samples
  #  ref_var: the reference variable
  #  center: the reference group

  #center_sample = Metadat[Metadat[,ref_var] == center, id_var]

  # Should keep the order of distance matrix the same as Metadat
  #colnames(dist) = Metadat[,id_var]

  # Distance between samples and center samples
  dist <- as.matrix(dist)
  diag(dist) = NA

  dist_center = dist[,which(groups_ind == ref)]
  rm(dist)
  # find the average distance
  res = apply(dist_center, 1, median, na.rm = TRUE)
  res = as.matrix(res)

  attr(res, "filter") = "Reference"
  class(res) = "filter"
  return(res)

}


#' Gaussian density filter
#'
#' \code{filter_gaussian} computes the Gaussian density filter.
#'
#' The Gaussian density filter is defined as \eqn{f(x_i;\sigma) = C\sum_{j=1}^n
#' exp(-\frac{\|x_i-x_j\|^2}{2\sigma^2})}, where \eqn{C} is the normalizing
#' constant and \eqn{\sigma} is the scalar controling the sensitivity of the
#' gaussian kernel.
#'
#' @param dist The distance matrix.
#' @param sigma A scalar controling the sensitivity of the gaussian kernel.
#' @param ... Further arguments.
#'
#' @return A matrix object of filter values.
#' @export
#'
#' @examples
#' tp_data = chicken_generator(1)
#' tp_dist = dist(tp_data[,-1])
#' filter_gaussian(dist=tp_dist, sigma=1)
#'
filter_gaussian = function(dist, sigma = 1,...) {
  # sigma: control the sensitivity of gaussian kernel
  dist = as.matrix(dist)
  dist = exp(-dist^2/(2*sigma^2))
  ff = rowSums(dist)
  ff = ff/sum(ff)

  ff = as.matrix(ff)
  attr(ff, "filter") = "Gaussian"
  class(ff) = "filter"
  return(ff)
}


#' Coordinate projection filter
#'
#' \code{filter_coordinate} computes the Coordinate projection filter.
#'
#' Coordinate projection filter is defined as \eqn{f(x_i;k) = x_{ik}}.
#'
#'
#' @inheritParams filter_pca
#' @param k A scalar or vector deciding to which coordinate(s) the data matrix is
#'   projected.
#' @param ... Further arguments.
#'
#' @return A matrix object of filter values.
#' @export
#'
#' @examples
#' tp_data = chicken_generator(1)
#' filter_coordinate(tp_data[,-1],1)
#'
filter_coordinate = function(dat, k, ...) {
  # k: the coordinate to be projected

  res = dat[,k,drop=FALSE]
  res = as.matrix(res)

  attr(res, "filter") = "Coordinate"
  class(res) = "filter"
  return(res)
}


#' Distance to measure filter
#'
#' \code{filter_dtm} computes the Distance to measure filter.
#'
#' Distance to measure filter function is defined as \eqn{f(x_i;k) =
#' \frac{1}{k}\sum_{j=2}^{k+1} d^p(x_i,x_{(j)})^{1/p}}.
#'
#'
#' @param dist The distance matrix.
#' @param k A numeric scalar deciding the number of neighbors to be considered
#' @param p An optional numeric scalar denoting the exponent.
#' @param ... Further arguments.
#'
#' @return A matrix object of filter values.
#' @export
#'
#' @examples
#' tp_data = chicken_generator(1)
#' tp_dist = dist(tp_data[,-1])
#' filter_dtm(dist = tp_dist, k=3)
filter_dtm = function(dist, k, p=2, ...) {
  # k: the number of neighbors to be considered
  dist = as.matrix(dist)
  dist_sort = t(apply(dist, 1, sort))
  dist_sort = dist_sort[,2:k+1]^p
  dist_sort = rowMeans(dist_sort)

  res = (dist_sort)^(1/p)
  res = as.matrix(res)

  class(res) = "filter"
  attr(res, "filter") = "DTM"
  return(res)
}


#' Eccentricity filter
#'
#' \code{filter_eccen} computes the Eccentricity filter.
#'
#' Eccentricity filter function is defined as \eqn{f(x_i;p) =
#' (\frac{1}{n}\sum_{j=1}^n d^p(x_i,x_j))^{1/p}}.
#'
#' @param dist The distance matrix.
#' @param p A numeric scalar denoting the exponent.
#' @param ... Further arguments.
#'
#' @return A matrix object of filter values.
#' @export
#'
#' @examples
#' tp_data = chicken_generator(1)
#' tp_dist = dist(tp_data[,-1])
#' filter_eccen(dist = tp_dist, p = 2)
#'
filter_eccen = function(dist, p, ...) {
  # p: exponent
  dist = as.matrix(dist)
  res = rowMeans(dist^p)^(1/p)
  res = as.matrix(res)

  attr(res, "filter") = "Eccentricity"
  class(res) = "filter"
  return(res)
}


#' Which filter function produces the filter?
#'
#' \code{filter_name} returns the name of the filter function that produces the
#' provided filter.
#'
#' @param filter_vec The filter vector.
#'
#' @return The name of the filter.
#' @export
#'
#' @examples
#' tp_data = chicken_generator(1)
#' tp_dist = dist(tp_data[,-1])
#' res = filter_eccen(dist = tp_dist, p = 2)
#' filter_name(res)
#' ## "Eccentricity"
filter_name = function(filter_vec) {
  if(!("filter" %in% names(attributes(filter_vec)))) {
    stop("Filter vector not from standard filter functions.")
  }

  return(attr(filter_vec, "filter"))
}

#' Print filter function name
#'
#' @param x a filter object
#'
#' @return the name of the filter function producing the filter
#' @export
#'
#' @seealso \code{\link{filter_pca}}
print.filter = function(x) {
  cat("This is the filter vecter returned by the ",
      attributes(x)$filter, " filter function.\n", sep = "")
}
