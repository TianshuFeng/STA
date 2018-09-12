filter_pca = function(dat, k = 1, ...) {
  require(rARPACK)
  # dat: the numeric dataset matrix, rach row represents a sample vector
  # ...: pass to cov function
  cov_temp = cov(dat,...)
  cat("Cov complete")
  eigen_vector = eigs_sym(cov_temp ,k)$vectors

  return(as.matrix(dat) %*% eigen_vector)
}


filter_Linf = function(dist, ...) {
  # dist: distance matrix
  return(apply(dist, 1, max))
}


filter_ref = function(Metadat, dist,
                      id_var = "bcr_patient_barcode",
                      ref_var = "acronym",
                      center = "OV", ...) {
  #  Find the average distance between samples under center level and other samples
  #  ref_var: the reference variable
  #  center: the reference group

  center_sample = Metadat[Metadat[,ref_var] == center, id_var]

  # Should keep the order of distance matrix the same as Metadat
  colnames(dist) = Metadat[,id_var]

  diag(dist) = NA

  # Distance between samples and center samples
  dist_center = dist[,center_sample]
  rm(dist)
  # find the average distance
  return(apply(dist_center, 1, median, na.rm = TRUE))

}


filter_gaussian = function(dist, sigma = 1,...) {
  # sigma: control the sensitivity of gaussian kernel
  dist = exp(-dist^2/(2*sigma^2))
  ff = rowSums(dist)
  ff = ff/sum(ff)
  return(ff)
}


filter_coordinate = function(dat, k, ...) {
  # k: the coordinate to be projected

  return(dat[,k])
}


filter_dtm = function(dist, k, p=2, ...) {
  # k: the number of neighbors to be considered

  dist_sort = t(apply(dist, 1, sort))
  dist_sort = dist_sort[,2:k+1]^p
  dist_sort = rowMeans(dist_sort)
  return((dist_sort)^(1/p))
}


filter_eccen = function(dist, p, ...) {
  # p: exponent

  return(rowMeans(dist^p)^(1/p))
}
