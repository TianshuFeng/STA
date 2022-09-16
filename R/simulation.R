#' Chicken foot simulation data generator
#'
#' \code{chicken_generator} generates chicken foot simulation data.
#'
#' \code{chicken_generator} generates chicken foot simulation data with 450 data
#' points in total. It consists of three toes and one shank. Each of the toes
#' contains 105 data points and the shank contains 135 data points.
#'
#' @param seed A scaler determining the seed to be used.
#'
#' @return A dataframe of the simulated data. The first column denotes the
#'   groups of the points and the last three columns are the coordinates of data
#'   points.
#' @export
#'
#' @examples
#' chicken_generator(1)
#'
chicken_generator <- function(seed = 1) {
  library(MASS)
  set.seed(seed)

  x_1 <- 0:20
  y_1 <- 2 * x_1
  d1 <- cbind(1, x_1, y_1)

  x_2 <- -x_1
  y_2 <- y_1
  d2 <- cbind(2, x_2, y_2)

  x_3 <- 0
  y_3 <- 0:(-20)
  d3 <- cbind(4, x_3, y_3)

  x_4 <- 0
  y_4 <- 0:20 * 2
  d4 <- cbind(3, x_4, y_4)

  d <- rbind(d1, d2, d3, d4)

  d <- cbind(d, 0)

  dd <- NULL
  for (i in 1:nrow(d)) {
    tmp <- mvrnorm(n = 5,
                   mu = as.vector(d[i, 2:4]),
                   Sigma = 0.7 *
                     diag(3))
    tmp <- cbind(d[i, 1], tmp)
    dd <- rbind(dd, (tmp))
  }

  dd <-
    rbind(dd, cbind(4, mvrnorm(
      n = 30,
      mu = c(0,-0.5, 0),
      Sigma = diag(c(16, 8, 1))
    )))

  dd <- as.data.frame(dd)
  colnames(dd) <- c("Group", "X", "Y", "Z")

  dd$Y <- dd$Y - min(dd$Y)
  dd$X <- dd$X - min(dd$X)
  dd$Z <- dd$Z - min(dd$Z)

  dd$X <- dd$X / max(dd$X + 0.01)
  dd$Y <- dd$Y / max(dd$Y + 0.01)
  dd$Z <- dd$Z / max(dd$Z + 0.01) / 5

  dd$Group[dd$Group == 1] <- "Toe 1"
  dd$Group[dd$Group == 2] <- "Toe 2"
  dd$Group[dd$Group == 3] <- "Toe 3"
  dd$Group[dd$Group == 4] <- "Shank"

  return(dd)
}
