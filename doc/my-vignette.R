## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE----------------------------------------------------------
#  library(SemiMapper)
#  simu_data <-  chicken_generator(seed = 1)
#  simu_data_mapper <- mapper.kmeans(dat = simu_data[,2:4],
#                                 filter_values = simu_data$Y,
#                                 num_intervals = 10,
#                                 percent_overlap = 70)
#  simple_visNet(simu_data_mapper, filter = simu_data$Y, color_filter = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  simple_visNet(simu_data_mapper, color_filter = FALSE, groups_ind = simu_data$Group)

## ----eval=FALSE----------------------------------------------------------
#  simu_color <- auto_set_colorcode(simu_data$Group, palette = "Set1")
#  check_color_code(simu_color)
#  simple_visNet(simu_data_mapper, color_filter = FALSE,
#                groups_ind = simu_data$Group,color_code = simu_color)

## ----eval = FALSE--------------------------------------------------------
#  tp_data <- chicken_generator(1)
#  ff <- filter_coordinate(tp_data[,-1], 2)
#  tp_data_mapper <- mapper.kmeans(dat = tp_data[,2:4],
#                                 filter_values = ff,
#                                 num_intervals = 10,
#                                 percent_overlap = 70)
#  network_visualization(tp_data_mapper, groups_ind = tp_data$Group, dat = tp_data[,2:4],
#                        folder = "Exp_network")

## ----eval=FALSE----------------------------------------------------------
#  add_analysis_js <- paste0('Node Index:<b>',
#                             1:length(tp_data_mapper$points_in_vertex),
#                             '</b><br>')
#  network_visualization(tp_data_mapper, groups_ind = tp_data$Group, dat = tp_data[,2:4],
#                        folder = "Exp_network", add_analysis_js = add_analysis_js)
#  

## ----eval = FALSE--------------------------------------------------------
#  tp_data <- chicken_generator(1)
#  tp_dist <- dist(tp_data[,-1])
#  ff <- filter_eccen(dist = tp_dist, p = 2)
#  tp_data_mapper <- mapper.kmeans(dat = tp_data[,2:4],
#                                 filter_values = ff,
#                                 num_intervals = 10,
#                                 percent_overlap = 70)
#  network_visualization(tp_data_mapper, groups_ind = tp_data$Group, dat = tp_data[,2:4],
#                        folder = "Exp_network")

## ----eval = FALSE--------------------------------------------------------
#  tp_data <- chicken_generator(1)
#  tp_dist <- dist(tp_data[,-1])
#  ff <- filter_ref(dist=tp_dist, groups_ind=tp_data$Group, ref = "Shank")
#  tp_data_mapper <- mapper.kmeans(dat = tp_data[,2:4],
#                                 filter_values = ff,
#                                 num_intervals = 10,
#                                 percent_overlap = 70)
#  network_visualization(tp_data_mapper, groups_ind = tp_data$Group, dat = tp_data[,2:4],
#                        folder = "Exp_network")

## ----eval=FALSE----------------------------------------------------------
#  filter_eccen = function(dist, p, ...) {
#  
#    dist <- as.matrix(dist)
#    res<- rowMeans(dist^p)^(1/p)
#    res <- as.matrix(res)
#  
#    attr(res, "filter") <- "Eccentricity"
#    class(res) <- "filter"
#    return(res)
#  }

## ----eval=FALSE----------------------------------------------------------
#  tp_data <- chicken_generator(1)
#  tp_data_mapper <- mapper.kmeans(dat = tp_data[,2:4],
#                                  filter_values = tp_data$Y,
#                                  num_intervals = 10,
#                                  percent_overlap = 70)
#  tp_shannon <- mapper_shannon_index(tp_data_mapper, tp_data$Group)
#  tp_shannon

## ----eval=FALSE----------------------------------------------------------
#  tp_data <- chicken_generator(1)
#  tp_data_mapper<- mapper.kmeans(dat = tp_data[,2:4],
#                                 filter_values = tp_data$Y,
#                                 num_intervals = 10,
#                                 percent_overlap = 70)
#  tp_spread <- spread_measure(tp_data_mapper, tp_data$Group)
#  tp_spread

## ------------------------------------------------------------------------
library(SemiMapper)
filter_names <- c("Coordinate","Eccen","LInf","Ref","DTM","Gaussian","PCA")
res_filter <- data.frame(Filter = filter_names,
                        weighted_shannon = c(1.389, 1.221, 1.453, 1.158, 1.345, 1.399, 1.349),
                        spread_measure = c(2.767, 3.101, 2.607, 4.001, 4.119, 3.957, 2.034))
res <- pareto_opt(res_filter)
print(res)

## ------------------------------------------------------------------------
# Illustration of Pareto frontier
library(ggplot2)
res <- pareto_opt(res_filter, top = nrow(res_filter))
res$.level[res$.level>1] <- "Others"
res$.level[res$.level==1] <- "Frontier"
class(res) <- "data.frame"
gp <- ggplot(res, aes(x = weighted_shannon, y = spread_measure,
                      color = factor(.level),
                      label = Filter)) +
  geom_point(size = 3) +
  geom_step(aes(group = .level), data = res[res$.level=="Frontier",],
            direction = "vh", size = 2) +
  geom_label(aes(fill = factor(.level)), colour = "white", fontface = "bold") +
  labs(fill='') + labs(color='') + xlab("Shannon index") +
  ylab("Spread measure")


## ----fig2, fig.height = 2.5, fig.width = 5, fig.align = "center", echo=FALSE----
gp

## ----eval = FALSE--------------------------------------------------------
#  tp_data <- chicken_generator(1)
#  tp_dist <- dist(tp_data[,-1])
#  fe <- filter_eccen(dist = tp_dist, p = 2)
#  fc <- filter_coordinate(tp_data[,-1], 2)
#  fg <- filter_gaussian(dist=tp_dist, sigma=1)
#  res_eval <- filter_evaluate(fe,fc,fg,
#                              dat = tp_data, group_ind = tp_data$Group,
#                              num_intervals = 10, percent_overlap = 70)

## ----eval=FALSE----------------------------------------------------------
#  res <- pareto_opt(res_eval)
#  print(res)

