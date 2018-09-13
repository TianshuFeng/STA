
#' Remove NULL nodes in the graph
#'
#' @param obj_mapper The Mapper object
#'
#' @return A Mapper object
#'
null_remover = function(obj_mapper) {

  null_idx = which(sapply(obj_mapper$points_in_vertex,is.null),
                   arr.ind=TRUE)
  if(length(null_idx) > 0){
    obj_mapper$adjacency = obj_mapper$adjacency[-null_idx, -null_idx]

    obj_mapper$points_in_vertex =
      obj_mapper$points_in_vertex[-null_idx]

    obj_mapper$num_vertices = length(obj_mapper$points_in_vertex)

    obj_mapper$level_of_vertex = obj_mapper$level_of_vertex[-null_idx]
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
color_map_Spectral = function(x){
  require(RColorBrewer)
  color_temp = colorRamp(brewer.pal(11,"Spectral"))(x)
  color_hex = rgb(color_temp[,1], color_temp[,2], color_temp[,3], maxColorValue = 255)
  return(color_hex)
}


savepdf <- function(file, folder = "", width=10, height=10)
{
  dir.create(file.path(folder, "figures"), showWarnings = FALSE)
  fname <- paste(folder, "/figures/",file,".png",sep="")
  png(fname)
  par(mgp=c(0,0,0), tcl=-0.4, mar=c(0,0,0,0), oma = c(0,0,0,0), bg=NA)

}

