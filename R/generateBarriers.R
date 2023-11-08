#' Generate barriers from sf file and permeability given by user
#'
#' @name generateBarriers
#' @description This R function generates a barrier dataset from any input
#' barrier shapefile and user-defined list of corresponding permeability values
#' (or one value to be applied to the entire list)
#' @param barriers Barriers provided by user in sf format.
#' @param p_list Either a list of permeability values (that should correspond
#' to the length of barriers) or a single value (to be applied to all barriers)
#' @return a matrix of six columns:
#' @useDynLib abmFences
#' @export
#'

generateBarriers <- function(barriers, p_list) {
  # inputs are a shapefile of barrier data and a list of accompanying
  # permeabilities for each one. If all are the same for one barrier type,
  # then we only have to run it once for the whole dataset
  LP <- length(p_list)
  if (class(barriers)[1] == 'list') NB <- length(barriers)
  else NB <- nrow(barriers)
  if (LP == 1) {
    message("NOTE: setting all barrier permeabilities to ", p_list, '.')
    barrier.df <- .segmentBarrier(barriers, p_list)
  } else if (LP != NB) {
    message("WARNING! length of permeability vector must be 1 or #barriers")
    return()
  } else {
    message("Changing permeability based on barrier.")
    for (i in 1:NB) {
      b <- .segmentBarrier(barriers[i,], p_list[i], i)
      if (i == 1) barrier.df <- b
      else barrier.df <- rbind(barrier.df, b)
    }
  }
  return(barrier.df)
}

# helper function not to export
.segmentBarrier <- function(fence, p_cross, id=0) {
  data <- fences %>% st_segments() %>% st_coordinates()
  NR <- nrow(data)
  my.fence <- matrix(data=0, nrow=NR/2, ncol=6)
  my.fence[,1:2] <- data[seq(1, NR, by=2), 1:2]
  my.fence[,3:4] <- data[seq(2, NR, by=2), 1:2]
  my.fence[,5] <- p_cross
  my.fence[,6] <- id
  colnames(my.fence) = c('x', 'y', 'xend', 'yend', 'perm', 'id')
  return(my.fence)
}
