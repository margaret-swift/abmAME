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

generateBarriers <- function(barriers, rast, p_list, ...) {
  # inputs are a shapefile of barrier data and a list of accompanying
  # permeabilities for each one. If all are the same for one barrier type,
  # then we only have to run it once for the whole dataset
  LP <- length(p_list)
  if (class(barriers)[1] == 'list') NB <- length(barriers)
  else NB <- nrow(barriers)
  if (LP != NB) {
    message("WARNING! length of permeability vector must be 1 or #barriers")
    return()
  } else {
    if (LP == 1) {
      message("NOTE: setting all barrier permeabilities to ", p_list, '.')
      p <- p_list
    }
    for (i in 1:NB) {
      if (LP > 1) p = p_list[i]
      b <- barriers[[i]]
      segs <- .segmentBarriers(b, i, ...)
      mat  <- .createMatrix(segs, p)
      lookup<-.createLookup(segs, rast)

      if (i == 1) {
        barrier.df <- mat
        lookup.df <- lookup
      } else {
        barrier.df <- rbind(barrier.df, mat)
        lookup.df <- rbind(lookup.df, lookup)
      }
    }
  }
  message(dim(barrier.df))
  message('barriers')
  print(head(barrier.df))
  message('lookup')
  print(head(lookup))
  # only keep barriers that are within the environmental raster.
  message('Removing barriers outside environment...')
  inx.keep <- which(barrier.df[,6] %in% lookup[,3])
  barrier.df <- barrier.df[inx.keep,]
  return(list(barrier.df, lookup.df))
}

# helper functions not to export
.segmentBarriers <- function(b, id=0, ...) {
  data <- b %>% st_segments(...)
  NR = nrow(data)
  dig <- floor(log10(NR)) + 1
  data$INX <- id + (10 ^ -dig) * (1:NR)
  return(data)
}

.createMatrix <- function(data, p) {
  coords <- data %>% st_coordinates()
  NR <- nrow(coords)
  # fill in barrier matrix
  b.mat <- matrix(data=0, nrow=NR/2, ncol=6)
  b.mat[,1:2] <- coords[seq(1, NR, by=2), 1:2]
  b.mat[,3:4] <- coords[seq(2, NR, by=2), 1:2]
  b.mat[,5] <- p
  b.mat[,6] <- data$INX
  colnames(b.mat) = c('x', 'y', 'xend', 'yend', 'perm', 'id')
  b.mat
}

.createLookup <- function(data, rast) {
  tab <- terra::extract(rast, data)
  tab$INX <- data$INX[tab$ID]
  inx.na <- is.na(tab[,2])
  tab <- tab[!inx.na,]
  tab <- tab[order(tab[,2]),]
  tab <- as.matrix(tab)
  return(tab)
}
