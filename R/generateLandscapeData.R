#' Generate barriers and landscape data matrices from sf file and permeability given by user
#'
#' @name generateLandscapeData
#' @description This R function generates a barrier dataset from any input
#' barrier shapefile and user-defined list of corresponding permeability values
#' (or one value to be applied to the entire list)
#' @param barriers Barriers provided by user in sf format.
#' @param rasters Landscape raster data
#' @param p_list Either a list of permeability values (that should correspond
#' to the length of barriers) or a single value (to be applied to all barriers)
#' @param crs a CRS for reprojection
#' @return a list with: landscape data matrices, landscape data extent and resolution, and barrier data frame
#' @useDynLib abmAME
#' @export
generateLandscapeData <- function(rasters, barriers, perms=1, crs=32735){
  require(terra)
  require(sf)
  # local functions
  .projectMe <- function(obj) {
    mycrs = crs(obj, describe=T)$code
    if (mycrs != crs) {
      if (class(obj)[1] == "SpatRaster") obj = project(obj, paste0("EPSG:", crs))
      if (class(obj)[1] == "sf") obj = st_transform(obj, crs)
    }
    return(obj)
  }

  message("== RASTERS ==")
  message("reprojecting rasters if necessary...")
  rast.reproj = lapply(rasters, .projectMe)

  message("acquiring extent and resolution of rasters...")
  rast <- rast.reproj[[1]]
  ext <- as.list(ext(rast))
  envExt <-c(unlist(ext, use.names=FALSE), res(rast))

  message('creating raster matrices for landscape variables...')
  matrices <- lapply(rast.reproj, .createMat)
  checks <- lapply(matrices, .checkMat)
  if ( sum(unlist(checks)) != length(checks) ) {
    stop("All the landscape layers should be numeric matricies, with values between -99.9 and 1")
  }

  message("== BARRIERS ==")
  message("reprojecting barriers...")
  barr.reproj <- lapply(barriers, .projectMe)

  message("cropping barriers to raster extent...")
  .cropMe <- function(x) {
    st_agr(x) = "constant"
    st_crop(x, unlist(ext))
  }
  barr.crop <- sapply(barr.reproj, .cropMe)

  message('generating barrier data')
  barr.df <- .generateBarriers(barr.crop, perms, progress=FALSE)

  mydata <- list(
    landscape_data = matrices,
    landscape_meta = envExt,
    barrier_data = barr.df
  )
  return(mydata)
}


# ******************************************************************************
# LOCAL FUNCTIONS (DO NOT EXPORT)
# ******************************************************************************
.createMat <- function(rast) {
  # transforming raster to matrix with size and extent info
  mat <- as.matrix(rast, wide=TRUE)
  mat[is.nan(mat)] <- -99.9
  mat
}
# CHECK INPUTS
.checkMat <- function(mat) {
  minCheck <- min(mat, na.rm=TRUE) >= -99.9
  maxCheck <- max(mat, na.rm=TRUE) <= 1
  flag = all(minCheck, maxCheck)
  flag
}
.segmentBarriers <- function(b, id=0, ...) {
  data <- st_segments(b, ...)
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
.generateBarriers <- function(barriers, p_list, ...) {
  # inputs are a shapefile of barrier data and a list of accompanying
  # permeabilities for each one. If all are the same for one barrier type,
  # then we only have to run it once for the whole dataset
  LP <- length(p_list)
  if (class(barriers)[1] == 'list') NB <- length(barriers)
  else NB <- nrow(barriers)
  if (LP != NB) {
    message("WARNING! length of permeability vector must be 1 or # of barriers")
    return()
  } else {
    if (LP == 1) p <- p_list
    df = NULL
    for (i in 1:NB) {
      if (LP > 1) p = p_list[i]
      b <- barriers[[i]]
      if (nrow(b) > 0) { # if cropping hasn't removed all data
        segs <- .segmentBarriers(b, i, progress=FALSE)
        mat  <- .createMatrix(segs, p)
        if (is.null(df)) df <- mat
        else df <- rbind(df, mat)
      }
    }
  }
  return(df)
}
