#' Transform plot ID to coordinates
#' 
#' The function transforms the plot ID (that consits of 6 digits) to
#' Spatialpoints in Swiss grid projection (old or new) or WGS projection. The
#' latter can be used to plot the points on a leaflet (see examples).
#' 
#' @param coordID Vector with plot coordinates
#' @param Z7 Logical 
#' @param projection What projection should be used for the coordinates. One of 
#' \enumerate{ 
#'     \item \code{CH-old}: Old Swiss-grid
#'     \item \code{CH-new}: New Swiss-grid
#'     \item \code{WGS}: WGS grid (e.g. to be used for leaflets, see example)}
#' @return Coordinates of sampling locations as \code{SpatialPoints} of package sp.
#' @export
#'
#' @examples
#' require(leaflet)
#' tmp <- coordID2coord(coordID = c(645260), Z7 = FALSE, projection = "WGS")
#' leaflet() %>% 
#'   addTiles() %>% 
#'   addCircles(data = tmp)
#' 
coordID2coord <- function(coordID, Z7 = TRUE, projection = "CH-old") {
  
  if (!is.character(projection)) stop("'Projection' should be one of 'CH-old', 'CH-new' or 'WGS'")
  proj <- match(projection, c("CH-old", "CH-new", "WGS"))
  if (is.na(proj)) stop("'Projection' should be one of 'CH-old', 'CH-new' or 'WGS'")
  
  if (Z7) {
    res <- data.frame(
      x = as.integer(substr(coordID, 1, 3))*1000 + 500,
      y = as.integer(substr(coordID, 4, 6))*1000 + 500)
  }
  if (!Z7) {
    res <- data.frame(
      x = as.integer(substr(coordID, 1, 3))*1000,
      y = as.integer(substr(coordID, 4, 6))*1000)
  }
  coordinates(res) <- ~ x + y
  proj4string(res) <- CRS("+init=epsg:21781") 
  spTransform(res, c(CRS("+init=epsg:21781"), CRS("+init=epsg:2056"), CRS("+init=epsg:4326"))[[proj]])
}
