#' Transform site ID to coordinates
#' 
#' The function transforms the site ID (that consits of 6 digits) to the
#' SpatialPolygons of squares (with a given sidelength) in Swiss grid projection
#' (old or new) or WGS projection. The latter can be used to plot the points on
#' a leaflet (see examples).
#' 
#' @param coordID Vector with site coordinates
#' @param projection What projection should be used for the coordinates. One of 
#' \enumerate{ 
#'     \item \code{CH-old}: Old Swiss-grid
#'     \item \code{CH-new}: New Swiss-grid
#'     \item \code{WGS}: WGS grid (e.g. to be used for leaflets, see example)}
#' @param sidelength The sidelength (in m) of the squares.
#' @return Shapefile (as \code{SpatialPolygons} of package sp) of the sampling squares.
#' @export
#'
#' @examples
#' require(leaflet)
#' tmp <- makesq(coordID = c(645260), projection = "WGS")
#' leaflet() %>% 
#'   addTiles() %>% 
#'   addPolygons(data = tmp)
#'
makesq <- function(coordID, projection = "CH-old", sidelength = 1000) {
  
  if (!is.character(projection)) stop("'Projection' should be one of 'CH-old', 'CH-new' or 'WGS'")
  proj <- match(projection, c("CH-old", "CH-new", "WGS"))
  if (is.na(proj)) stop("'Projection' should be one of 'CH-old', 'CH-new' or 'WGS'")
  
  p <- list()
  for(i in 1:length(coordID)) {
    x=as.integer(substr(coordID[i], 1, 3))*1000
    y=as.integer(substr(coordID[i], 4, 6))*1000
    sq <- data.frame(x=c(x, x+sidelength, x+sidelength, x, x),
                     y=c(y, y, y+sidelength,  y+sidelength, y))
    p[[i]] <- Polygons(list(Polygon(sq)), ID = i)
  }
  res <- SpatialPolygons(p, proj4string = CRS("+init=epsg:21781"))
  res <- spTransform(res, c(CRS("+init=epsg:21781"), CRS("+init=epsg:2056"), CRS("+init=epsg:4326"))[[proj]])
  return(res)
}





