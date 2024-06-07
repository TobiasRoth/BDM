#' Transform plot ID to coordinates
#' 
#' The function transforms the plot ID (that consits of 6 digits) to a tibble of
#' class sf with a geometry that contains the points in Swiss grid projection
#' (old or new) or WGS projection. The latter can be used to plot the points on
#' a leaflet (see examples).
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
    res <- tibble(
      aID_STAO = coordID,
      x = paste0(str_sub(coordID, 1, 3), "500"),
      y = paste0(str_sub(coordID, 4, 6), "500")
    ) %>% 
      st_as_sf(coords = c("x", "y"), crs = 21781)
  }
  if (!Z7) {
    res <- tibble(
      aID_STAO = coordID,
      x = paste0(str_sub(coordID, 1, 3), "000"),
      y = paste0(str_sub(coordID, 4, 6), "000")
    ) %>% 
      st_as_sf(coords = c("x", "y"), crs = 21781)
  }

  st_transform(res, crs = c(21781, 2056, 4326)[proj])
}

