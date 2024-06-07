#' Transform site ID to a shape of a square
#'
#' The function transforms the site ID (that consits of 6 digits and points to
#' the SW corner of the square) to a tibble with class sf and a polygon as a
#' geometry. The polygon is a square with a given sidelength in meter. The
#' projection is eiter a the Swiss projection (old or new) or WGS projection.
#' The latter can be used to plot the points on a leaflet (see examples).
#'
#' @param coordID Vector or tibble with site site-IDs. If a tibble, the column
#'   with the site-IDs should be named as "aID_STAO".
#' @param projection What projection should be used for the coordinates. One of
#' \enumerate{
#'     \item \code{CH-old}: Old Swiss-grid
#'     \item \code{CH-new}: New Swiss-grid
#'     \item \code{WGS}: WGS grid (e.g. to be used for leaflets, see example)}
#' @param sidelength The sidelength (in m) of the squares.
#' @return Tibble of class \code{sf} with a geometry that describes the sampling
#'   squares.
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
  
  # Match and control projection
  if (!is.character(projection)) stop("'Projection' should be one of 'CH-old', 'CH-new' or 'WGS'")
  proj <- match(projection, c("CH-old", "CH-new", "WGS"))
  if (is.na(proj)) stop("'Projection' should be one of 'CH-old', 'CH-new' or 'WGS'")
  
  # Calculate polygons from vector of site-IDs
  if(is.vector(coordID)) {
    res <- tibble(
      aID_STAO = coordID,
      x = paste0(str_sub(coordID, 1, 3), "500"),
      y = paste0(str_sub(coordID, 4, 6), "500")
    ) %>% 
      st_as_sf(coords = c("x", "y"), crs = 21781) %>% 
      st_buffer(dist = sidelength / 2, endCapStyle = "SQUARE") %>% 
      st_transform(crs = c(21781, 2056, 4326)[proj])
  }
  
  # Add polygons to table
  if(!is.vector(coordID)) {
    tmp <- coordID %>% 
      mutate(
        x = paste0(str_sub(aID_STAO, 1, 3), "500"),
        y = paste0(str_sub(aID_STAO, 4, 6), "500")
      ) %>% 
      st_as_sf(coords = c("x", "y"), crs = 21781) %>% 
      st_buffer(dist = sidelength / 2, endCapStyle = "SQUARE") %>% 
      st_transform(crs = c(21781, 2056, 4326)[proj])
    st_geometry(coordID) <- st_geometry(tmp)
    res <- coordID
  }

  res
}




