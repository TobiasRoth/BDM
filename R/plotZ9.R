#' Plot Z9 sites on a leaflet
#'
#' The function transforms the site ID (that consits of 6 digits) to spatial
#' data and plots them on a leaflet.
#'
#' @param coordID Vector with site coordinates
#' @param ortho Whether orthofotos should be used as background image.
#' @return A leaflet with ploted Z9 sites.
#' @export
#'
#' @examples
#' plotZ9(645260)
plotZ9 <- function(coordID,  ortho = TRUE) {
  if(ortho) tilprov <- "Esri.WorldImagery"
  if(!ortho) tilprov <- "OpenStreetMap.CH"
  q <- gBuffer(coordID2coord(coordID, Z7 = FALSE), width = sqrt(10 / pi), quadsegs = 100)
  q <- spTransform(q, CRS("+init=epsg:4326"))
  map <- leaflet() %>% 
    addProviderTiles(tilprov) %>% 
    addPolygons(data = q, color = "red", weight = 3, opacity = 1, fillColor = "red", label = as.character(coordID))
  return(map)
}
