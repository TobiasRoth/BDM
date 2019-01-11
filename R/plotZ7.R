#' Plot Z7 sites on a leaflet
#'
#' The function transforms the site ID (that consits of 6 digits) to spatial
#' data and plots them on a leaflet.
#'
#' @param coordID Vector with site coordinates
#' @param ortho Whether orthofotos should be used as background image.
#' @return A leaflet with ploted Z7 sites.
#' @export
#'
#' @examples
#' plotZ7(645260)
plotZ7 <- function(coordID, ortho = FALSE) {
  if(ortho) tilprov <- "Esri.WorldImagery"
  if(!ortho) tilprov <- "OpenStreetMap.CH"
  q <- makesq(coordID, projection = "WGS")
  map <- leaflet() %>% 
    addProviderTiles(tilprov) %>% 
    addPolygons(data = q, label = as.character(coordID))
  return(map)
}
