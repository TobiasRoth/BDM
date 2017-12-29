#' Plot Z7 sites on a leaflet
#'
#' The function transforms the site ID (that consits of 6 digits) to spatial
#' data and plots them on a leaflet.
#'
#' @param coordID Vector with site coordinates
#' @return A leaflet with ploted Z7 or Z9 sites.
#' @export
#'
#' @examples
#' require(leaflet)
#' plotZ7(645260)
plotZ7 <- function(coordID, Z7 = TRUE) {
  q <- makesq(coordID, projection = "WGS")
  map <- leaflet() %>% 
    addTiles() %>% 
    addPolygons(data = q)
  return(map)
}
