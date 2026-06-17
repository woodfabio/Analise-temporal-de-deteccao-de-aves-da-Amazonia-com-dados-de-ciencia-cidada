# --------------------------------------
# FUNCTION make_grid
# required packages: terra sf
# description: draws grid of hexagons with given diameter or area
#              within the boundaries of a SpatVector geometry 
# inputs: SpatVector, hexagon diameter (distance between oposite sides 
#         in degrees) or hexagon area
# outputs: SpatVector grid of hexagons 
########################################
make_grid <- function(x, cell_diameter, cell_area) {
  
  if (missing(cell_diameter)) {
    if (missing(cell_area)) {
      stop("Must provide cell_diameter or cell_area")
    } else {
      cell_diameter <- sqrt(2 * cell_area / sqrt(3))
    }
  }
  
  # create hexagonal grid over extent of x
  # cell diameter is the distance between hexagon sides in 
  # units of the x coordinate system (degrees)
  g <- st_make_grid(x,cellsize=cell_diameter,square=FALSE)
  # convert sfc_polygon to spatVector
  g <- vect(g)
  # keep only what intersects with x
  g <- terra::intersect(g,x)
  # insert a column with hexagon number as id1
  g$id1 <- seq(1:length(g))
  
  return(g)
} # end of function make_grid
# --------------------------------------
# make_grid()
