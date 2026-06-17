## Prepare universe of analysis for birds range
## For 3 municipalities only

# clean up
rm(list=ls())
# load packages
library(geodata)
library(terra)
library(sf)

# source functions
source("Functions/MakeGrid.R")
source("barracudar/SetUpLog.R")

## A. Outline study area without national borders ----

# read shapefiles of all municipalities of Amazonas state into a spatVector object
path <- "CleanedData/"
f <- c(paste(path,"AM_Municipios_2024.shp",sep=""))

# turn into list of SpatVectors
vv <- sapply(f, vect)

# unite all SpatVector elements of list into single SpatVector
area <- vect(vv)

# get only 3 municipalities
area_selected <- subset(area, NM_MUN %in% c("Manaus","Presidente Figueiredo","Rio Preto da Eva"), NSE=TRUE)

# unify r municipalities
area_united <- aggregate(area_selected)

# save map with borders
saveRDS(wrap(area_united),file="DataObjects/MunicipalitiesMapWithBorders.rds")

# aggregate geometries
area_final <- aggregate(area_united, dissolve=TRUE)
is.valid(area_final, TRUE, TRUE)

## B. Remove large inner water bodies ----

waterRast <- rast("CleanedData/2024_water_water_surface_1-3-5_768f8768-15ba-4d03-9fbe-7cc65bbcff24.tif")

# let's crop only the area of the 3 municipalities
waterRast <- crop(waterRast, area_final)
crs(waterRast) <- crs(area_final)
waterRast <- mask(waterRast, area_final)
plot(waterRast)

# the following line takes some time
water_bodies <- as.polygons(waterRast, dissolve = TRUE)
myshp <- writeVector(water_bodies, "water_bodies.shp", filetype = "ESRI Shapefile", overwrite = TRUE)

# check if it's alrigth (TRUE if alright)
myshp

plot(water_bodies)

# Make sure coordinate systems match
crs(water_bodies) <- crs(area_final)

# the following line may take a long time to run
area_final <- erase(area_final, water_bodies)
plot(area_final)

# save map outline
saveRDS(wrap(area_final),file="DataObjects/MapOutline_area_final.rds")

## C. Clean up ----

# Clean unnecessary large objects from workspace
rm(list = c("f","path"))


## D. Make the grid ----

# make hexagonal grid - ignore the 'x is invalid' message
system.time(hex <- make_grid(x=area_final, cell_diameter = 0.1))

# make sure the coordinate system is right
crs(hex) <- crs(area_final)

# visualize
plot(hex)

# save map grid
saveRDS(wrap(hex),file="DataObjects/area_final_grid.rds")

## E. List and count neighbors of each grid cell ----

neighs <- adjacent(x=hex, type="queen") # indices (col 2) of each (col 1) hexagon's neighbors
nNeigh <- table(neighs[,1])   # vector with number of neighbors per hexagon (name)

# show number of neighbors table
table(nNeigh)

# check that neighbor count makes sense by mapping number of neighbors as
# hexagon color
cores <- rep(0,length(nNeigh))
cores[which(nNeigh==6)] <- "pink"
cores[which(nNeigh==5)] <- "lightblue"
cores[which(nNeigh==4)] <- "red"
cores[which(nNeigh==3)] <- "blue"
cores[which(nNeigh==2)] <- "gray"
cores[which(nNeigh==1)] <- "black"

#plot(hex_area_final,col=cores)
plot(hex,col=cores)

# check that all hexagons in hex are accounted for (all well if TRUE)
#dim(hex_area_final)[1]==dim(nNeigh)
dim(hex)[1]==dim(nNeigh)

# bundle and save neighborhood information
bundle <- list(neighs=neighs,nneigh=nNeigh)
saveRDS(object=bundle,
        file="DataObjects/neighborhoods_area_final.rds") 

## F. Final clean-up ----
rm(list=c("bundle", "neighs", "nNeigh"))

# Forest Coverage

# load raster layer
AMrast <- rast("CleanedData/2024_coverage_lclu_1-3-5_205b2356-657e-4212-be77-c67edea889b6.tif")

# let's crop only the area of the 3 municipalities
forestRast <- crop(AMrast, area_final)
crs(forestRast) <- crs(area_final)
forestRast <- mask(forestRast, area_final)
plot(forestRast)

# check values in the Raster
terra::freq(AMrast)

# notice that almost all values are "value 3" (forest coverage)

# compute hexagon area
HexArea <- expanse(hex,unit="km")

# compute proportion of hexagon with forest coverage (FC)
meanFC <- extract(forestRast==3,hex,fun=mean,na.rm=TRUE)

# calculate forest coverage area as product of hexagon area and
# meanFC. Units are square kilometers
FCArea <- HexArea*meanFC[,2]

# Check that forest coverage numbers make sense
cores <- rep(0,length(FCArea))
cores[which(FCArea<20)] <- "white"
cores[which(FCArea>=20 & FCArea<40)] <- "#8CF36B"
cores[which(FCArea>=40 & FCArea<60)] <- "#4AF616"
cores[which(FCArea>=60 & FCArea<80)] <- "#0F9C08"
cores[which(FCArea>=80)] <- "#0F5701"
plot(hex, col=cores, main="Forest Coverage (km^2)",border="darkgray",lwd=0.1)
add_legend("topright",border="darkgray",
           legend=c("<20","21-40","41-60","61-80",">80"),
           fill=c("white","#8CF36B","#4AF616","#0F9C08","#0F5701"),
           cex=0.95,bty="n")
plot(area_final,add=TRUE)

# Standardize forest coverage Area
muFC <- mean(FCArea)
sdFC <- sd(FCArea)
stFC <- (FCArea - muFC) / sdFC

# check resolution of the forest coverage raster objects to be combined
res(forestRast)
#Adjust RESOLUTION to the raster with better resolution

# Make masks with the conditions that we want
FCpresent <- forestRast == 3

writeRaster(FCpresent,
            "Outputs/FCpresent.tif",
            overwrite = TRUE)

# bundle and save covariate information
bundle <- list(ForCover_km2=FCArea, StdForCover=stFC)

saveRDS(object=bundle,
        file="DataObjects/covariates_area_final.rds")
