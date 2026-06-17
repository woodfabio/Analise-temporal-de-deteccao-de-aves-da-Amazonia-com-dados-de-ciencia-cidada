## Spatialize eBird data

# clean up
rm(list=ls())

# load packages
library(terra)

# library(dplyr)
library(log4r)

# source functions
source("barracudar/SetUpLog.R")

# Load data ----
# spatial objects
area <- unwrap(readRDS("DataObjects/MapOutline_area_final.rds"))
areaWB <- readRDS("DataObjects/MunicipalitiesMapWithBorders.rds")
hex <- readRDS("DataObjects/area_final_grid.rds")

# neighborhood info
# get the neighbors for each site (neighs) and the number of neighbors per 
# site (nneigh) from "DataObjects/area_final_neighborhoods.rds"
bundle <- readRDS("DataObjects/neighborhoods_area_final.rds")
neighs <- bundle[[1]]
nneigh <- bundle[[2]]
rm(bundle)

# observations from eBird
sps_names <- c("Sittasomus griseicapillus", 
               "Dendrocolaptes certhia", 
               "Pachyramphus marginatus",
               "Formicarius colma",
               "Automolus ochrolaemus",
               "Perissocephalus tricolor",
               "Pithys albifrons",
               "Grallaria varia",
               "Trogon viridis",
               "Thamnomanes ardesiacus")

eventTable_names <- c("DataObjects/Processed_S_griseicapillus_eBird.rds",
                      "DataObjects/Processed_D_certhia_eBird.rds",
                      "DataObjects/Processed_P_marginatus_eBird.rds",
                      "DataObjects/Processed_F_colma_eBird.rds",
                      "DataObjects/Processed_A_ochrolaemus_eBird.rds",
                      "DataObjects/Processed_P_tricolor_eBird.rds",
                      "DataObjects/Processed_P_albifrons_eBird.rds",
                      "DataObjects/Processed_G_varia_eBird.rds",
                      "DataObjects/Processed_T_viridis_eBird.rds",
                      "DataObjects/Processed_T_ardesiacus_eBird.rds")

final_eventTable_names <- c("DataObjects/Final_Processed_S_griseicapillus_eBird.rds",
                            "DataObjects/Final_Processed_D_certhia_eBird.rds",
                            "DataObjects/Final_Processed_P_marginatus_eBird.rds",
                            "DataObjects/Final_Processed_F_colma_eBird.rds",
                            "DataObjects/Final_Processed_A_ochrolaemus_eBird.rds",
                            "DataObjects/Final_Processed_P_tricolor_eBird.rds",
                            "DataObjects/Final_Processed_P_albifrons_eBird.rds",
                            "DataObjects/Final_Processed_G_varia_eBird.rds",
                            "DataObjects/Final_Processed_T_viridis_eBird.rds",
                            "DataObjects/Final_Processed_T_ardesiacus_eBird.rds")

ebpoints_names <- c("DataObjects/eBpointsSgriseicapillus.rds",
                    "DataObjects/eBpointsDcerthia.rds",
                    "DataObjects/eBpointsPmarginatus.rds",
                    "DataObjects/eBpointsFcolma.rds",
                    "DataObjects/eBpointsAochrolaemus.rds",
                    "DataObjects/eBpointsPtricolor.rds",
                    "DataObjects/eBpointsPalbifrons.rds",
                    "DataObjects/eBpointsGvaria.rds",
                    "DataObjects/eBpointsTviridis.rds",
                    "DataObjects/eBpointsTardesiacus.rds")

# choose focal species from 1 to 10
focal_sps <- 1

eventTable <- readRDS(eventTable_names[focal_sps])

# Define focal species name
fsps <- sps_names[focal_sps]

# Protocol breakdown
fprot <- table(eventTable$PROTOCOL.NAME[which(eventTable$FOCAL.SPECIES.SEEN==1)])

cat("There are",sum(eventTable$FOCAL.SPECIES.SEEN),"records with the following
    protocols:",(if (length(fprot[which(names(fprot)=="Area")])==0) 0 else fprot[which(names(fprot)=="Area")]),"Area,",
                (if (length(fprot[which(names(fprot)=="Stationary")])==0) 0 else fprot[which(names(fprot)=="Stationary")]),"Stationary, and",
                (if (length(fprot[which(names(fprot)=="Traveling")])==0) 0 else fprot[which(names(fprot)=="Traveling")]),"Traveling")

# Spatialize species observations ----
eBird_points <- data.frame(x = eventTable$LONGITUDE,
                           y = eventTable$LATITUDE )

# transform to spatial points object
eBird_points <- vect(as.matrix(eBird_points))

# add sampling event identifier to table
eBird_points$SampleID <- eventTable$SAMPLING.EVENT.IDENTIFIER

# make sure coordinate reference systems match
crs(eBird_points) <- crs(area)

# this command may take a long time to run
system.time(eBird_points <- terra::crop(eBird_points,area))

# filter rows of eventTable according to cropped eBird_points
eventTable <- eventTable[which(eventTable$SAMPLING.EVENT.IDENTIFIER %in% eBird_points$SampleID),]

# add column of focal species observation to eBird_points
eBird_points$withSps <- eventTable$FOCAL.SPECIES.SEEN

# Visualize ----
plot(areaWB,col="lightgray")
points(eBird_points, cex=0.3)
points(eBird_points[which(eBird_points$withSps==1)],cex=0.5,col="red")
lines(area, col="white")

# Cross observations with hexagonal grid ----
system.time(hexeBird <- extract(hex,eBird_points))
eventTable$GRID.CELL <- hexeBird[,2]
eBird_points$GRID.CELL <- hexeBird[,2]

# Save spatialized eventTable ----
# Save evenTable with events that are within the universe of
# analysis (hex) and with an extra column indicating the corresponding
# hexagon in that same universe, or SpatVector hex
saveRDS(object=eventTable,
        file=final_eventTable_names[focal_sps]) 
saveRDS(object=eBird_points,
        file=ebpoints_names[focal_sps])
