## Import clean and combine eBird data
## Porto Alegre, February 28, 2026

# clean up
rm(list=ls())
# load packages
library(dplyr)
library(pracma)
library(ggplot2)
library(log4r)

# source functions
source("Functions/CleanEbird.R")
source("barracudar/SetUpLog.R")

# A. Set up filtering values ----
cwmt <- 4    # set the maximum route length as given multiple of cell widths

# this is a cell with multiplier
gcwdd <- 0.1 # grid cell width in decimal degrees, match with Code1
gcwkm <- round(gcwdd*100,digits=0) # approximate grid cell width in km
maxrl <- cwmt*gcwkm # obtain maximum route length 
maxar <- 100 # set maximum effort area in ha
mindate <- as.Date("2002-01-01")
maxdate <- as.Date("2026-12-31")
rm(list=c("cwmt","gcwdd","gcwkm"))

# B. Load and clean eBird data ----

# load eBird data from 3 municipalities
# (notice that this data is not available in this Github repository since it
# needs permission from eBird to be used, instead, in this repository we have
# the already processed data in files for each species, so you can run line
# 58 to 103 to choose the focal species and then run starting from line 122)

manaus <- read.delim("OriginalData/eBird raw data/all_sps/ebd_BR-AM-038_200201_202612_smp_relJan-2026/ebd_BR-AM-038_200201_202612_smp_relJan-2026.txt", quote ="")
pf <- read.delim("OriginalData/eBird raw data/all_sps/ebd_BR-AM-048_200201_202612_smp_relJan-2026/ebd_BR-AM-048_200201_202612_smp_relJan-2026.txt", quote ="")
rpe <- read.delim("OriginalData/eBird raw data/all_sps/ebd_BR-AM-049_200201_202612_smp_relJan-2026/ebd_BR-AM-049_200201_202612_smp_relJan-2026.txt", quote ="")

datesmanaus <- as.Date(manaus$OBSERVATION.DATE) # filter dates
rowkeepasmanaus <- which(datesmanaus>=mindate & datesmanaus<=maxdate)

datespf <- as.Date(pf$OBSERVATION.DATE) # filter dates
rowkeepaspf<- which(datespf>=mindate & datespf<=maxdate)

datesrpe <- as.Date(rpe$OBSERVATION.DATE) # filter dates
rowkeepasrpe <- which(datesrpe>=mindate & datesrpe<=maxdate)

manaus <- manaus[rowkeepasmanaus,]
pf <- pf[rowkeepaspf,]
rpe <- rpe[rowkeepasrpe,]

# bind 3 municipalities toghether. We will use the separated municipalities objects later for plotting.
ebas <- rbind(manaus,pf,rpe)

# objects to choose which species to analyze

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

files_names <- c("DataObjects/ProcEbDataSg.rds",
                 "DataObjects/ProcEbDataDc.rds",
                 "DataObjects/ProcEbDataPm.rds",
                 "DataObjects/ProcEbDataFc.rds",
                 "DataObjects/ProcEbDataAo.rds",
                 "DataObjects/ProcEbDataPt.rds",
                 "DataObjects/ProcEbDataPa.rds",
                 "DataObjects/ProcEbDataGv.rds",
                 "DataObjects/ProcEbDataTv.rds",
                 "DataObjects/ProcEbDataTa.rds")

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

plots_names <- c("Plots/FigureA_TossedObs_S_griseicapillus.pdf",
                 "Plots/FigureA_TossedObs_D_certhia.pdf",
                 "Plots/FigureA_TossedObs_P_marginatus.pdf",
                 "Plots/FigureA_TossedObs_F_colma.pdf",
                 "Plots/FigureA_TossedObs_A_ochrolaemus.pdf",
                 "Plots/FigureA_TossedObs_P_tricolor.pdf",
                 "Plots/FigureA_TossedObs_P_albifrons.pdf",
                 "Plots/FigureA_TossedObs_G_varia.pdf",
                 "Plots/FigureA_TossedObs_T_viridis.pdf",
                 "Plots/FigureA_TossedObs_T_ardesiacus.pdf")

# choose focal species from 1 to 10
focal_sps <- 1

eb <- clean_ebird(ebdata=ebas,fsps=sps_names[focal_sps],maxrl=maxrl, maxar=maxar)
ov <- eb[[2]] # extract observation vector, with number of records at each stage
eb <- eb[[1]] # extract eBird data

obskpt <- ov

# save to Data Objects folder

bundle <- list(eBdata=eb,obskept=ov)
saveRDS(object=bundle,
        file=files_names[focal_sps]) 
rm(bundle)


# C. Load Filtered eBird data and process it down to table of events ----

# load eBird data processed above
restored_bundle <- readRDS(file=files_names[focal_sps])
eBdata <- restored_bundle[[1]]
obskpt <- restored_bundle[[2]]
rm(restored_bundle)

# Build event table and sort it
# get row number for the first instance of each sampling event
firstrow <- which(!duplicated(eBdata$SAMPLING.EVENT.IDENTIFIER))

# build table
eventTable <- eBdata[firstrow,]

# sort
eventTable <- eventTable[order(eventTable$SAMPLING.EVENT.IDENTIFIER),]

# do all observations of each event have the same coordinates?
# use dplyr tricks to handle the huge data frame and pracma for timing
tic()
result <- eventTable %>%
  group_by(SAMPLING.EVENT.IDENTIFIER) %>%
  summarize(SameCoords = n() == 1 || all(LATITUDE == first(LATITUDE) & LONGITUDE == first(LONGITUDE)))

# Check if all SameCoords are TRUE
if(all(result$SameCoords)) {
  cat("All values of SAMPLING.EVENT.IDENTIFIER have the same coordinates.")
} else {
  cat("Some values of SAMPLING.EVENT.IDENTIFIER have different coordinates.")
}
toc()
rm(result)

# FOCAL SPECIES SIGHTINGS

# get rows of eBird data that have observations of the focal sps
rowsfsps <- which(eBdata$SCIENTIFIC.NAME==sps_names[focal_sps])

# get events where focal species was seen
evtsfsps <- eBdata$SAMPLING.EVENT.IDENTIFIER[rowsfsps]

# GATHER INFO PER EVENT USING VECTORS WITH AS MANY ROWS AS EVENTS IN THE TABLE

# build vector that says whether the focal species was seen in the corresponding
# sampling event
fsps1 <- as.numeric(eventTable$SAMPLING.EVENT.IDENTIFIER %in% evtsfsps)

# build vector with number of records per event
nrecs <- table(eBdata$SAMPLING.EVENT.IDENTIFIER) # checked that event order ok

# build vector with number of species per event (to check if number of records
# and number of species is the same)
nsps <- aggregate(SCIENTIFIC.NAME ~ SAMPLING.EVENT.IDENTIFIER, eBdata, function(x) length(unique(x)))$SCIENTIFIC.NAME
cat("There are", sum(nrecs!=nsps),"lists with duplicated species records.")

# ADD NEW COLUMNS TO THE EVENT TABLE DATAFRAME
eventTable$FOCAL.SPECIES.SEEN <- fsps1
eventTable$N.SPS.DETECTED <- nsps

# CLEAN OUT COLUMNS THAT NO LONGER ADD ANY INFORMATION
removeColumns <- c("SCIENTIFIC.NAME",
                   "EXOTIC.CODE",
                   "ALL.SPECIES.REPORTED",
                   "APPROVED")
eventTable <- eventTable[,-which(colnames(eventTable) %in% removeColumns)]

# SAVE EVENTS TABLE DATAFRAME IN RDATA
saveRDS(object=eventTable,
        file=eventTable_names[focal_sps])

# D. Count and plot tossed observations ----

## Barplot of eBird observation filtering
nf <- 6 # number of filters
np <- 3 # number of municipalities

namp <- c("Manaus","Presidente Figueiredo","Rio Preto da Eva") # Municipalities abbreviation
namf <- c("Incomplete","Wrong Protocol","Multiple input",
          "Long route","Large area","Good to use") # filters

ovmanaus <- clean_ebird(ebdata=manaus,fsps=sps_names[focal_sps],maxrl=maxrl, maxar=maxar)[[2]]
ovpf <- clean_ebird(ebdata=pf,fsps=sps_names[focal_sps],maxrl=maxrl, maxar=maxar)[[2]]
ovrpe <- clean_ebird(ebdata=rpe,fsps=sps_names[focal_sps],maxrl=maxrl, maxar=maxar)[[2]]

obskpt <- rbind(ovmanaus,ovpf,ovrpe)
row.names(obskpt) <- c("Manaus", "Presidente Figueiredo", "Rio Preto da Eva")

Bd <- matrix(0, nrow = nf, ncol = np, byrow = T,
             dimnames = list(namf,namp))

Bd[,1] <- 100*obskpt["Manaus",-c(1:2)] / sum(obskpt["Manaus",3])
Bd[,2] <- 100*obskpt["Presidente Figueiredo",-c(1:2)] / sum(obskpt["Presidente Figueiredo",3])
Bd[,3] <- 100*obskpt["Rio Preto da Eva",-c(1:2)] / sum(obskpt["Rio Preto da Eva",3])

Bd<-Bd-rbind(Bd[2:nf,],rep(0,3))
Bd<-Bd[nrow(Bd):1,]

# plot while saving
pdf(file=plots_names[focal_sps],
    width = 8.8, height = 6.5)
p<-barplot(height = (Bd),
           xlab = "Municipalities",
           ylab = "% eBird Records", ylim = c(0,115),
           legend.text = rownames(Bd),
           args.legend = list(x = "bottomright",inset=c(0.05,0.05), bg="white"),
           col = c("#4A9BDC","#32C7A9","#81BB42","#E9BF35","#FE801A","#DF2E28")
)
text(x = p, y = 100 + 6, labels = as.character(obskpt[,1]))
dev.off()

# final clean-up

rm(list=c("Bd","eBdata","eventTable","obskpt",
          "evtsfsps","firstrow","fsps1","maxrl","maxar","namf",
          "namp","nf","np","nrecs","nsps","p","rowsfsps",
          "maxdate","mindate"))
