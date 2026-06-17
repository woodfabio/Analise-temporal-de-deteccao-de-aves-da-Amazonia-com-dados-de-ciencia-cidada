# occupancy model
# to run with the whole dataset

# Clean up
rm(list=ls())
# load packages
library(terra)

# Load data ----
# Geography
ngstuff <- readRDS("DataObjects/neighborhoods_area_final.rds")
neighs <- ngstuff[[1]]
nneigh <- ngstuff[[2]]
hex <- unwrap(readRDS("DataObjects/area_final_grid.rds"))
Municip_wborders <- readRDS("DataObjects/MunicipalitiesMapWithBorders.rds")

# Observations

# cut off date
# choose one of the following dates

cutoffdate <- as.Date("2022-09-16")
cutoffdate <- as.Date("2020-01-01")

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

outputs_names_previous <- c("Outputs/temporal_analysis/Model_S_griseicapillus_output_1.rds",
                            "Outputs/temporal_analysis/Model_D_certhia_output_1.rds",
                            "Outputs/temporal_analysis/Model_P_marginatus_output_1.rds",
                            "Outputs/temporal_analysis/Model_F_colma_output_1.rds",
                            "Outputs/temporal_analysis/Model_A_ochrolaemus_output_1.rds",
                            "Outputs/temporal_analysis/Model_P_tricolor_output_1.rds",
                            "Outputs/temporal_analysis/Model_P_albifrons_output_1.rds",
                            "Outputs/temporal_analysis/Model_G_varia_output_1.rds",
                            "Outputs/temporal_analysis/Model_T_viridis_output_1.rds",
                            "Outputs/temporal_analysis/Model_T_ardesiacus_output_1.rds")

outputs_names_posterior <- c("Outputs/temporal_analysis/Model_S_griseicapillus_output_2.rds",
                             "Outputs/temporal_analysis/Model_D_certhia_output_2.rds",
                             "Outputs/temporal_analysis/Model_P_marginatus_output_2.rds",
                             "Outputs/temporal_analysis/Model_F_colma_output_2.rds",
                             "Outputs/temporal_analysis/Model_A_ochrolaemus_output_2.rds",
                             "Outputs/temporal_analysis/Model_P_tricolor_output_2.rds",
                             "Outputs/temporal_analysis/Model_P_albifrons_output_2.rds",
                             "Outputs/temporal_analysis/Model_G_varia_output_2.rds",
                             "Outputs/temporal_analysis/Model_T_viridis_output_2.rds",
                             "Outputs/temporal_analysis/Model_T_ardesiacus_output_2.rds")

# choose focal species from 1 to 10
focal_sps <- 1

eventTable <- readRDS(final_eventTable_names[focal_sps])
eventTable <- eventTable[!(is.na(eventTable$N.SPS.DETECTED) |     # exclude NA in n sps
                           is.na(eventTable$DURATION.MINUTES) |   # exclude NA in duration
                           eventTable$DURATION.MINUTES > 24*60),] # exclude long duration

# get first slice of data (based on cut off date)
eventTable <- eventTable[eventTable$OBSERVATION.DATE <= cutoffdate,]

# or get second slice of data (also based on cut off date)
eventTable <- eventTable[eventTable$OBSERVATION.DATE > cutoffdate,]

# Visualize observations
seen <- table(eventTable$GRID.CELL,eventTable$FOCAL.SPECIES.SEEN)
seen <- seen[which(seen[,2]>0),2]
cellids <- hex$id1
cres <- rep("white",length(nneigh))
cres[which(cellids %in% as.numeric(names(which(seen==1))))] <- "black"
cres[which(cellids %in% as.numeric(names(which(seen>1 & seen<5))))] <- "darkred"
cres[which(cellids %in% as.numeric(names(which(seen==5))))] <- "red"
cres[which(cellids %in% as.numeric(names(which(seen>5 & seen<=100))))] <- "orange"
cres[which(cellids %in% as.numeric(names(which(seen>100))))] <- "yellow"
plot(hex, col=cres, main="Number of RQ observations",border="darkgray",lwd=0.1)
add_legend("topright", border="darkgray",
           legend=c("0","1","2-4","5","6-100",">100"),
           fill=c("white","black","darkred","red","orange","yellow"),
           cex=0.95,bty="n")
plot(Municip_wborders,add=TRUE)

# Prepare observation vetting material -----
# get list with hexagon number and eBird sampling event identifier for 
# hexagons with 1, 2, 3, 4, and 5 RQ eBird sigthings
# this is meant to help validate observations from hexagons with only a few
eBird_ids <- list()
for(i in 1:5) {
  hexns <- as.numeric(names(seen[which(seen==i)]))
  eBird_ids[[i]] <- eventTable[which(eventTable$GRID.CELL %in% hexns & eventTable$FOCAL.SPECIES.SEEN),
                               #which(colnames(eventTable)) %in% c(1:3,8,17)]
                               which(colnames(eventTable) %in% c("GLOBAL.UNIQUE.IDENTIFIER",
                                                                  "COUNTY",
                                                                  "SAMPLING.EVENT.IDENTIFIER",
                                                                  "GRID.CELL"))]
  
}
names(eBird_ids) <- c("1 obs","2 obs","3 obs","4 obs","5 obs")

# Upload and setup Covariates ----
covs <- readRDS("DataObjects/covariates_area_final.rds")

stFC <- covs$StdForCover
rm(covs)

# Load libraries ----
library(nimble)
library(pracma)

# source functions

# Write NIMBLE model ----
Model2_RQmap <- nimbleCode({
  
  ## Priors
  
  # Beta0
  beta0 <- logit(mean_b0)
  mean_b0 ~ dunif(0,1)

  # Beta1
  beta1 ~ dnorm(0, sd=10)
  
  # Detection priors
  
  # Alpha0
  alpha0 <- logit(mean_a0)
  mean_a0 ~ dunif(0,1)
  
  # Alpha2
  alpha2 ~ dnorm(0,sd=10)
  
  ## Likelihood
  # Cell-specific occupancy
  for(i in 1:ncell) { # loop over grid cells
    
    z[i] ~ dbern(psi[i]) # True occupancy z at site i
    logit(psi[i]) <- beta0 + beta1*ForCover[i]
     
  }
  
  # list-specific detection of the focal species
  for(j in 1:nobs1) { # loop over observations
    
    logit(p1[j]) <- alpha0 + alpha2*tspent[j]
    zp1[j] <- p1[j]*z[cellID[j]]
    y1[j] ~ dbern(zp1[j])
    
  }
  
})

# Prepare detection covariates ----

# Standardize tspt
rawtspt <- eventTable$DURATION.MINUTES
mutspt <- mean(rawtspt)
sdtspt <- sd(rawtspt)
standardized_tspt <- (rawtspt - mutspt) / sdtspt

# Set up Nimble objects----
# constants
str(constantsM2 <- list(ncell = length(nneigh),
                        nobs1 = dim(eventTable)[1],
                        ForCover = stFC,
                        tspent = standardized_tspt,
                        cellID = eventTable$GRID.CELL))
# data 
str(dataM2 <- list(y1 = eventTable$FOCAL.SPECIES.SEEN)) # vector of 1 and 0
# and inits 
# z will be a vector of 1 and 0 with as many elements as sites
# z = 1 for those sites where the focal species was seen, z = 0 otherwise 
zobs <- rep(0,length(nneigh))
zobs[unique(eventTable$GRID.CELL[which(eventTable$FOCAL.SPECIES.SEEN==1)])] <- 1
str(initsM2 <- list(beta0 = rnorm(1),
                    beta1 = rnorm(1),
                    alpha0 = runif(1,min=0, max=3), 
                    alpha2 = runif(1,min=0, max=3),
                    z = zobs))


# Build and Compile Nimble model ----
tic() # create model object
M2Mod <- nimbleModel(code = Model2_RQmap, name = "M2_covsCAR",
                     constants = constantsM2,
                     data = dataM2,
                     inits = initsM2)
toc() # 8240 seconds (5324 in August 2025) (28834 in Outubro 2025)
# 3537 segundos no Simões em outubro de 2025 
# M0cMod$initializeInfo() returns [Note] All model variables are initialized.
# M0cMod$calculate() returns ...

tic() # configure an MCMC algorithm
conf <- configureMCMC(M2Mod, monitors = c("beta0","beta1",
                                          "alpha0","alpha2",
                                          "z", "psi","p1"),
                      enableWAIC = TRUE)
toc() # 27 seconds (or 12 in 10/25)

tic() # build an MCMC object
MCMC_M2 <- buildMCMC(conf)
toc() # 593 seconds (178 in 10/25)

tic() # compile model object to C++
cM2Mod <- compileNimble(M2Mod)
toc() # 43 seconds (or a bit less in 08/25 )

tic()
cMCMC_M2 <- compileNimble(MCMC_M2, project = cM2Mod)
toc() #  370 seconds (115 in 10/25   )


# Run the beast ----
tic()
output <- runMCMC(cMCMC_M2, nchains = 3, niter = 100000, nburnin = 20000,
                  thin = 100, summary = TRUE, WAIC = TRUE) # test
toc()

# Save output ----

# previous time slice
saveRDS(object=output,
        file=outputs_names_previous[focal_sps])

# posterior time slice
saveRDS(object=output,
        file=outputs_names_posterior[focal_sps])