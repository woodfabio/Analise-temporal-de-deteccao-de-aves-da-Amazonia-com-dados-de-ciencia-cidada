## Script that takes output from code 4's NIMBLE run and draws a map
## of psi and/or z

# Clean up
rm(list=ls())

# Load libraries ----
library(terra)
library(basetheme)
library(basicMCMCplots)
library(nimble)
library(MCMCvis)

# define species objects ----

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

# Load data ----
out1 <- readRDS(outputs_names_previous[focal_sps])
out2 <- readRDS(outputs_names_posterior[focal_sps])

hex <- unwrap(readRDS("DataObjects/area_final_grid.rds"))
borders <- readRDS("DataObjects/MunicipalitiesMapWithBorders.rds")

# Traceplots ----
chainsPlot(out1$samples, var=c("beta0",
                               "alpha0",
                               "alpha1",
                               "alpha2"))

chainsPlot(out2$samples, var=c("beta0",
                               "alpha0",
                               "alpha1",
                               "alpha2"))

# Compare prior and posterior of some parameters ----
SimsPriors <- cbind(beta1=rnorm(n=10000,mean=0,sd=100),
                    alpha0=logit(runif(n=10000,min=0,max=1)),
                    alpha2=runif(n=10000,min=0,max=3))

MCMCtrace(out1$samples,
          params=c('beta1','alpha0','alpha2'),
          ISB=FALSE, exact=TRUE, priors=SimsPriors,Rhat=TRUE, n.eff=FALSE,
          post_zm=FALSE,pdf=FALSE,plot=TRUE)

MCMCtrace(out2$samples,
          params=c('beta1','alpha0','alpha2'),
          ISB=FALSE, exact=TRUE, priors=SimsPriors,Rhat=TRUE, n.eff=FALSE,
          post_zm=FALSE,pdf=FALSE,plot=TRUE)

# check detection values ----
alpha0_samples_1 <- c(out1$samples$chain1[,"alpha0"], out1$samples$chain2[,"alpha0"], out1$samples$chain3[,"alpha0"])
alpha0_samples_2 <- c(out2$samples$chain1[,"alpha0"], out2$samples$chain2[,"alpha0"], out2$samples$chain3[,"alpha0"])

mean(plogis(alpha0_samples_1))
quantile(plogis(alpha0_samples_1), c(0.025, 0.975))

mean(plogis(alpha0_samples_2))
quantile(plogis(alpha0_samples_2), c(0.025, 0.975))

# Maps ----

# out1
par(mfrow = c(1, 2))
p1 <- out1$summary$all.chains[5:4152,"Mean"]
psi <- out1$summary$all.chains[4153:4628,"Mean"]
z <- out1$summary$all.chains[4629:5104,"Mean"]
plot(hex,col=num2col(x=psi))
plot(hex,col=num2col(z))

colrs <- rep(0,length(psi))
colrs[which(psi<0.2)] <- "white"
colrs[which(psi>=0.2 & psi<0.4)] <- "yellow"
colrs[which(psi>=0.4 & psi<0.6)] <- "orange"
colrs[which(psi>=0.6 & psi<0.8)] <- "red"
colrs[which(psi>=0.8)] <- "darkred"

plot(hex,col=colrs,main="Occupancy (psi)",border="darkgray",)
add_legend("topright",border="darkgray",
           legend=c("< 0.2","0.2-0.4","0.4-0.6","0.6-0.8",">0.8"),
           fill=c("white","yellow","orange","red","darkred"),
           cex=0.7,bty="n")
plot(borders,add=TRUE,border="darkgray",)

colrs <- rep(0,length(z))
colrs[which(z<0.2)] <- "white"
colrs[which(z>=0.2 & z<0.4)] <- "yellow"
colrs[which(z>=0.4 & z<0.6)] <- "orange"
colrs[which(z>=0.6 & z<0.8)] <- "red"
colrs[which(z>=0.8)] <- "darkred"

plot(hex,col=colrs,main="Mean presence (z)",border="darkgray")
add_legend("topright",
           legend=c("< 0.2","0.2-0.4","0.4-0.6","0.6-0.8",">0.8"),
           fill=c("white","yellow","orange","red","darkred"),
           cex=0.7,bty="n")
plot(borders,add=TRUE,border="darkgray",)

colrs <- rep(0,length(psi))
colrs[which(p1<0.2)] <- "white"
colrs[which(p1>=0.2 & p1<0.4)] <- "yellow"
colrs[which(p1>=0.4 & p1<0.6)] <- "orange"
colrs[which(p1>=0.6 & p1<0.8)] <- "red"
colrs[which(p1>=0.8)] <- "darkred"

plot(hex,col=colrs,main="Detection (p)",border="darkgray",)
add_legend("topright",border="darkgray",
           legend=c("< 0.2","0.2-0.4","0.4-0.6","0.6-0.8",">0.8"),
           fill=c("white","yellow","orange","red","darkred"),
           cex=0.7,bty="n")
plot(borders,add=TRUE,border="darkgray",)


# out2
par(mfrow = c(1, 2))
p1 <- out2$summary$all.chains[5:4152,"Mean"]
psi <- out2$summary$all.chains[4153:4628,"Mean"]
z <- out2$summary$all.chains[4629:5104,"Mean"]
plot(hex,col=num2col(x=psi))
plot(hex,col=num2col(z))

colrs <- rep(0,length(psi))
colrs[which(psi<0.2)] <- "white"
colrs[which(psi>=0.2 & psi<0.4)] <- "yellow"
colrs[which(psi>=0.4 & psi<0.6)] <- "orange"
colrs[which(psi>=0.6 & psi<0.8)] <- "red"
colrs[which(psi>=0.8)] <- "darkred"

plot(hex,col=colrs,main="Occupancy (psi)",border="darkgray",)
add_legend("topright",border="darkgray",
           legend=c("< 0.2","0.2-0.4","0.4-0.6","0.6-0.8",">0.8"),
           fill=c("white","yellow","orange","red","darkred"),
           cex=0.7,bty="n")
plot(borders,add=TRUE,border="darkgray",)

colrs <- rep(0,length(psi))
colrs[which(z<0.2)] <- "white"
colrs[which(z>=0.2 & z<0.4)] <- "yellow"
colrs[which(z>=0.4 & z<0.6)] <- "orange"
colrs[which(z>=0.6 & z<0.8)] <- "red"
colrs[which(z>=0.8)] <- "darkred"

plot(hex,col=colrs,main="Mean presence (z)",border="darkgray")
add_legend("topright",
           legend=c("< 0.2","0.2-0.4","0.4-0.6","0.6-0.8",">0.8"),
           fill=c("white","yellow","orange","red","darkred"),
           cex=0.7,bty="n")
plot(borders,add=TRUE,border="darkgray",)

colrs <- rep(0,length(p1))
colrs[which(psi<0.2)] <- "white"
colrs[which(psi>=0.2 & psi<0.4)] <- "yellow"
colrs[which(psi>=0.4 & psi<0.6)] <- "orange"
colrs[which(psi>=0.6 & psi<0.8)] <- "red"
colrs[which(psi>=0.8)] <- "darkred"

plot(hex,col=colrs,main="Detection (p)",border="darkgray",)
add_legend("topright",border="darkgray",
           legend=c("< 0.2","0.2-0.4","0.4-0.6","0.6-0.8",">0.8"),
           fill=c("white","yellow","orange","red","darkred"),
           cex=0.7,bty="n")
plot(borders,add=TRUE,border="darkgray",)


# Plots ----

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

# detection

# cut off date on 2022

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

outs_p_2022 <- data.frame(matrix(ncol = 6, nrow = 10))
cols <- c("Mean 1", 
          "Mean 2", 
          "Quantile low 1", 
          "Quantile up 1", 
          "Quantile low 2", 
          "Quantile up 2")
colnames(outs_p_2022) <- cols

means1 <- rep(NA, 10)
means2 <- rep(NA, 10)
quantileslow1 <- rep(NA, 10)
quantilesup1 <- rep(NA, 10)
quantileslow2 <- rep(NA, 10)
quantilesup2 <- rep(NA, 10)
for (i in 1:10) {
  output <- readRDS(file = outputs_names_previous[i])
  alpha0_samples_1 <- c(output$samples$chain1[,"alpha0"], 
                        output$samples$chain2[,"alpha0"], 
                        output$samples$chain3[,"alpha0"])
  means1[i] <- mean(plogis(alpha0_samples_1))
  quantileslow1[i] <- quantile(plogis(alpha0_samples_1), 0.025)
  quantilesup1[i] <- quantile(plogis(alpha0_samples_1), 0.975)
  
  output <- readRDS(file = outputs_names_posterior[i])
  alpha0_samples_2 <- c(output$samples$chain1[,"alpha0"], 
                        output$samples$chain2[,"alpha0"], 
                        output$samples$chain3[,"alpha0"])
  means2[i] <- mean(plogis(alpha0_samples_2))
  quantileslow2[i] <- quantile(plogis(alpha0_samples_2), 0.025)
  quantilesup2[i] <- quantile(plogis(alpha0_samples_2), 0.975)
}

outs_p_2022$`Mean 1` <- means1
outs_p_2022$`Mean 2` <- means2
outs_p_2022$`Quantile low 1` <- quantileslow1
outs_p_2022$`Quantile up 1` <- quantilesup1
outs_p_2022$`Quantile low 2` <- quantileslow2
outs_p_2022$`Quantile up 2` <- quantilesup2

par(mar = c(5, 15, 4, 2) + 0.1) 
plot(outs_p_2022$`Mean 1`, c(1:10)+0.1,
     xlab = "Média do intercepto da detecção (alpha0) com data de corte em 16/09/2022",
     xlim = c(0, 0.25),
     yaxt = "n",
     ylab = "",
     pch = 1)
title(ylab = "Espécies", line = 12)
points(outs_p_2022$`Mean 2`, c(1:10)-0.1,
     yaxt = "n",
     ylab = "",
     pch = 19)
axis(2, at = c(1:10), labels = sps_names, las = 1, font.axis = 3)
abline(h = c(1:10), lty = 2, col = "gray")

segments(x0 = outs_p_2022$`Quantile low 1`, y0 = c(1:10)+0.1, 
         x1 = outs_p_2022$`Quantile up 1`, y1 = c(1:10)+0.1)

segments(x0 = outs_p_2022$`Quantile low 2`, y0 = c(1:10)-0.1, 
         x1 = outs_p_2022$`Quantile up 2`, y1 = c(1:10)-0.1)

legend("bottomright",
       legend=c("antes", "depois"),
       pch=c(1, 19))  

# cut off date on 2020

outputs_names_previous_2020 <- c("Outputs/temporal_analysis/2020_cutoffdate/Model_S_griseicapillus_output_2020_1.rds",
                            "Outputs/temporal_analysis/2020_cutoffdate/Model_D_certhia_output_2020_1.rds",
                            "Outputs/temporal_analysis/2020_cutoffdate/Model_P_marginatus_output_2020_1.rds",
                            "Outputs/temporal_analysis/2020_cutoffdate/Model_F_colma_output_2020_1.rds",
                            "Outputs/temporal_analysis/2020_cutoffdate/Model_A_ochrolaemus_output_2020_1.rds",
                            "Outputs/temporal_analysis/2020_cutoffdate/Model_P_tricolor_output_2020_1.rds",
                            "Outputs/temporal_analysis/2020_cutoffdate/Model_P_albifrons_output_2020_1.rds",
                            "Outputs/temporal_analysis/2020_cutoffdate/Model_G_varia_output_2020_1.rds",
                            "Outputs/temporal_analysis/2020_cutoffdate/Model_T_viridis_output_2020_1.rds",
                            "Outputs/temporal_analysis/2020_cutoffdate/Model_T_ardesiacus_output_2020_1.rds")

outputs_names_posterior_2020 <- c("Outputs/temporal_analysis/2020_cutoffdate/Model_S_griseicapillus_output_2020_2.rds",
                             "Outputs/temporal_analysis/2020_cutoffdate/Model_D_certhia_output_2020_2.rds",
                             "Outputs/temporal_analysis/2020_cutoffdate/Model_P_marginatus_output_2020_2.rds",
                             "Outputs/temporal_analysis/2020_cutoffdate/Model_F_colma_output_2020_2.rds",
                             "Outputs/temporal_analysis/2020_cutoffdate/Model_A_ochrolaemus_output_2020_2.rds",
                             "Outputs/temporal_analysis/2020_cutoffdate/Model_P_tricolor_output_2020_2.rds",
                             "Outputs/temporal_analysis/2020_cutoffdate/Model_P_albifrons_output_2020_2.rds",
                             "Outputs/temporal_analysis/2020_cutoffdate/Model_G_varia_output_2020_2.rds",
                             "Outputs/temporal_analysis/2020_cutoffdate/Model_T_viridis_output_2020_2.rds",
                             "Outputs/temporal_analysis/2020_cutoffdate/Model_T_ardesiacus_output_2020_2.rds")

outs_p_2020 <- data.frame(matrix(ncol = 6, nrow = 10))
cols <- c("Mean 1", 
          "Mean 2", 
          "Quantile low 1", 
          "Quantile up 1", 
          "Quantile low 2", 
          "Quantile up 2")
colnames(outs_p_2020) <- cols

means1 <- rep(NA, 10)
means2 <- rep(NA, 10)
quantileslow1 <- rep(NA, 10)
quantilesup1 <- rep(NA, 10)
quantileslow2 <- rep(NA, 10)
quantilesup2 <- rep(NA, 10)
for (i in 1:10) {
  output <- readRDS(file = outputs_names_previous_2020[i])
  alpha0_samples_1 <- c(output$samples$chain1[,"alpha0"], 
                        output$samples$chain2[,"alpha0"], 
                        output$samples$chain3[,"alpha0"])
  means1[i] <- mean(plogis(alpha0_samples_1))
  quantileslow1[i] <- quantile(plogis(alpha0_samples_1), 0.025)
  quantilesup1[i] <- quantile(plogis(alpha0_samples_1), 0.975)
  
  output <- readRDS(file = outputs_names_posterior_2020[i])
  alpha0_samples_2 <- c(output$samples$chain1[,"alpha0"], 
                        output$samples$chain2[,"alpha0"], 
                        output$samples$chain3[,"alpha0"])
  means2[i] <- mean(plogis(alpha0_samples_2))
  quantileslow2[i] <- quantile(plogis(alpha0_samples_2), 0.025)
  quantilesup2[i] <- quantile(plogis(alpha0_samples_2), 0.975)
}

outs_p_2020$`Mean 1` <- means1
outs_p_2020$`Mean 2` <- means2
outs_p_2020$`Quantile low 1` <- quantileslow1
outs_p_2020$`Quantile up 1` <- quantilesup1
outs_p_2020$`Quantile low 2` <- quantileslow2
outs_p_2020$`Quantile up 2` <- quantilesup2

par(mar = c(5, 15, 4, 2) + 0.1) 
plot(outs_p_2020$`Mean 1`, c(1:10)+0.1,
     xlab = "Média do intercepto da detecção (alpha0) com data de corte em 01/01/2020",
     xlim = c(0, 0.25),
     yaxt = "n",
     ylab = "",
     pch = 1)
title(ylab = "Espécies", line = 12)
points(outs_p_2020$`Mean 2`, c(1:10)-0.1,
       yaxt = "n",
       ylab = "",
       pch = 19)
axis(2, at = c(1:10), labels = sps_names, las = 1, font.axis = 3)
abline(h = c(1:10), lty = 2, col = "gray")

segments(x0 = outs_p_2020$`Quantile low 1`, y0 = c(1:10)+0.1, 
         x1 = outs_p_2020$`Quantile up 1`, y1 = c(1:10)+0.1)

segments(x0 = outs_p_2020$`Quantile low 2`, y0 = c(1:10)-0.1, 
         x1 = outs_p_2020$`Quantile up 2`, y1 = c(1:10)-0.1)

legend("bottomright",
       legend=c("antes", "depois"),
       pch=c(1, 19))


# occurence

# cut off date on 2022

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

outs_psi_2022 <- data.frame(matrix(ncol = 6, nrow = 10))
cols <- c("Mean 1", 
          "Mean 2", 
          "Quantile low 1", 
          "Quantile up 1", 
          "Quantile low 2", 
          "Quantile up 2")
colnames(outs_psi_2022) <- cols

means1 <- rep(NA, 10)
means2 <- rep(NA, 10)
quantileslow1 <- rep(NA, 10)
quantilesup1 <- rep(NA, 10)
quantileslow2 <- rep(NA, 10)
quantilesup2 <- rep(NA, 10)
for (i in 1:10) {
  output <- readRDS(file = outputs_names_previous[i])
  beta0_samples_1 <- c(output$samples$chain1[,"beta0"], 
                        output$samples$chain2[,"beta0"], 
                        output$samples$chain3[,"beta0"])
  means1[i] <- mean(plogis(beta0_samples_1))
  quantileslow1[i] <- quantile(plogis(beta0_samples_1), 0.025)
  quantilesup1[i] <- quantile(plogis(beta0_samples_1), 0.975)
  
  output <- readRDS(file = outputs_names_posterior[i])
  beta0_samples_2 <- c(output$samples$chain1[,"beta0"], 
                        output$samples$chain2[,"beta0"], 
                        output$samples$chain3[,"beta0"])
  means2[i] <- mean(plogis(beta0_samples_2))
  quantileslow2[i] <- quantile(plogis(beta0_samples_2), 0.025)
  quantilesup2[i] <- quantile(plogis(beta0_samples_2), 0.975)
}

outs_psi_2022$`Mean 1` <- means1
outs_psi_2022$`Mean 2` <- means2
outs_psi_2022$`Quantile low 1` <- quantileslow1
outs_psi_2022$`Quantile up 1` <- quantilesup1
outs_psi_2022$`Quantile low 2` <- quantileslow2
outs_psi_2022$`Quantile up 2` <- quantilesup2

par(mar = c(5, 15, 4, 2) + 0.1) 
plot(outs_psi_2022$`Mean 1`, c(1:10)+0.1,
     xlab = "Média do intercepto da detecção (beta0) com data de corte em 16/09/2022",
     xlim = c(0, 1),
     yaxt = "n",
     ylab = "",
     pch = 1)
title(ylab = "Espécies", line = 12)
points(outs_psi_2022$`Mean 2`, c(1:10)-0.1,
       yaxt = "n",
       ylab = "",
       pch = 19)
axis(2, at = c(1:10), labels = sps_names, las = 1, font.axis = 3)
abline(h = c(1:10), lty = 2, col = "gray")

segments(x0 = outs_psi_2022$`Quantile low 1`, y0 = c(1:10)+0.1, 
         x1 = outs_psi_2022$`Quantile up 1`, y1 = c(1:10)+0.1)

segments(x0 = outs_psi_2022$`Quantile low 2`, y0 = c(1:10)-0.1, 
         x1 = outs_psi_2022$`Quantile up 2`, y1 = c(1:10)-0.1)

legend("bottomleft",
       legend=c("antes", "depois"),
       pch=c(1, 19))  

# cut off date on 2020

outputs_names_previous_2020 <- c("Outputs/temporal_analysis/2020_cutoffdate/Model_S_griseicapillus_output_2020_1.rds",
                                 "Outputs/temporal_analysis/2020_cutoffdate/Model_D_certhia_output_2020_1.rds",
                                 "Outputs/temporal_analysis/2020_cutoffdate/Model_P_marginatus_output_2020_1.rds",
                                 "Outputs/temporal_analysis/2020_cutoffdate/Model_F_colma_output_2020_1.rds",
                                 "Outputs/temporal_analysis/2020_cutoffdate/Model_A_ochrolaemus_output_2020_1.rds",
                                 "Outputs/temporal_analysis/2020_cutoffdate/Model_P_tricolor_output_2020_1.rds",
                                 "Outputs/temporal_analysis/2020_cutoffdate/Model_P_albifrons_output_2020_1.rds",
                                 "Outputs/temporal_analysis/2020_cutoffdate/Model_G_varia_output_2020_1.rds",
                                 "Outputs/temporal_analysis/2020_cutoffdate/Model_T_viridis_output_2020_1.rds",
                                 "Outputs/temporal_analysis/2020_cutoffdate/Model_T_ardesiacus_output_2020_1.rds")

outputs_names_posterior_2020 <- c("Outputs/temporal_analysis/2020_cutoffdate/Model_S_griseicapillus_output_2020_2.rds",
                                  "Outputs/temporal_analysis/2020_cutoffdate/Model_D_certhia_output_2020_2.rds",
                                  "Outputs/temporal_analysis/2020_cutoffdate/Model_P_marginatus_output_2020_2.rds",
                                  "Outputs/temporal_analysis/2020_cutoffdate/Model_F_colma_output_2020_2.rds",
                                  "Outputs/temporal_analysis/2020_cutoffdate/Model_A_ochrolaemus_output_2020_2.rds",
                                  "Outputs/temporal_analysis/2020_cutoffdate/Model_P_tricolor_output_2020_2.rds",
                                  "Outputs/temporal_analysis/2020_cutoffdate/Model_P_albifrons_output_2020_2.rds",
                                  "Outputs/temporal_analysis/2020_cutoffdate/Model_G_varia_output_2020_2.rds",
                                  "Outputs/temporal_analysis/2020_cutoffdate/Model_T_viridis_output_2020_2.rds",
                                  "Outputs/temporal_analysis/2020_cutoffdate/Model_T_ardesiacus_output_2020_2.rds")

outs_psi_2020 <- data.frame(matrix(ncol = 6, nrow = 10))
cols <- c("Mean 1", 
          "Mean 2", 
          "Quantile low 1", 
          "Quantile up 1", 
          "Quantile low 2", 
          "Quantile up 2")
colnames(outs_psi_2020) <- cols

means1 <- rep(NA, 10)
means2 <- rep(NA, 10)
quantileslow1 <- rep(NA, 10)
quantilesup1 <- rep(NA, 10)
quantileslow2 <- rep(NA, 10)
quantilesup2 <- rep(NA, 10)
for (i in 1:10) {
  output <- readRDS(file = outputs_names_previous_2020[i])
  beta0_samples_1 <- c(output$samples$chain1[,"beta0"], 
                       output$samples$chain2[,"beta0"], 
                       output$samples$chain3[,"beta0"])
  means1[i] <- mean(plogis(beta0_samples_1))
  quantileslow1[i] <- quantile(plogis(beta0_samples_1), 0.025)
  quantilesup1[i] <- quantile(plogis(beta0_samples_1), 0.975)
  
  output <- readRDS(file = outputs_names_posterior_2020[i])
  beta0_samples_2 <- c(output$samples$chain1[,"beta0"], 
                       output$samples$chain2[,"beta0"], 
                       output$samples$chain3[,"beta0"])
  means2[i] <- mean(plogis(beta0_samples_2))
  quantileslow2[i] <- quantile(plogis(beta0_samples_2), 0.025)
  quantilesup2[i] <- quantile(plogis(beta0_samples_2), 0.975)
}

outs_psi_2020$`Mean 1` <- means1
outs_psi_2020$`Mean 2` <- means2
outs_psi_2020$`Quantile low 1` <- quantileslow1
outs_psi_2020$`Quantile up 1` <- quantilesup1
outs_psi_2020$`Quantile low 2` <- quantileslow2
outs_psi_2020$`Quantile up 2` <- quantilesup2

par(mar = c(5, 15, 4, 2) + 0.1) 
plot(outs_psi_2020$`Mean 1`, c(1:10)+0.1,
     xlab = "Média do intercepto da detecção (beta0) com data de corte em 01/01/2020",
     xlim = c(0, 1),
     yaxt = "n",
     ylab = "",
     pch = 1)
title(ylab = "Espécies", line = 12)
points(outs_psi_2020$`Mean 2`, c(1:10)-0.1,
       yaxt = "n",
       ylab = "",
       pch = 19)
axis(2, at = c(1:10), labels = sps_names, las = 1, font.axis = 3)
abline(h = c(1:10), lty = 2, col = "gray")

segments(x0 = outs_psi_2020$`Quantile low 1`, y0 = c(1:10)+0.1, 
         x1 = outs_psi_2020$`Quantile up 1`, y1 = c(1:10)+0.1)

segments(x0 = outs_psi_2020$`Quantile low 2`, y0 = c(1:10)-0.1, 
         x1 = outs_psi_2020$`Quantile up 2`, y1 = c(1:10)-0.1)

legend("bottomleft",
       legend=c("antes", "depois"),
       pch=c(1, 19))  

