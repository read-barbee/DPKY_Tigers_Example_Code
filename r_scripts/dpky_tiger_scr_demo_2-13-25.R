#### DPKY Tiger SCR models NIMBLE ####

# Author: Read Barbee

# Date:2025-02-13 

# Purpose: Estimate annual density and abundance of tigers in DPKY


################################ Libraries #################################
library(tidyverse)
library(nimbleSCR)
library(doParallel)
library(coda)
library(ggmcmc)
library(MCMCvis)
library(sf)
library(terra)
library(loo)
library(mapview)


################################ Helper functions #################################

split_logProb <- function(mcmc_list) {
  # Function to separate logProb_y columns
  split_list <- map(mcmc_list, ~ {
    col_names <- attr(.x, "dimnames")[[2]]
    
    # Identify columns that contain "logProb_y"
    logProb_cols <- grepl("logProb_y", col_names)
    
    # Separate into two matrices: one for logProb_y and one for others
    logProb_chain <- .x[, logProb_cols]
    non_logProb_chain <- .x[, !logProb_cols]
    
    # Update dimnames for both subsets
    attr(logProb_chain, "dimnames")[[2]] <- col_names[logProb_cols]
    attr(non_logProb_chain, "dimnames")[[2]] <- col_names[!logProb_cols]
    
    list(logProb = logProb_chain, non_logProb = non_logProb_chain)
  })
  
  # Extract the logProb chains and non-logProb chains
  logProb_mcmc_list <- map(split_list, "logProb")
  non_logProb_mcmc_list <- map(split_list, "non_logProb")
  
  # Restore the class to "mcmc.list" for both
  class(logProb_mcmc_list) <- "mcmc.list"
  class(non_logProb_mcmc_list) <- "mcmc.list"
  
  return(list(logProb_mcmc_list = logProb_mcmc_list,
              non_logProb_mcmc_list = non_logProb_mcmc_list))
}

#select dataframes for models with covs of interest for model averaging
select_dataframes <- function(df_list, column_name) {
  new_list <- map(df_list, function(df) {
    if (column_name %in% colnames(df)) {
      return(df)
    } else {
      return(NULL)
    }})
  
  output <- Filter(Negate(is.null), new_list)#new_list[!is.null(new_list)]
  return(output)
}

################################ User defined parameters #################################

set.seed(777)

#initial values for parameters 
sigma_init = 40
g0_init = 0.05
psi_init = 0.5

pi_init <- 0.5 # probability of augmented individual's sex

#number of individuals to augment population with (possible but never detected)
naug <- 500


#mcmc parameters
chains <- 3
burnin <- 1000
iter = 10000 


###############################################################################

#1. Import data

###############################################################################
#year
year <- "2020"

grid <- "og_grid"

#capture data (simulated)
ca <- read_csv("data/secr/2020_CA_demo.csv")

#trap data
td <- read_csv("data/secr/2020_TD.csv")

#protected area polygon
DPKY <- st_read("data/secr/gis_layers/TLNP-HABprj.shp")

#habitat mask from SECR analysis
DPKY_mask <- read_csv("data/secr/DPKY_mask_70_day_clip_new_9-15-24.csv")


###############################################################################

#2. Create the state space matrix with habitat mask

###############################################################################

#convert the habitat mask to an sf object
DPKY_mask_sf <- st_as_sf(DPKY_mask, coords = c(x = "x", y ="y"), crs = 32648)


#add column denoting whether grid points are inside or outside the protected area polygon. If they are all inside, this still creates the necessary vector of 1s for rasterization.
habitat_mask <- DPKY_mask_sf %>%
  mutate(mask = as.integer(st_intersects(DPKY_mask_sf, DPKY, sparse = FALSE)))

#Make a template raster to transfer mask values to
tmp_rast <- rast(ext(habitat_mask), resolution = c(500, 500), crs = "EPSG:32648", vals = 999)

#rasterize the mask points
habitat_rast <- habitat_mask %>% 
  rasterize(y = tmp_rast, field = "mask") 


#convert the raster to a matrix
habitat.mx <- habitat_rast %>% 
  as.matrix(wide = TRUE)

#replace NA values with 0
habitat.mx[is.na(habitat.mx)] <- 0

#Crop the DPKY polygon to the area of the mask
DPKY_crop <- st_crop(DPKY, DPKY_mask_sf) #mask_hull

#Calculate the area of inference for density estimates in km2
area_dpky_crop <- as.numeric(st_area(DPKY_crop)) /1e6

#Define the limits of the state space based on the habitat matrix
y.max <- dim(habitat.mx)[1]
x.max <- dim(habitat.mx)[2]

xlim <- c(0, x.max)
ylim <- c(0, y.max)


###############################################################################

#3. Create the trap location matrix. Index number of traps and trap effort

###############################################################################

#convert the detector locations to an sf object for plotting
trap_coords <- st_as_sf(td %>% select(X,Y), coords = c(x = "X", y = "Y"), crs = 32648)

plot(habitat_rast)
points(trap_coords)


#retrieve the raw detector coordinates
tp_coords <-  st_coordinates(trap_coords)

ext <- ext(habitat_rast)  # Extent of the raster: xmin, xmax, ymin, ymax
res <- res(habitat_rast)

# Calculate scaling factors
x_scale_tp <- (ncol(habitat_rast) - 1) / (ext[2] - ext[1])
y_scale_tp <- (nrow(habitat_rast) - 1) / (ext[4] - ext[3])

# Scale the coordinates to match the raster matrix
scaled_coords_tp <- cbind(
  x = (tp_coords[, "X"] - ext[1]) * x_scale_tp + 1,
  y = (tp_coords[, "Y"] - ext[3]) * y_scale_tp + 1
)

# Extract values from the matrix using scaled coordinates
# Ensure the coordinates are within bounds of the matrix dimensions
valid_coords_tp <- scaled_coords_tp[
  scaled_coords_tp[, "x"] >= 1 &
    scaled_coords_tp[, "x"] <= ncol(habitat.mx) &
    scaled_coords_tp[, "y"] >= 1 &
    scaled_coords_tp[, "y"] <= nrow(habitat.mx), 
]



X <-  scaled_coords_tp #valid_coords_tp

# K = the number of days in which each detector was active
K <- td %>% 
  select(-c(`#STATION`, X, Y)) %>% 
  rowSums() 

#J = the number of detectors
J <- nrow(td)


#plot the original habitat matrix, the scaled matrix, and the scaled coordinates to make sure everything lines up
plot(rast(habitat.mx))
points(X, col = "red")




###############################################################################

#4. Create the augmented capture history matrix

###############################################################################

#rename capture data columns
colnames(ca) <- c("Session", "ID", "Occasion", "Station", "Sex")

#index the number of unique captures
ncaps <- nrow(ca)

#index the number of detected individuals
ninds <- length(unique(ca$ID))


###### ADD THIS
# Get station number from td
station_num <- td %>%
  rename(Station = `#STATION`) %>% 
  select(Station) %>%
  mutate(Station_num = 1:n())


#condense capture data to the number of times each individual was detected at each station. Add numeric ids for numbers and stations for matrix indexing.
capssum <- ca %>%
  group_by(ID) %>%
  mutate(IDnum = cur_group_id()) %>% #add individual indexes
  as.data.frame() %>%
  left_join(station_num, by="Station") %>%
  group_by(IDnum, Station_num) %>%
  mutate(Count = n()) %>% #count the rows for each combination of individual and station 
  slice(1) %>% #remove duplicate rows caused by independent counts for each detection occasion
  as.data.frame() %>%  
  select(-Occasion)


### Create the raw detection history matrix ###

#  nrow = number of identified individuals
#  ncol = number of trap locations

y <- matrix(0, nrow = ninds, ncol = J)
for(i in 1:nrow(capssum)){
  y[capssum$IDnum[i], capssum$Station_num[i]] <- capssum$Count[i]
}

#Add a row of zeros to the detection matrix for every augmented individual 
yaug <- rbind(y, matrix(0, naug, J))

#record the augmented population size
M <- nrow(yaug)


#create sex covariate
sex <- capssum %>%
  group_by(IDnum) %>%
  slice(1) %>%
  select(IDnum, Sex) %>%
  mutate(sexcov = case_when(Sex == "M" ~ 0,
                            Sex == "F" ~ 1,
                            .default = NA))


###############################################################################

#5. Generate initial values for activity centers to be estimated for each individual

###############################################################################

#sample random points as initial activity centers for each individual
r_extent <- st_as_sfc(st_bbox(habitat_rast))

cells_with_1 <- which(values(habitat_rast) == 1)
r_coords <- xyFromCell(habitat_rast, cells_with_1)

r_sampled_coords <- r_coords[sample(nrow(r_coords), nrow(yaug)), ]

sxy_init <- st_as_sf(as.data.frame(r_sampled_coords), coords = c("x", "y"), crs = crs(habitat_rast))


sxy_init_coords <- sxy_init %>% 
  st_coordinates()


# Calculate scaling factors
x_scale_sx <- (ncol(habitat_rast) - 1) / (ext[2] - ext[1])
y_scale_sx <- (nrow(habitat_rast) - 1) / (ext[4] - ext[3])

# Scale the coordinates to match the raster matrix
scaled_coords_sx <- cbind(
  x = (sxy_init_coords[, "X"] - ext[1]) * x_scale_sx + 1,
  y = (sxy_init_coords[, "Y"] - ext[3]) * y_scale_sx + 1
)

# Extract values from the matrix using scaled coordinates
# Ensure the coordinates are within bounds of the matrix dimensions
valid_coords_sx <- scaled_coords_sx[
  scaled_coords_sx[, "x"] >= 1 &
    scaled_coords_sx[, "x"] <= ncol(habitat.mx) &
    scaled_coords_sx[, "y"] >= 1 &
    scaled_coords_sx[, "y"] <= nrow(habitat.mx), 
]

# Extract the matrix values at the valid coordinates
#values_at_coords <- habitat.mx[cbind(valid_coords[, "y"], valid_coords[, "x"])]

plot(rast(habitat.mx))
points(X, col ="red")
points(valid_coords_sx, col = "blue")


###############################################################################

#6. Define model inputs: Data, constants, and initial values

###############################################################################

#define vector to constrain activity centers to suitable habitat
ok <- rep(1, nrow(yaug)) 

#define vector of individual's latent state in the population (included or not)
zi <- c(rep(1, nrow(y)), rep(0, naug)) 

#define vector of sexes for augmented population
sex_dat = c(sex$sexcov, rep(NA, naug))

sex_init = c( rep(NA, length(sex$sexcov)), rbinom(naug, 1, 0.5))


constants <- list(M = M, 
                  K = K, 
                  J = J, 
                  area = area_dpky_crop, 
                  y.max = y.max, 
                  x.max = x.max
)

data <- list(y = yaug, 
             X = X, 
             xlim = xlim, 
             ylim = ylim, 
             habitat.mx = habitat.mx,
             ones = ok,
             sex = sex_dat)

#modify dHabitat function to return very small log likelihood values for ACs outside of habitat mask instead of -Inf. This will help in calculating WAIC. Not working
dHabitatMaskrb <- nimbleFunction(
  run = function (x = double(0), 
                  s = double(1), 
                  xmax = double(0), 
                  xmin = double(0), 
                  ymax = double(0), 
                  ymin = double(0), 
                  habitatMask = double(2), 
                  log = integer(0)) 
  {
    returnType(double(0))
    if (s[1] < xmin) 
      return(1e-10)
    if (s[1] > xmax) 
      return(1e-10)
    if (s[2] < ymin) 
      return(1e-10)
    if (s[2] > ymax) 
      return(1e-10)
    test <- 1 - (habitatMask[trunc(s[2]) + 1, trunc(s[1]) + 1] == 0)
    test <- max(test, 1e-10)  # Avoid 0 probabilities
    
    if (log) 
      return(log(test))
    else return(test)
  }
)


rHabitatMaskrb <- nimbleFunction(
  run = function (n = double(0), 
                  s = double(1), 
                  xmax = double(0), 
                  xmin = double(0), 
                  ymax = double(0), 
                  ymin = double(0), 
                  habitatMask = double(2)) 
  {
    returnType(double(0))
    if (s[1] < xmin) 
      return(1e-10)
    if (s[1] > xmax) 
      return(1e-10)
    if (s[2] < ymin) 
      return(1e-10)
    if (s[2] > ymax) 
      return(1e-10)
    if (habitatMask[trunc(s[2]) + 1, trunc(s[1]) + 1] == 0) 
      return(1e-10)
    else return(1)
  }
)

registerDistributions("dHabitatMaskrb")



###############################################################################

#7a. M0: g0 ~ 1, sigma ~ 1

###############################################################################

inits_m0 <- list(s = valid_coords_sx,
                 sigma = sigma_init, 
                 g0 = g0_init, 
                 psi = psi_init, 
                 z = zi)


m0 <- nimbleCode({ 
  
  # Priors 
  psi ~ dbeta(1,1)
  g0 ~ dunif(0.01,1)
  sigma ~ dunif(0.01,100)
  
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    
    ## habitat constraint: 1 if the randomly generated AC falls in suitable habitat, -Inf if it does not.
    ones[i] ~ dHabitatMaskrb( s = s[i,1:2],
                              xmin = xlim[1],
                              xmax = xlim[2],
                              ymin = ylim[1],
                              ymax = ylim[2],
                              habitat = habitat.mx[1:y.max,1:x.max])
    
    for(j in 1:J){
      dist[i,j] <- sqrt((s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2)
      g[i,j] <- g0 * exp(-dist[i,j]^2 / (2 * sigma^2)) * z[i]
      y[i,j] ~ dpois(g[i,j] * K[j]) # K is different per trap so indexed
      #loglik[i,j]<-dnorm(y[i,j], mu, sd=sigma, log=1)
    } 
  }
  N <- sum(z[])
  D <- N/area
  #D_100 <- D * 100
}) 


Rmodel_m0 <- nimbleModel(m0, constants, data, inits_m0)

conf_m0 <- configureMCMC(Rmodel_m0, 
                         monitors = c("D", "N", "g0", "sigma", "psi", "logProb_y"), 
                         enableWAIC = TRUE,
                         print = FALSE)
conf_m0$removeSamplers("s")
ACnodes_m0 <- paste0("s[", 1:M, ", 1:2]")
for(node in ACnodes_m0) {
  conf_m0$addSampler(target = node,
                     type = "RW_block",
                     control = list(adaptScaleOnly = TRUE),
                     silent = TRUE)
}
Rmcmc_m0 <- buildMCMC(conf_m0)

Cmodel_m0 <- compileNimble(Rmodel_m0)
Cmcmc_m0 <- compileNimble(Rmcmc_m0, project = Rmodel_m0)
MCMC_runtime_m0 <- system.time(
  samples_m0 <- runMCMC(Cmcmc_m0, 
                        nchains = chains, 
                        nburnin=burnin, 
                        niter = iter, 
                        WAIC = TRUE,
                        samplesAsCodaMCMC = TRUE)
)

mcmc.out_m0 <- coda::mcmc.list(samples_m0$samples)

log_prob_split_m0 <- split_logProb(mcmc.out_m0)

main_list_m0 <- log_prob_split_m0$non_logProb_mcmc_list
log_list_m0 <- log_prob_split_m0$logProb_mcmc_list

ll_mat_m0 <- as.matrix(log_list_m0)

loo_m0 <- loo::loo(ll_mat_m0, cores = 3) #takes awhile

output_m0 <- MCMCsummary(main_list_m0) %>% rownames_to_column("parameter")

#object for plotting
S_m0 <- ggs(main_list_m0) %>% 
  mutate(Parameter = fct_relevel(Parameter, "N", "D", "g0", "sigma", "psi"))

waic_tab_m0 <- tibble(model = "m0",
                      waic = samples_m0$WAIC$WAIC,
                      pwaic = samples_m0$WAIC$pWAIC,
                      lppd = samples_m0$WAIC$lppd,
                      elpd_loo = loo_m0$estimates["elpd_loo",1],
                      p_loo =  loo_m0$estimates["p_loo",1],
                      looic = loo_m0$estimates["looic",1])


#output the posterior summaries
#write_csv(output_m0, paste0("m0_70_day_baseline_ext_",year, "_results.csv"))

###############################################################################

#7b. M1: g0 ~ SEX, sigma ~ 1

###############################################################################

inits_m1 <- list(s = valid_coords_sx,
                 sigma = sigma_init, 
                 g0 = rep(g0_init, 2),
                 psi = psi_init, 
                 pi = pi_init,
                 z = zi,
                 sex = sex_init)



m1 <- nimbleCode({ 
  
  # Priors 
  psi ~ dbeta(1,1)
  pi ~ dbeta(1,1) 
  
  g0[1] ~ dunif(0.01,1)
  g0[2] ~ dunif(0.01,1)
  sigma ~ dunif(0.01,100)
  
  for(i in 1:M){
    z[i] ~ dbern(psi)
    sex[i] ~ dbern(pi)
    sex2[i] <- sex[i] + 1
    
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    
    ## habitat constraint: 1 if the randomly generated AC falls in suitable habitat, -Inf if it does not.
    ones[i] ~ dHabitatMaskrb( s = s[i,1:2],
                              xmin = xlim[1],
                              xmax = xlim[2],
                              ymin = ylim[1],
                              ymax = ylim[2],
                              habitat = habitat.mx[1:y.max,1:x.max])
    
    for(j in 1:J){
      dist[i,j] <- sqrt((s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2)
      g[i,j] <- g0[sex2[i]] * exp(-dist[i,j]^2 / (2 * sigma^2)) * z[i]
      y[i,j] ~ dpois(g[i,j] * K[j]) # K is different per trap so indexed
    } 
  }
  
  N <- sum(z[])
  D <- N/area
  #D_100 <- D * 100
}) 

Rmodel_m1 <- nimbleModel(m1, constants, data, inits_m1)

conf_m1 <- configureMCMC(Rmodel_m1, 
                         monitors = c("D", "N", "g0", "sigma", "psi", "pi", "logProb_y"),
                         enableWAIC = TRUE,
                         print = FALSE)
conf_m1$removeSamplers("s")
ACnodes_m1 <- paste0("s[", 1:M, ", 1:2]")
for(node in ACnodes_m1) {
  conf_m1$addSampler(target = node,
                     type = "RW_block",
                     control = list(adaptScaleOnly = TRUE),
                     silent = TRUE)
}
Rmcmc_m1 <- buildMCMC(conf_m1)

Cmodel_m1 <- compileNimble(Rmodel_m1)
Cmcmc_m1 <- compileNimble(Rmcmc_m1, project = Rmodel_m1)
MCMC_runtime_m1 <- system.time(
  samples_m1 <- runMCMC(Cmcmc_m1, 
                        nchains = chains, 
                        nburnin=burnin, 
                        niter = iter, 
                        WAIC = TRUE,
                        samplesAsCodaMCMC = TRUE)
)


mcmc.out_m1 <- coda::mcmc.list(samples_m1$samples)

log_prob_split_m1 <- split_logProb(mcmc.out_m1)

main_list_m1 <- log_prob_split_m1$non_logProb_mcmc_list
log_list_m1 <- log_prob_split_m1$logProb_mcmc_list

ll_mat_m1 <- as.matrix(log_list_m1)

loo_m1 <- loo::loo(ll_mat_m1, cores = 3) #takes awhile, about 2.5 min

output_m1 <- MCMCsummary(main_list_m1) %>% rownames_to_column("parameter")

#object for plotting
S_m1 <- ggs(main_list_m1) %>% 
  mutate(Parameter = fct_relevel(Parameter, "N", "D", "g0[1]", "g0[2]", "sigma", "psi", "pi")) %>% 
  mutate(Parameter = fct_recode(Parameter, "g0[Male]" = "g0[1]", "g0[Female]" = "g0[2]"))

waic_tab_m1 <- tibble(model = "m1",
                      waic = samples_m1$WAIC$WAIC,
                      pwaic = samples_m1$WAIC$pWAIC,
                      lppd = samples_m1$WAIC$lppd,
                      elpd_loo = loo_m1$estimates["elpd_loo",1],
                      p_loo =  loo_m1$estimates["p_loo",1],
                      looic = loo_m1$estimates["looic",1])


#output the posterior summaries
#write_csv(output_m1, paste0("m1_70_day_baseline_ext_", year, "_results.csv"))


###############################################################################

#7c. M2: g0 ~ 1, sigma ~ SEX

###############################################################################

inits_m2 <- list(s = valid_coords_sx,
                 sigma = rep(sigma_init, 2), 
                 g0 = g0_init,
                 psi = psi_init, 
                 pi = pi_init,
                 z = zi,
                 sex = sex_init)

m2 <- nimbleCode({ 
  
  # Priors 
  psi ~ dbeta(1,1)
  pi ~ dbeta(1,1) 
  
  g0 ~ dunif(0.01,1)
  sigma[1] ~ dunif(0.01,100)
  sigma[2] ~ dunif(0.01,100)
  
  for(i in 1:M){
    z[i] ~ dbern(psi)
    sex[i] ~ dbern(pi)
    sex2[i] <- sex[i] + 1
    
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    
    ## habitat constraint: 1 if the randomly generated AC falls in suitable habitat, -Inf if it does not.
    ones[i] ~ dHabitatMaskrb( s = s[i,1:2],
                              xmin = xlim[1],
                              xmax = xlim[2],
                              ymin = ylim[1],
                              ymax = ylim[2],
                              habitat = habitat.mx[1:y.max,1:x.max])
    
    for(j in 1:J){
      dist[i,j] <- sqrt((s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2)
      g[i,j] <- g0 * exp(-dist[i,j]^2 / (2 * sigma[sex2[i]]^2)) * z[i]
      y[i,j] ~ dpois(g[i,j] * K[j]) # K is different per trap so indexed
    } 
  }
  
  N <- sum(z[])
  D <- N/area
  #D_100 <- D * 100
}) 

Rmodel_m2 <- nimbleModel(m2, constants, data, inits_m2)

conf_m2 <- configureMCMC(Rmodel_m2, 
                         monitors = c("D", "N", "g0", "sigma", "psi", "pi", "logProb_y"),
                         enableWAIC = TRUE,
                         print = FALSE)
conf_m2$removeSamplers("s")
ACnodes_m2 <- paste0("s[", 1:M, ", 1:2]")
for(node in ACnodes_m2) {
  conf_m2$addSampler(target = node,
                     type = "RW_block",
                     control = list(adaptScaleOnly = TRUE),
                     silent = TRUE)
}
Rmcmc_m2 <- buildMCMC(conf_m2)

Cmodel_m2 <- compileNimble(Rmodel_m2)
Cmcmc_m2 <- compileNimble(Rmcmc_m2, project = Rmodel_m2)
MCMC_runtime_m2 <- system.time(
  samples_m2 <- runMCMC(Cmcmc_m2, 
                        nchains = chains, 
                        nburnin=burnin, 
                        niter = iter, 
                        WAIC = TRUE,
                        samplesAsCodaMCMC = TRUE)
)


mcmc.out_m2 <- coda::mcmc.list(samples_m2$samples)

log_prob_split_m2 <- split_logProb(mcmc.out_m2)

main_list_m2 <- log_prob_split_m2$non_logProb_mcmc_list
log_list_m2 <- log_prob_split_m2$logProb_mcmc_list

ll_mat_m2 <- as.matrix(log_list_m2)

loo_m2 <- loo::loo(ll_mat_m2, cores = 3) #takes awhile, about 2.5 min

output_m2 <- MCMCsummary(main_list_m2) %>% rownames_to_column("parameter")

#object for plotting
S_m2 <- ggs(main_list_m2) %>% 
  mutate(Parameter = fct_relevel(Parameter, "N", "D", "g0", "sigma[1]", "sigma[2]", "psi", "pi")) %>% 
  mutate(Parameter = fct_recode(Parameter, "sigma[Male]" = "sigma[1]", "sigma[Female]" = "sigma[2]"))

waic_tab_m2 <- tibble(model = "m2",
                      waic = samples_m2$WAIC$WAIC,
                      pwaic = samples_m2$WAIC$pWAIC,
                      lppd = samples_m2$WAIC$lppd,
                      elpd_loo = loo_m2$estimates["elpd_loo",1],
                      p_loo =  loo_m2$estimates["p_loo",1],
                      looic = loo_m2$estimates["looic",1])


#output the posterior summaries
#write_csv(output_m2, paste0("m2_70_day_baseline_ext_", year, "_results.csv"))


###############################################################################

#7d. M3: g0 ~ SEX, sigma ~ SEX

###############################################################################

inits_m3 <- list(s = valid_coords_sx,
                 sigma = rep(sigma_init, 2), 
                 g0 = rep(g0_init, 2),
                 psi = psi_init, 
                 pi = pi_init,
                 z = zi,
                 sex = sex_init)


m3 <- nimbleCode({ 
  
  # Priors 
  psi ~ dbeta(1,1)
  pi ~ dbeta(1,1) 
  
  g0[1] ~ dunif(0.01,1)
  g0[2] ~ dunif(0.01,1)
  sigma[1] ~ dunif(0.01,100)
  sigma[2] ~ dunif(0.01,100)
  
  for(i in 1:M){
    z[i] ~ dbern(psi)
    sex[i] ~ dbern(pi)
    sex2[i] <- sex[i] + 1
    
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    
    ## habitat constraint: 1 if the randomly generated AC falls in suitable habitat, -Inf if it does not.
    ones[i] ~ dHabitatMaskrb( s = s[i,1:2],
                              xmin = xlim[1],
                              xmax = xlim[2],
                              ymin = ylim[1],
                              ymax = ylim[2],
                              habitat = habitat.mx[1:y.max,1:x.max])
    
    for(j in 1:J){
      dist[i,j] <- sqrt((s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2)
      g[i,j] <- g0[sex2[i]] * exp(-dist[i,j]^2 / (2 * sigma[sex2[i]]^2)) * z[i]
      y[i,j] ~ dpois(g[i,j] * K[j]) # K is different per trap so indexed
    } 
  }
  
  N <- sum(z[])
  D <- N/area
  #D_100 <- D * 100
}) 

Rmodel_m3 <- nimbleModel(m3, constants, data, inits_m3)

conf_m3 <- configureMCMC(Rmodel_m3, 
                         monitors = c("D", "N", "g0", "sigma", "psi", "pi", "logProb_y"),
                         enableWAIC = TRUE,
                         print = FALSE)
conf_m3$removeSamplers("s")
ACnodes_m3 <- paste0("s[", 1:M, ", 1:2]")
for(node in ACnodes_m3) {
  conf_m3$addSampler(target = node,
                     type = "RW_block",
                     control = list(adaptScaleOnly = TRUE),
                     silent = TRUE)
}
Rmcmc_m3 <- buildMCMC(conf_m3)

Cmodel_m3 <- compileNimble(Rmodel_m3)
Cmcmc_m3 <- compileNimble(Rmcmc_m3, project = Rmodel_m3)
MCMC_runtime_m3 <- system.time(
  samples_m3 <- runMCMC(Cmcmc_m3, 
                        nchains = chains, 
                        nburnin=burnin, 
                        niter = iter, 
                        WAIC = TRUE,
                        samplesAsCodaMCMC = TRUE)
)


mcmc.out_m3 <- coda::mcmc.list(samples_m3$samples)

log_prob_split_m3 <- split_logProb(mcmc.out_m3)

main_list_m3 <- log_prob_split_m3$non_logProb_mcmc_list
log_list_m3 <- log_prob_split_m3$logProb_mcmc_list

ll_mat_m3 <- as.matrix(log_list_m3)

loo_m3 <- loo::loo(ll_mat_m3, cores = 3) #takes awhile, about 2.5 min


output_m3 <- MCMCsummary(main_list_m3) %>% rownames_to_column("parameter")

#object for plotting
S_m3 <- ggs(main_list_m3) %>% 
  mutate(Parameter = fct_relevel(Parameter, "N", "D", "g0[1]", "g0[2]", "sigma[1]", "sigma[2]", "psi", "pi")) %>% 
  mutate(Parameter = fct_recode(Parameter, 
                                "g0[Male]" = "g0[1]",
                                "g0[Female]" = "g0[2]",
                                "sigma[Male]" = "sigma[1]", 
                                "sigma[Female]" = "sigma[2]"))

waic_tab_m3 <- tibble(model = "m3",
                      waic = samples_m3$WAIC$WAIC,
                      pwaic = samples_m3$WAIC$pWAIC,
                      lppd = samples_m3$WAIC$lppd,
                      elpd_loo = loo_m3$estimates["elpd_loo",1],
                      p_loo =  loo_m3$estimates["p_loo",1],
                      looic = loo_m3$estimates["looic",1])


#output the posterior summaries
#write_csv(output_m3, paste0("m3_70_day_baseline_ext_", year, "_results.csv"))


###############################################################################

#8. Model diagnostics

###############################################################################

mod_list <- list(S_m0, 
                 S_m1,
                 S_m2,
                 S_m3)

names(mod_list) <- c("tiger_scr_m0",
                     "tiger_scr_m1",
                     "tiger_scr_m2",
                     "tiger_scr_m3") 




for(i in 1:length(mod_list)){
  
  ################################ Convergence diagnostics (Rhat) #################################
  
  mod <- mod_list[[i]]
  mod_name <- names(mod_list)[i]
  
  
  #traceplots
  tps <- ggs_traceplot(mod)
  
  #ggsave(plot = tP_m3, filename = paste0(mod_name, "_70_day_baseline_ext_", year, "_traceplots.png"), width = 12, height = 8)
  
  ggs_diagnostics(mod)
  
  #note: simulated data may destabilize the N and D posteriors somewhat
  ggs_density(mod) 
  
  ggs_effective(mod)
  
  #distribution of density samples
  ggs_ppmean(S_m3, family="D", outcome = S_m3 %>% filter(Parameter == "D") %>% pull(value) %>% mean(na.rm=T))
  
  
  ################################ Autocorrelation plots #################################
  
  auto_corr_plot <- ggs_autocorrelation(mod)
  
  # ggsave(plot =  auto_corr_plot, paste0(mod_name, "_70_day_baseline_ext_", year, "_autocorrelation.png"), width = 10, height = 5, dpi = 300)
  
}
################################ Prior sensitivity analysis (m3)  #################################



#weakly informative prior for psi
m3_alt1 <- nimbleCode({ 
  
  # Priors 
  psi ~ dbeta(2,2)
  pi ~ dbeta(1,1) 
  
  g0[1] ~ dunif(0.01,1)
  g0[2] ~ dunif(0.01,1)
  sigma[1] ~ dunif(0.01,100)
  sigma[2] ~ dunif(0.01,100)
  
  for(i in 1:M){
    z[i] ~ dbern(psi)
    sex[i] ~ dbern(pi)
    sex2[i] <- sex[i] + 1
    
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    
    ## habitat constraint: 1 if the randomly generated AC falls in suitable habitat, -Inf if it does not.
    ones[i] ~ dHabitatMaskrb( s = s[i,1:2],
                              xmin = xlim[1],
                              xmax = xlim[2],
                              ymin = ylim[1],
                              ymax = ylim[2],
                              habitat = habitat.mx[1:y.max,1:x.max])
    
    for(j in 1:J){
      dist[i,j] <- sqrt((s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2)
      g[i,j] <- g0[sex2[i]] * exp(-dist[i,j]^2 / (2 * sigma[sex2[i]]^2)) * z[i]
      y[i,j] ~ dpois(g[i,j] * K[j]) # K is different per trap so indexed
    } 
  }
  
  N <- sum(z[])
  D <- N/area
  #D_100 <- D * 100
}) 

Rmodel_m3_alt1 <- nimbleModel(m3_alt1, constants, data, inits_m3)

conf_m3_alt1 <- configureMCMC(Rmodel_m3_alt1, 
                              monitors = c("D", "N", "g0", "sigma", "psi", "pi", "logProb_y"),
                              enableWAIC = TRUE,
                              print = FALSE)
conf_m3_alt1$removeSamplers("s")
ACnodes_m3_alt1 <- paste0("s[", 1:M, ", 1:2]")
for(node in ACnodes_m3_alt1) {
  conf_m3_alt1$addSampler(target = node,
                          type = "RW_block",
                          control = list(adaptScaleOnly = TRUE),
                          silent = TRUE)
}
Rmcmc_m3_alt1 <- buildMCMC(conf_m3_alt1)

Cmodel_m3_alt1 <- compileNimble(Rmodel_m3_alt1)
Cmcmc_m3_alt1 <- compileNimble(Rmcmc_m3_alt1, project = Rmodel_m3_alt1)
MCMC_runtime_m3_alt1 <- system.time(
  samples_m3_alt1 <- runMCMC(Cmcmc_m3_alt1, 
                             nchains = chains, 
                             nburnin=burnin, 
                             niter = iter, 
                             WAIC = TRUE,
                             samplesAsCodaMCMC = TRUE)
)


mcmc.out_m3_alt1 <- coda::mcmc.list(samples_m3_alt1$samples)

log_prob_split_m3_alt1 <- split_logProb(mcmc.out_m3_alt1)

main_list_m3_alt1 <- log_prob_split_m3_alt1$non_logProb_mcmc_list
log_list_m3_alt1 <- log_prob_split_m3_alt1$logProb_mcmc_list

ll_mat_m3_alt1 <- as.matrix(log_list_m3_alt1)

loo_m3_alt1 <- loo::loo(ll_mat_m3_alt1, cores = 3) #takes awhile, about 2.5 min


output_m3_alt1 <- MCMCsummary(main_list_m3_alt1) %>% rownames_to_column("parameter")

#object for plotting
S_m3_alt1 <- ggs(main_list_m3_alt1) %>% 
  mutate(Parameter = fct_relevel(Parameter, "N", "D", "g0[1]", "g0[2]", "sigma[1]", "sigma[2]", "psi", "pi")) %>% 
  mutate(Parameter = fct_recode(Parameter, 
                                "g0[Male]" = "g0[1]",
                                "g0[Female]" = "g0[2]",
                                "sigma[Male]" = "sigma[1]", 
                                "sigma[Female]" = "sigma[2]"))



#lognormal priors for sigma
m3_alt2 <- nimbleCode({ 
  
  # Priors 
  psi ~ dbeta(1,1)
  pi ~ dbeta(1,1) 
  
  g0[1] ~ dunif(0.01,1)
  g0[2] ~ dunif(0.01,1)
  sigma[1] ~ dlnorm(3, 1)
  sigma[2] ~ dlnorm(3, 1)
  
  for(i in 1:M){
    z[i] ~ dbern(psi)
    sex[i] ~ dbern(pi)
    sex2[i] <- sex[i] + 1
    
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    
    ## habitat constraint: 1 if the randomly generated AC falls in suitable habitat, -Inf if it does not.
    ones[i] ~ dHabitatMaskrb( s = s[i,1:2],
                              xmin = xlim[1],
                              xmax = xlim[2],
                              ymin = ylim[1],
                              ymax = ylim[2],
                              habitat = habitat.mx[1:y.max,1:x.max])
    
    for(j in 1:J){
      dist[i,j] <- sqrt((s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2)
      g[i,j] <- g0[sex2[i]] * exp(-dist[i,j]^2 / (2 * sigma[sex2[i]]^2)) * z[i]
      y[i,j] ~ dpois(g[i,j] * K[j]) # K is different per trap so indexed
    } 
  }
  
  N <- sum(z[])
  D <- N/area
  #D_100 <- D * 100
}) 

Rmodel_m3_alt2 <- nimbleModel(m3_alt2, constants, data, inits_m3)

conf_m3_alt2 <- configureMCMC(Rmodel_m3_alt2, 
                              monitors = c("D", "N", "g0", "sigma", "psi", "pi", "logProb_y"),
                              enableWAIC = TRUE,
                              print = FALSE)
conf_m3_alt2$removeSamplers("s")
ACnodes_m3_alt2 <- paste0("s[", 1:M, ", 1:2]")
for(node in ACnodes_m3_alt2) {
  conf_m3_alt2$addSampler(target = node,
                          type = "RW_block",
                          control = list(adaptScaleOnly = TRUE),
                          silent = TRUE)
}
Rmcmc_m3_alt2 <- buildMCMC(conf_m3_alt2)

Cmodel_m3_alt2 <- compileNimble(Rmodel_m3_alt2)
Cmcmc_m3_alt2 <- compileNimble(Rmcmc_m3_alt2, project = Rmodel_m3_alt2)
MCMC_runtime_m3_alt2 <- system.time(
  samples_m3_alt2 <- runMCMC(Cmcmc_m3_alt2, 
                             nchains = chains, 
                             nburnin=burnin, 
                             niter = iter, 
                             WAIC = TRUE,
                             samplesAsCodaMCMC = TRUE)
)


mcmc.out_m3_alt2 <- coda::mcmc.list(samples_m3_alt2$samples)

log_prob_split_m3_alt2 <- split_logProb(mcmc.out_m3_alt2)

main_list_m3_alt2 <- log_prob_split_m3_alt2$non_logProb_mcmc_list
log_list_m3_alt2 <- log_prob_split_m3_alt2$logProb_mcmc_list

ll_mat_m3_alt2 <- as.matrix(log_list_m3_alt2)

loo_m3_alt2 <- loo::loo(ll_mat_m3_alt2, cores = 3) #takes awhile, about 2.5 min


output_m3_alt2 <- MCMCsummary(main_list_m3_alt2) %>% rownames_to_column("parameter")

#object for plotting
S_m3_alt2 <- ggs(main_list_m3_alt2) %>% 
  mutate(Parameter = fct_relevel(Parameter, "N", "D", "g0[1]", "g0[2]", "sigma[1]", "sigma[2]", "psi", "pi")) %>% 
  mutate(Parameter = fct_recode(Parameter, 
                                "g0[Male]" = "g0[1]",
                                "g0[Female]" = "g0[2]",
                                "sigma[Male]" = "sigma[1]", 
                                "sigma[Female]" = "sigma[2]"))


samples_og <- S_m3 %>% 
  mutate(prior = "Uninformative", .before = everything())

samples_alt1 <- S_m3_alt1 %>% 
  mutate(prior = "Psi beta(2,2)", .before = everything())

samples_alt2 <- S_m3_alt2 %>% 
  mutate(prior = "Sigma lognormal(3,1)", .before = everything())



sens_dat <- list(samples_og ,
                 samples_alt1,
                 samples_alt2)


sens_dat <- bind_rows(sens_dat) %>% 
  rename(par = Parameter) %>% 
  mutate(prior = factor(prior),
         par = factor(par)) 



prior_sens_plot <- ggplot(sens_dat) +
  geom_density(aes(x = value, fill = prior)) +
  facet_wrap(vars(par), scales = "free") +
  theme(plot.background = element_rect(fill = "white"))

#ggsave(plot = prior_sens_plot, paste0("secr_m3_2021_prior_sens.png"), width = 10, height = 5, dpi = 300)





###############################################################################

#9. Create model selection table

###############################################################################
mod_weights <- loo::loo_model_weights(list(loo_m0, 
                                           loo_m1,
                                           loo_m2,
                                           loo_m3))

mod_weights_df <- mod_weights %>% 
  as.data.frame() %>% 
  rownames_to_column("model") %>% 
  mutate(model =fct_recode(model, 
                           "m0" = "model1",
                           "m1" = "model2",
                           "m2" = "model3",
                           "m3" = "model4")) %>% 
  rename(weight = x)

waic_tab <- bind_rows(waic_tab_m0,
                      waic_tab_m1,
                      waic_tab_m2,
                      waic_tab_m3) %>% 
  arrange(looic) %>% 
  mutate(dwaic = waic - waic[1],
         dlooic = looic - looic[1],
         pmps1 = exp(-0.5 * (waic - min(waic)))) %>% 
  mutate(pmps2 = sum(pmps1)/pmps1) %>% 
  rename(pmps = pmps2) %>% 
  left_join(mod_weights_df, by = "model") %>% 
  select(model, looic, dlooic, waic, dwaic, weight, pwaic, elpd_loo, p_loo, lppd, pmps)

#write_csv(waic_tab, paste0("mod_sel_70_day_baseline_ext_", year, ".csv"))



###############################################################################

#10. Compute model-averaged parameters

###############################################################################

sample_list <- list(main_list_m0,
                    main_list_m1,
                    main_list_m2,
                    main_list_m3)

names(sample_list) <- c("m0", "m1", "m2", "m3")

sample_list_df <- map(sample_list, 
                      function(x){x %>% 
                          map(as.data.frame) %>% 
                          imap( ~ mutate(.x, source = .y, .before = D)) %>%
                          bind_rows()})


#export samples for later use

for(i in 1:length(sample_list_df)){
  name <- names(sample_list_df)[i]
  write_csv(sample_list_df[[i]], paste0("tiger_nimble_secr_samples_baseline_ext_", name, "_", year, ".csv"))
}



cov_names <- c("N", "D", "g0", "sigma","g0[1]", "g0[2]", "sigma[1]", "sigma[2]", "psi", "pi")

mod_names <-names(sample_list_df)

avg_covs <- list()

for(i in 1:length(cov_names)){
  cov_i <- cov_names[i]
  
  cov_mods <- select_dataframes(sample_list_df, cov_i)
  
  cov_mod_names <- names(cov_mods)
  
  weight_vals <- mod_weights_df %>%
    filter(model %in% cov_mod_names) %>%
    pull(weight) %>% as.numeric()
  
  weighted_samples_i <- list()
  for(j in 1:length(cov_mods)){
    weighted_samples_i[[j]] <- cov_mods[[j]][,cov_i]*weight_vals[j]
    #print(j)
  }
  
  avg_covs[[i]] <-  Reduce('+', weighted_samples_i)
  
  print(i)
}

names(avg_covs) <- cov_names

avg_covs_df <- as.data.frame(avg_covs)

avg_covs_mcmc <- coda::as.mcmc(avg_covs_df)

mcmc_list <- coda::mcmc.list(avg_covs_mcmc)

summary <- MCMCvis::MCMCsummary(mcmc_list) %>% rownames_to_column("parameter")

#write_csv(summary, paste0("tiger_nimble_secr_mod_avg_coeffs_baseline_ext_", year, ".csv"))





