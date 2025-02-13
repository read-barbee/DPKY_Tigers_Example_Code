#### DPKY Prey Occupancy Models ####

# Author: Read Barbee

# Date:2025-02-13 

# Purpose: Fit occupancy models to prey data and generate model-averaged coefficients


################################ Libraries #################################
library(tidyverse)
library(ubms)
library(sf)
library(mapview)
library(ggspatial)
library(terra)
library(bayesplot)


################################ Run this code to use the old version of the ubms fitList function that returns model weights #################################
setClass("ubmsFitList", slots = c(models = "list"))

setGeneric("fitList", function(...){
  unmarked::fitList(...)
})

#' Create a List of ubmsFit Models
#'
#' Create a list of ubmsFit models
#'
#' @param ... \code{ubmsFit} model objects, preferably named, or a list
#'  of such models
#'
#' @return An object of class \code{ubmsFitList} containing the list of models
#'
#' @aliases fitList fitList,list-method
#' @export
setMethod("fitList", "ubmsFit", function(...){
  mods <- list(...)
  mod_names <- names(mods)
  if(is.null(mod_names)) mod_names <- rep("", length(mods))
  mod_names[mod_names==""] <- NA
  obj_names=sapply(substitute(list(...))[-1], deparse)
  mod_names[is.na(mod_names)] <- obj_names[is.na(mod_names)]
  names(mods) <- mod_names
  fitList(mods)
})

setMethod("fitList", "list", function(...){
  mods <- list(...)[[1]]
  if(!inherits(mods[[1]],"ubmsFit", )) return(unmarked::fitList(fits=mods))
  if(is.null(names(mods))) names(mods) <- paste0("mod", 1:length(mods))
  new("ubmsFitList", models=mods)
})

#' Model Selection For a List of ubmsFit Models
#'
#' Construct a model selection table from a \code{ubmsFitList}
#'
#' @param object An object of class \code{ubmsFitList}
#' @param ... Currently ignored
#'
#' @return A \code{data.frame} of model fit information with one row per
#'  model in the input \code{ubmsFitList}. Models are ranked in descending
#'  order by expected log pointwise predictive density (\code{elpd}).
#' @seealso \code{\link[loo]{loo}}, \code{\link[loo]{loo_compare}}
#'
#' @importFrom loo loo_compare loo_model_weights
#' @importFrom unmarked modSel
#' @export
setMethod("modSel", "ubmsFitList", function(object, ...){
  #loos <- lapply(object@models, loo, ...)
  loos <- lapply(object@models, function(x) x@loo)
  elpd <- sapply(loos, function(x) x$estimates[1])
  p_loo <- sapply(loos, function(x) x$estimates[2])
  compare <- loo::loo_compare(loos)[names(elpd),]
  wts <- as.vector(loo::loo_model_weights(loos))
  out <- data.frame(elpd=elpd, nparam=p_loo, elpd_diff=compare[,1],
                    se_diff=compare[,2], weight=wts)
  out[order(out$elpd_diff, decreasing=TRUE),]
})




################################ Helper functions #################################

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


################################ User-defined parameters #################################
chains <- 3
iter = 10000
cores = chains
###############################################################################

#1. Import data

###############################################################################


#camera activity matrix
cam_act <- read_csv("data/occupancy/camera_activity/cam_act_daily_2021.csv")


#prey covariates
covs <- read_csv("data/occupancy/covariates/prey_covariates_demo.csv") 



#convert cam activity df to matrix
cam_act_mat <- cam_act %>% as.matrix()


effort_mat_binary <- cam_act_mat
effort_mat_binary[is.na(effort_mat_binary)] <- 0


#Species detection histories
gaur_demo <- read_csv("data/occupancy/detection_histories_demo/prey/gaur_det_hist_demo.csv")

banteng_demo <- read_csv("data/occupancy/detection_histories_demo/prey/banteng_det_hist_demo.csv")

serow_demo <- read_csv("data/occupancy/detection_histories_demo/prey/serow_det_hist_demo.csv")

muntjac_demo <- read_csv("data/occupancy/detection_histories_demo/prey/muntjac_det_hist_demo.csv")

sambar_demo <- read_csv("data/occupancy/detection_histories_demo/prey/sambar_det_hist_demo.csv")

wild_pig_demo <- read_csv("data/occupancy/detection_histories_demo/prey/wild_pig_det_hist_demo.csv")

mouse_deer_demo <- read_csv("data/occupancy/detection_histories_demo/prey/mouse_deer_det_hist_demo.csv")

det_hists <- list(gaur_demo,
                  banteng_demo,
                  serow_demo,
                  muntjac_demo,
                  sambar_demo,
                  wild_pig_demo,
                  mouse_deer_demo)

names(det_hists) <- c("gaur",
                      "banteng",
                      "serow",
                      "muntjac",
                      "sambar",
                      "wild_pig",
                      "mouse_deer")



###############################################################################

#2. Make unmarked frames

###############################################################################

#Note: the daily activity matrix is added to the dataframe for completeness here, but a precalculated site-level count of trap days is used as the detection covariate. 

umfs <- list()

for (i in 1:length(det_hists)){
  species_i <- names(det_hists)[i]
  y_i <- det_hists[[species_i]]
  
  umfs[[i]] <- unmarkedFrameOccu(y = y_i,
                         siteCovs = covs,
                         obsCovs = list(effort = effort_mat_binary))
}

names(umfs) <- names(det_hists)




###############################################################################

#3. Fit candidate models for each species

###############################################################################

#choose which prey species you want to run the analysis for by setting the i value
species_i <-names(det_hists)[i]
  
  umf <- umfs[[species_i]]
  
  #canidate models 
  all_covs <- stan_occu(~site_effort ~elevation + ndvi + dist_water + dist_roads + dist_fe, umf, chains = chains, iter = iter, cores = cores)
  
  environment <- stan_occu(~site_effort ~elevation + ndvi + dist_water, umf, chains = chains, iter = iter, cores = cores)
  
  disturbance <- stan_occu(~site_effort ~dist_roads + dist_fe, umf, chains = chains, iter = iter, cores = cores)
  
  block <- stan_occu(~site_effort ~block, umf, chains = chains, iter = iter, cores = cores)
  
  all_covs_block <- stan_occu(~site_effort ~elevation + ndvi + dist_water + dist_roads + dist_fe + block + (elevation*block) + (ndvi*block) + (dist_water*block) + (dist_roads*block) + (dist_fe*block), umf, chains = chains, iter = iter, cores = cores)
  
  environment_block <- stan_occu(~site_effort ~elevation + ndvi + dist_water + block + (elevation*block) + (ndvi*block) + (dist_water*block), umf, chains = chains, iter = iter, cores = cores)
  
  disturbance_block <- stan_occu(~site_effort ~dist_roads + dist_fe + block + (dist_roads*block) + (dist_fe*block), umf, chains = chains, iter = iter, cores = cores)
  
  mod_list <- list(all_covs, 
                   environment,
                   disturbance,
                   block,
                   all_covs_block,
                   environment_block,
                   disturbance_block)
  
  names(mod_list) <- c(paste0(species_i, "_all_covs"), 
                       paste0(species_i, "_environment"),
                       paste0(species_i, "_disturbance"),
                       paste0(species_i, "_block"),
                       paste0(species_i, "_all_covs_block"),
                       paste0(species_i, "_environment_block"),
                       paste0(species_i, "_disturbance_block")) %>% 
    str_replace_all(" ", "_")
  
###############################################################################
  
#4. Calculate model diagnostics
  
###############################################################################
  
  
  for(i in 1:length(mod_list)){
    
    ################################ Convergence diagnostics (Rhat) #################################
    
    mod <- mod_list[[i]]
    mod_name <- names(mod_list)[i]
    
    mod_summary <- summary(mod, submodel = "state") %>% 
      rownames_to_column("parameter")
    
    
    #write_csv(mod_summary, paste0("mod_name, "_summary.csv"))
    
    
    #traceplots
    tps <- traceplot(mod, pars = c("beta_state", "beta_det"))
    
    #ggsave(plot = tps, paste0(mod_name, "_traceplots.png"), width = 10, height = 5, dpi = 300)
    
    ################################ Posterior predictive checks #################################
    
    #ubms ppc
    ubms_ppc <- gof(mod, draws=1000, quiet=TRUE)
    
    ubms_ppc_df <- enframe(c(ubms_ppc@statistic, ubms_ppc@estimate, ubms_ppc@post_pred_p))
    
    ubms_ppc_df <- tibble(mackenzie_bailey_chi_square_estimate = ubms_ppc@estimate,
                          post_pred_p = ubms_ppc@post_pred_p)
    
    #write_csv(ubms_ppc_df, paste0(mod_name, "_mb_chi_square.csv"))
    
    #calculate observed success proportion for each site across all sampling occasions
    y_obs <- mod@data@y %>% 
      as.data.frame() %>% 
      rowwise() %>% 
      summarize(proportion = mean(c_across(everything()), na.rm=T)) %>% 
      ungroup() %>% 
      pull(proportion)
    
    
    # Extract posterior predictive samples
    y_rep <- posterior_predict(mod, "y", draws=1000)
    
    n_draws <- 1000
    n_sites <- 187
    n_obs_per_site <- 243
    
    # Reshape the matrix to 3D: [n_draws, n_obs_per_site, n_sites]
    y_array <- array(y_rep, dim = c(n_draws, n_obs_per_site, n_sites))
    
    # Calculate row-wise means for each site
    post_preds <- apply(y_array, c(1, 3), mean, na.rm=T)  # Result: [n_draws, n_sites]
    
    
    #plot posterior predictive checks
    
    ppc_dens <- pp_check(y_obs, yrep = post_preds, fun = "dens_overlay") +
      theme(plot.background = element_rect(fill = "white"))
    
    ppc_ecdf <- pp_check(y_obs, yrep = post_preds, fun = "ecdf_overlay") +
      theme(plot.background = element_rect(fill = "white"))
    
    ppc_hist <- pp_check(y_obs, yrep = post_preds, fun = "stat") +
      theme(plot.background = element_rect(fill = "white"))
    
    
    
    #ggsave(plot = ppc_dens, paste0(mod_name, "_ppc_dens.png"), width = 7, height = 5, dpi = 300)
    
    #ggsave(plot = ppc_hist, paste0(mod_name, "_ppc_hist.png"), width = 7, height = 5, dpi = 300)
    
    
    ################################ Autocorrelation plots #################################
    
    mod_stanfit <-mod@stanfit
    
    #mcmc_rhat_data(rhats)
    
    acf_state <- mcmc_acf(mod_stanfit, regex_pars = "beta_state") +
      theme(plot.background = element_rect(fill = "white"))
    
    acf_det <- mcmc_acf(mod_stanfit, regex_pars = "beta_det") +
      theme(plot.background = element_rect(fill = "white"))
    
    
   # ggsave(plot = acf_state, paste0(mod_name, "_autocorrelation_state.png"), width = 25, height = 5, dpi = 300)
    
    #ggsave(plot = acf_det, paste0(mod_name, "_autocorrelation_det.png"), width = 7, height = 5, dpi = 300)
    
    ################################ Prior sensitivity analysis  #################################
    
    mod_form <- eval(mod@call$formula)
    mod_dat <- eval(mod@call$data)
    
    state_pars <- mod@submodels@submodels$state@formula %>% as.character()
    
    state_pars <- state_pars[2] %>% 
      str_split("\\+") %>% 
      unlist() %>% 
      str_squish()
    
    state_pars <- c("intercept_state", state_pars)
    
    det_pars = c("intercept_det", "site_effort")
    

#modify priors for inf_fit 1 and inf_fit2 as desired
    
    inf_fit1 <- stan_occu(mod_form, 
                          mod_dat,
                          prior_intercept_state = logistic(0, 0.5),
                          prior_coef_state = logistic(0, 0.5),
                          prior_intercept_det = logistic(0, 0.5),
                          prior_coef_det = logistic(0, 0.5),
                          chains = chains, 
                          iter = iter, 
                          cores = cores)
    
    inf_fit2 <- stan_occu(mod_form, 
                          mod_dat,
                          prior_intercept_state = logistic(0, 0.25),
                          prior_coef_state = logistic(0, 0.25),
                          prior_intercept_det = logistic(0, 0.25),
                          prior_coef_det = logistic(0, 0.25),
                          chains = chains, 
                          iter = iter, 
                          cores = cores)
    
    
    post_samples_orig_state <- ubms::extract(mod, pars = "beta_state" ,permuted = T) %>% 
      as.data.frame() 
    
    names(post_samples_orig_state) <- state_pars
    
    post_samples_orig_state <- post_samples_orig_state %>% 
      mutate(prior = "logistic(0,1)", .before = everything())
    
    
    post_samples_orig_det <- ubms::extract(mod, pars = "beta_det" ,permuted = T) %>% 
      as.data.frame() 
    
    names(post_samples_orig_det) <- det_pars
    
    post_samples_orig_det <- post_samples_orig_det %>% 
      mutate(prior = "logistic(0,1)", .before = everything())
    
    
    post_samples_inf1_state <- ubms::extract(inf_fit1, pars = "beta_state" ,permuted = T) %>% 
      as.data.frame() 
    
    names(post_samples_inf1_state) <- state_pars
    
    post_samples_inf1_state<- post_samples_inf1_state %>% 
      mutate(prior = "logistic(0,0.5)", .before = everything())
    
    
    post_samples_inf1_det <- ubms::extract(inf_fit1, pars = "beta_det" ,permuted = T) %>% 
      as.data.frame() 
    
    names(post_samples_inf1_det) <- det_pars
    
    post_samples_inf1_det <- post_samples_inf1_det %>% 
      mutate(prior = "logistic(0,0.5)", .before = everything())
    
    
    
    post_samples_inf2_state <- ubms::extract(inf_fit2, pars = "beta_state" ,permuted = T) %>% 
      as.data.frame() 
    
    names(post_samples_inf2_state) <- state_pars
    
    post_samples_inf2_state<- post_samples_inf2_state %>% 
      mutate(prior = "logistic(0,0.25)", .before = everything())
    
    
    post_samples_inf2_det <- ubms::extract(inf_fit2, pars = "beta_det" ,permuted = T) %>% 
      as.data.frame() 
    
    names(post_samples_inf2_det) <- det_pars
    
    post_samples_inf2_det <- post_samples_inf2_det %>% 
      mutate(prior = "logistic(0,0.25)", .before = everything())
    
    
    sens_dat <- list(post_samples_orig_state,
                     post_samples_orig_det,
                     post_samples_inf1_state,
                     post_samples_inf1_det,
                     post_samples_inf2_state,
                     post_samples_inf2_det)
    
    sens_dat <- map(sens_dat, function(x){x %>% pivot_longer(-prior, names_to = "par", values_to = "value")})
    
    sens_dat <- bind_rows(sens_dat) %>% 
      mutate(prior = factor(prior),
             par = factor(par)) %>% 
      mutate(prior = fct_recode(prior, "logistic(0,1) [Default]" = "logistic(0,1)")) %>% 
      mutate(par = fct_relevel(par, "intercept_state", after = 0)) %>%        
      mutate(par = fct_relevel(par, "intercept_det", "site_effort", after = Inf)) %>% 
      mutate(par = fct_relabel(par, function(x){x %>% 
          str_to_sentence() %>% 
          str_replace("_", " ")}))
    
    if(any(str_detect(levels(sens_dat$par), "Ndvi"))){
      sens_dat <- sens_dat %>% 
        mutate(par = fct_recode(par, "NDVI" = "Ndvi"))
    }
    
    prior_sens_plot <- ggplot(sens_dat) +
      geom_density(aes(x = value, fill = prior)) +
      facet_wrap(vars(par), scales = "free") +
      theme(plot.background = element_rect(fill = "white"))
    
    #ggsave(plot = prior_sens_plot, paste0(mod_name, "_prior_sens.png"), width = 10, height = 5, dpi = 300)
    
    
    
    print(paste0(i, "/", length(mod_list)))
  }
  
  
  
  
  
###############################################################################
  
#5. Model selection
  
###############################################################################
  
  fit_list <- fitList(mod_list)
  
  sel <- modSel(fit_list)
  
  
  elpd_loo <- vector()
  p_loo <- vector()
  looic <- vector()
  elpd_waic <- vector()
  p_waic <- vector()
  waic <- vector()
  
  for(i in 1:length(mod_list)){
    loo <- ubms::loo(mod_list[[i]])
    waic_l <- ubms::waic(mod_list[[i]])
    elpd_loo[i] <- loo$estimates["elpd_loo",1]
    p_loo[i] <- loo$estimates["p_loo",1]
    looic[i] <- loo$estimates["looic",1]
    elpd_waic[i] <- waic_l$estimates["elpd_waic",1]
    p_waic[i] <- waic_l$estimates["p_waic",1]
    waic[i] <- waic_l$estimates["waic",1]
    
    print(i)
  }
  
  sel2 <- tibble( mod = names(mod_list),
                  elpd_loo = elpd_loo,
                  p_loo = p_loo,
                  looic = looic,
                  elpd_waic = elpd_waic,
                  p_waic = p_waic,
                  waic = waic
                  )
  
  sel_final <- sel %>% 
    rownames_to_column("mod") %>% 
    as_tibble() %>% 
    left_join(sel2, by = "mod") %>% 
    mutate(dlooic = looic - min(looic),
           dwaic = waic - min(waic))
    
#}

species_name <- species_i %>% str_replace(coll(" "), coll("_"))

#write_csv(sel_final, paste0(species_name, "_mod_sel_ubms_site_eff.csv"))


###############################################################################

#6. Coefficient estimates

###############################################################################

for(i in 1:length(mod_list)){
  
mod_name <- names(mod_list)[i] 

state_summ <- summary(mod_list[[i]], submodel = "state") %>% 
  rownames_to_column("parameter") %>% 
  mutate(parameter = str_replace(parameter, coll("(Intercept)"), coll("psi_int")))

det_summ <- summary(mod_list[[i]], submodel = "det") %>% 
  rownames_to_column("parameter") %>% 
  mutate(parameter = str_replace(parameter, coll("(Intercept)"), coll("p_int")))

full_summ <- bind_rows(state_summ, det_summ)

write_csv(full_summ, paste0(species_name, "_",  mod_name, "_coef_summary_site_eff.csv"))

print(i)

}

  

###############################################################################

#7. Model average coefficients

###############################################################################
  weights <- tibble(mod = names(mod_list)) %>%
    left_join(sel_final %>% select(mod, weight), by = "mod")

  #extract marginal posterior samples

  #all covs
  all_covs_state <- ubms::extract(all_covs)$beta_state %>% as.data.frame()
  colnames(all_covs_state) <-  c("psi_int","elevation", "ndvi", "dist_water", "dist_roads", "dist_fe")
  all_covs_det <- ubms::extract(all_covs)$beta_det %>% as.data.frame()
  colnames(all_covs_det) <-  c("p_int", "effort")

  #environment
  environment_state <- ubms::extract(environment)$beta_state %>% as.data.frame()
  colnames(environment_state) <-  c("psi_int","elevation", "ndvi", "dist_water")
 environment_det <- ubms::extract(environment)$beta_det %>% as.data.frame()
  colnames(environment_det) <-  c("p_int", "effort")

  #disturbance
  disturbance_state <- ubms::extract(disturbance)$beta_state %>% as.data.frame()
  colnames(disturbance_state) <-  c("psi_int", "dist_roads", "dist_fe")
  disturbance_det <- ubms::extract(disturbance)$beta_det %>% as.data.frame()
  colnames(disturbance_det) <-  c("p_int", "effort")

  #block
  block_state <- ubms::extract(block)$beta_state %>% as.data.frame()
  colnames(block_state) <-  c("psi_int", "block")
  block_det <- ubms::extract(block)$beta_det %>% as.data.frame()
  colnames(block_det) <-  c("p_int", "effort")

  #all covs block
  all_covs_block_state <- ubms::extract(all_covs_block)$beta_state %>% as.data.frame()
  colnames(all_covs_block_state) <-  c("psi_int","elevation", "ndvi", "dist_water", "dist_roads", "dist_fe", "block","elevation:block", "ndvi:block", "dist_water:block", "dist_roads:block", "dist_fe:block")
  all_covs_block_det <- ubms::extract(all_covs_block)$beta_det %>% as.data.frame()
  colnames(all_covs_block_det) <-  c("p_int", "effort")

  #environment block
  environment_block_state <- ubms::extract(environment_block)$beta_state %>% as.data.frame()
  colnames(environment_block_state) <-  c("psi_int","elevation", "ndvi", "dist_water", "block", "elevation:block", "ndvi:block", "dist_water:block")
  environment_block_det <- ubms::extract(environment_block)$beta_det %>% as.data.frame()
  colnames(environment_block_det) <-  c("p_int", "effort")

  #disturbance block
  disturbance_block_state <- ubms::extract(disturbance_block)$beta_state %>% as.data.frame()
  colnames(disturbance_block_state) <-  c("psi_int", "dist_roads", "dist_fe", "block", "dist_roads:block", "dist_fe:block")
  disturbance_block_det <- ubms::extract(disturbance_block)$beta_det %>% as.data.frame()
  colnames(disturbance_block_det) <-  c("p_int", "effort")

  state_list <- list(all_covs_state,
                     environment_state,
                     disturbance_state,
                     block_state,
                     all_covs_block_state,
                     environment_block_state,
                     disturbance_block_state)

  names(state_list) <- list("all_covs_state",
                            "environment_state",
                            "disturbance_state",
                            "block_state",
                            "all_covs_block_state",
                            "environment_block_state",
                            "disturbance_block_state")


  cov_names_state <- colnames(all_covs_block_state)
  mod_names <-names(mod_list)

  avg_covs <- list()
  

  for(i in 1:length(cov_names_state)){
    cov_i <- cov_names_state[i]

    cov_mods <- select_dataframes(state_list, cov_i)

    cov_mod_names <- names(cov_mods) %>% str_remove(coll("_state"))
    
    cov_mod_names <- paste0(species_i, "_", cov_mod_names)

    weight_vals <- weights %>%
      filter(mod %in% cov_mod_names) %>%
      pull(weight)

    weighted_samples_i <- list()
    for(j in 1:length(cov_mods)){
    weighted_samples_i[[j]] <- cov_mods[[j]][,cov_i]*weight_vals[j]
    #print(j)
    }

    avg_covs[[i]] <-  Reduce('+', weighted_samples_i) 

    print(i)
  }

names(avg_covs) <- cov_names_state

avg_covs_df <- as.data.frame(avg_covs)

avg_covs_mcmc <- coda::as.mcmc(avg_covs_df)

mcmc_list <- coda::mcmc.list(avg_covs_mcmc)

summary <- MCMCvis::MCMCsummary(mcmc_list) %>% rownames_to_column("parameter")

#write_csv(summary, paste0(species_name, "_mod_avg_coeffs_site_eff.csv"))


# Detection

det_list <- list(all_covs_det,
                 environment_det,
                 disturbance_det,
                 block_det,
                 all_covs_block_det,
                 environment_block_det,
                 disturbance_block_det)

names(det_list) <- list("all_covs_det",
                        "environment_det",
                        "disturbance_det",
                        "block_det",
                        "all_covs_block_det",
                        "environment_block_det",
                        "disturbance_block_det")


cov_names_det <- colnames(all_covs_block_det)
avg_covs_det <- list()

for(i in 1:length(cov_names_det)){
  cov_i <- cov_names_det[i]
  
  cov_mods <- select_dataframes(det_list, cov_i)
  
  cov_mod_names <- names(cov_mods) %>% str_remove(coll("_det"))
  cov_mod_names <- paste0(species_i, "_", cov_mod_names)
  
  weight_vals <- weights %>%
    filter(mod %in% cov_mod_names) %>%
    pull(weight)
  
  weighted_samples_i <- list()
  for(j in 1:length(cov_mods)){
    weighted_samples_i[[j]] <- cov_mods[[j]][,cov_i]*weight_vals[j]
    #print(j)
  }
  
  avg_covs_det[[i]] <-  Reduce('+', weighted_samples_i) 
  
  print(i)
}

names(avg_covs_det) <- cov_names_det

avg_covs_det_df <- as.data.frame(avg_covs_det)

avg_covs_det_mcmc <- coda::as.mcmc(avg_covs_det_df)

mcmc_list_det <- coda::mcmc.list(avg_covs_det_mcmc)

summary_det <- MCMCvis::MCMCsummary(mcmc_list_det) %>% rownames_to_column("parameter")

#write_csv(summary_det, paste0(species_name, "_mod_avg_coeffs_site_eff_p.csv"))



###############################################################################

#8.Model average plots

###############################################################################

#object for plotting
plot_obj <- ggmcmc::ggs(mcmc_list)


plot_obj2 <- plot_obj %>% 
  filter(Parameter %in% c("psi_int",
                          "elevation",
                          "ndvi",
                          "dist_water",
                          "dist_roads",
                          "dist_fe")) %>% 
  mutate(Parameter =factor(Parameter, levels = c("psi_int", "elevation", "ndvi", "dist_water", "dist_roads", "dist_fe")))

avg_dists_plot <- ggmcmc::ggs_density(plot_obj2)
#ggsave(avg_dists_plot, file=paste0(species_name, "_mod_avg_coeffs_site_eff.png"), width = 12, height = 8)

