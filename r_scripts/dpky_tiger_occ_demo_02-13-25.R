#### DPKY Tiger Occupancy Models ####

# Author: Read Barbee

# Date:2025-02-10 

# Purpose: Fit occupancy models to tiger data and generate model-averaged coefficients

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

#tiger detection history (simulated)
det_hist <- read_csv("data/occupancy/detection_histories_demo/tiger/tiger_det_hist_demo.csv")

#camera activity matrix
cam_act <- read_csv("data/occupancy/camera_activity/cam_act_daily_2021.csv")


#tiger covariates
covs <- read_csv("data/occupancy/covariates/tiger_covariates_demo.csv") 


#convert cam activity df to matrix
cam_act_mat <- cam_act %>% as.matrix()


effort_mat_binary <- cam_act_mat
effort_mat_binary[is.na(effort_mat_binary)] <- 0




###############################################################################

#2. Make unmarked frame

###############################################################################

#Note: the daily activity matrix is added to the dataframe for completeness here, but a precalculated site-level count of trap days is used as the detection covariate. 

umf <-  unmarkedFrameOccu(y = det_hist,
                          siteCovs = covs,
                          obsCovs = list(effort = effort_mat_binary))


###############################################################################

#3. Fit candidate models

###############################################################################

#non-prey models
all_covs <- stan_occu(~site_effort ~elevation + ndvi + dist_water + dist_roads + dist_fe, umf, chains = chains, iter = iter, cores = cores)
environment <- stan_occu(~site_effort ~elevation + ndvi + dist_water, umf, chains = chains, iter = iter, cores = cores)
disturbance <- stan_occu(~site_effort ~dist_roads + dist_fe, umf, chains = chains, iter = iter, cores = cores)
block <- stan_occu(~site_effort ~block, umf, chains = chains, iter = iter, cores = cores)
all_covs_block <- stan_occu(~site_effort ~elevation + ndvi + dist_water + dist_roads + dist_fe + block + (elevation*block) + (ndvi*block) + (dist_water*block) + (dist_roads*block) + (dist_fe*block), umf, chains = chains, iter = iter, cores = cores)
environment_block <- stan_occu(~site_effort ~elevation + ndvi + dist_water + block + (elevation*block) + (ndvi*block) + (dist_water*block), umf, chains = chains, iter = iter, cores = cores)
disturbance_block <- stan_occu(~site_effort ~dist_roads + dist_fe + block + (dist_roads*block) + (dist_fe*block), umf, chains = chains, iter = iter, cores = cores)

#prey models
all_prey <- stan_occu(~site_effort ~ gaur + muntjac + sambar + wild_pig + mouse_deer , umf, chains = chains, iter = iter, cores = cores) 
all_prey_block <- stan_occu(~site_effort ~gaur + muntjac + sambar + wild_pig + mouse_deer + block + (gaur*block) + (muntjac*block) + (sambar*block) + (wild_pig*block) + (mouse_deer*block), umf, chains = chains, iter = iter, cores = cores)
main_prey <- stan_occu(~site_effort ~ sambar + wild_pig, umf, chains = chains, iter = iter, cores = cores)
main_prey_block <- stan_occu(~site_effort ~ sambar + wild_pig + block + (sambar*block) + (wild_pig*block), umf, chains = chains, iter = iter, cores = cores)
large_prey <- stan_occu(~site_effort ~ gaur + sambar, umf, chains = chains, iter = iter, cores = cores) 
large_prey_block <- stan_occu(~site_effort ~ gaur + sambar  + block + (gaur*block) + (sambar*block) , umf, chains = chains, iter = iter, cores = cores)


mod_list <- list(all_covs, 
                 environment,
                 disturbance,
                 block,
                 all_covs_block,
                 environment_block,
                 disturbance_block,
                 all_prey,
                 all_prey_block,
                 main_prey,
                 main_prey_block,
                 large_prey,
                 large_prey_block
                 )

names(mod_list) <- c("tiger_all_covs", 
                     "tiger_environment",
                     "tiger_disturbance",
                     "tiger_block",
                     "tiger_all_covs_block",
                     "tiger_environment_block",
                     "tiger_disturbance_block",
                     "tiger_all_prey",
                     "tiger_all_prey_block",
                     "tiger_main_prey",
                     "tiger_main_prey_block",
                     "tiger_large_prey",
                     "tiger_large_prey_block")




###############################################################################

#4. Calculate model diagnostics

###############################################################################



for(i in 1:length(mod_list)){

################################ Convergence diagnostics (Rhat) #################################
  
mod <- mod_list[[i]]
mod_name <- names(mod_list)[i]

mod_summary <- summary(mod, submodel = "state") %>% 
  rownames_to_column("parameter")


#write_csv(mod_summary, paste0(mod_name, "_summary.csv"))


#traceplots
tps <- traceplot(mod, pars = c("beta_state", "beta_det"))

#ggsave(plot = tps, paste0(mod_name, "_traceplots.png"), width = 7, height = 5, dpi = 300)

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

#pp_check(y_obs, yrep = post_preds, fun = "bars")

#ggsave(plot = ppc_dens, paste0(mod_name, "_ppc_dens.png"), width = 7, height = 5, dpi = 300)

#ggsave(plot = ppc_hist, paste0(mod_name, "_ppc_hist.png"), width = 7, height = 5, dpi = 300)


################################ Autocorrelation plots #################################

mod_stanfit <-mod@stanfit


acf_state <- mcmc_acf(mod_stanfit, regex_pars = "beta_state") +
  theme(plot.background = element_rect(fill = "white"))

acf_det <- mcmc_acf(mod_stanfit, regex_pars = "beta_det") +
  theme(plot.background = element_rect(fill = "white"))



#ggsave(plot = acf_state, paste0(mod_name, "_autocorrelation_state.png"), width = 12, height = 5, dpi = 300)

#ggsave(plot = acf_det, paste0(mod_name, "_autocorrelation_det.png"), width = 7, height = 5, dpi = 300)

################################ Prior sensitivity analysis  #################################

mod_form <- eval(mod@call$formula)
mod_dat <- eval(mod@call$data)

state_pars <- mod@submodels@submodels$state@formula %>% as.character()

state_pars <- state_pars[2] %>% 
  str_split("\\+") %>% 
  unlist() %>% 
  str_squish()

state_pars <- c("intercept_state", state_pars) #%>% str_to_sentence()

det_pars = c("intercept_det", "site_effort")


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

#ggsave(plot = prior_sens_plot, paste0(mod_name, "_prior_sens.png"), width = 7, height = 5, dpi = 300)



print(paste0(i, "/", length(mod_list)))
}
###############################################################################

#5. Model selection

###############################################################################


fit_list <- fitList(mod_list)

sel <- modSel(fit_list, simplify)


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

#write_csv(sel_final, paste0("tiger_occ_mod_sel_ubms_no_bs_site_eff2.csv"))


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
  
  #write_csv(full_summ, paste0("tiger_",  mod_name, "_coef_summary_no_bs.csv"))
  
  print(i)
  
}


###############################################################################

#7. Model average coefficients

###############################################################################
weights <- tibble(mod = names(mod_list)) %>%
  left_join(sel_final %>% select(mod, weight), by = "mod")

#extract marginal posterior samples

state_list <- list()
det_list <- list()

for(i in 1:length(mod_list)){
  
  mod = mod_list[[i]]
  mod_name <- names(mod_list)[i]
  mod_name_state <- paste0(mod_name, "_state")
  mod_name_det <- paste0(mod_name, "_det")

par_names <- as.character(mod@submodels@submodels$state@formula)[2] %>% 
  str_split(" \\+ ", simplify = TRUE) %>% as.vector()

par_names_all <- c("psi_int", par_names)

state_list[[mod_name_state]] <-  ubms::extract(mod)$beta_state %>% 
         as.data.frame()
  
colnames(state_list[[mod_name_state]]) <- par_names_all

det_list[[mod_name_det]] <- ubms::extract(mod)$beta_det %>% as.data.frame()
colnames(det_list[[mod_name_det]]) <-  c("p_int", "effort")

print(i)

}

### Calculate model averaged coefficients ###

cov_names_state <- c("psi_int", names(covs %>% select(-c(cell_id, site_effort, roads_hii, serow, banteng)))) 


mod_names <-names(mod_list)

avg_covs <- list()

for(i in 1:length(cov_names_state)){
  cov_i <- cov_names_state[i]
  
  cov_mods <- select_dataframes(state_list, cov_i)
  
  cov_mod_names <- names(cov_mods) %>% str_remove(coll("_state"))
  
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

#write_csv(summary, paste0("tiger_mod_avg_coeffs_no_bs_site_eff2.csv"))

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
#ggsave(avg_dists_plot, file=paste0("tiger_mod_avg_coeffs_no_bs_site_eff2.png"), width = 12, height = 8)

