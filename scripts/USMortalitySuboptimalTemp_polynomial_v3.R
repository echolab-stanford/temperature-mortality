################################################################################
# Last Updated: 1/31/2025
# Name: Khusel Avirmed
# Purpose: 
#     1. Create multiple plots regarding the number of deaths attributable to temperatures
#     2. Heterogeneity studies by income, type of deaths, age, and gender
################################################################################

# clean environment
rm(list = ls())
gc()
  

# import used packages
library(fixest)
library(glue) # Expressions enclosed by braces will be evaluated as R code
library(readr) 
library(dplyr)
library(tidyr)
library(ggplot2)
library(MetBrewer)
library(DescTools)
library(plotrix) # weighted histogram function
library(doParallel)
library(stringr)
library(patchwork)
library(lubridate)
# install.packages("arrow")
library(arrow)
library(data.table)

# setwd("~/BurkeLab Dropbox/projects/adaptation")
country_name <- "US" # Other options: EU, MEX
script_dir <- 'script' # contains helper functions or other scripts needed
data_dir <- 'data' # parent data dir
raw_dir <- paste(data_dir, country_name, "raw", sep = "/") # this includes raw data 
cleaned_dir <- paste(data_dir, country_name, "cleaned", sep = "/")
processed_dir <- 'processed'  # parent mid analysis processed data such as main df, bootstrapped df, etc
processed_country_dir <- paste(processed_dir, country_name, sep = "/") # country specific processed data dir
fig_dir <- "fig"
fig_country_dir <- paste(fig_dir, country_name, sep = "/") # country specific processed data dir

# ---- Section 0: Script Loading ----
source(paste(script_dir, "run_bootstrap_feols.R", sep = "/"))
source(paste(script_dir, "predict_bstrap_at_values.R", sep = "/"))
# source(paste(script_dir, "mortality_cost_analysis.R", sep = "/"))

# ---- Section 1: Data Loading ----
# mort <- read_parquet(paste(cleaned_dir, 'USAllCause_age_standardized_rates_county_month_1968_2020_by_age_groups_complete.pq', sep = '/')  )
mort <- read_parquet(paste(cleaned_dir, paste0(country_name, "FullPanel_complete.pq"), sep = "/"))
load(paste(processed_country_dir, paste0(country_name, "AnnualBinned_person_day_exposure.RData"), sep = "/"))

mort <- mort %>%
  rename(mean_tm = avg_temp, mean_income = pc_income, temp_d1 = pw_tmean_nopoly) %>%  
  filter(year < 2020)  
setDT(mort)
str(mort)

# --------------
# Constants 
# -------------------
pop_weight <- NULL
temp_range <- -50:50
Nboot <- 100
options(warn = -1)
set.seed(94305)
setFixest_nthreads(8, save = F)

# ------ 1. Model selections: 
# --------------- 1. model formulas
# --------------- 2. model run once
# --------------- 3. model bootstrap

# pooled model formula 
str(mort)
pooled_model_formula <- as.formula("rate_all ~ 
                                   temp_d1_l0 +
                                   temp_d2_l0 + 
                                   temp_d3_l0 + 
                                   temp_d4_l0 + 
                                   precip_d1_l0 + 
                                   precip_d2_l0 |  adm_id^year + adm_id^month")
pooled_coefs <- feols(pooled_model_formula, data = mort, weights = NULL, only.coef = TRUE, mem.clean = T)

# time_taken <- system.time({
#   pooled_bstrap <- run_bootstrap_feols(
#     sample_from = mort,
#     unit_column = "adm_id",
#     Nboot = 100,
#     model_list = list(pooled_model_formula),
#     weights = NULL,
#     model_name = "Pooled",
#     country_name = country_name,
#     processed_country_dir = processed_country_dir,
#     override = T
#   )
# })
# 
# print(time_taken)

# Introducing lag
pooled_formula_wlag_30days <- as.formula("rate_all ~
                                   temp_d1_l0 + temp_d1_l.[1:30] +
                                   temp_d2_l0 + temp_d2_l.[1:30] +
                                   temp_d3_l0 + temp_d3_l.[1:30] +
                                   temp_d4_l0 + temp_d4_l.[1:30] +
                                   precip_d1_l0 + precip_d1_l.[1:30] +
                                   precip_d2_l0 + precip_d2_l.[1:30] |  adm_id^year + adm_id^month")

pooled_coefs_wlag_30days <- feols(pooled_formula_wlag_30days, data = mort, weights = NULL, only.coef = TRUE, nthreads = 8)

# system.time({
#   pooled_bstrap_wlag_30days <- run_bootstrap_feols(
#     sample_from = mort,
#     unit_column = "adm_id",
#     Nboot = Nboot,
#     model_list = list(pooled_formula_wlag_30days),
#     weights = NULL,
#     model_name = "Pooled (lag = 30)",
#     country_name = country_name,
#     processed_country_dir = processed_country_dir,
#     override = T
#   )
# })

# ----------- 2. We find mmt and mmt adjusted temps to predict mortality
# first, find mmt and scale temps
find_mmt_and_scale_pred_temp <- function(coefs, 
                                         observed_temps, 
                                         degree = 4, 
                                         temp_range = -50:50, 
                                         mmt_quantiles = c(0.1, 0.9), 
                                         additional_grep = NULL, 
                                         lag.include = F) {
  
  # coefs = pooled_coefs_wlag
  # observed_temps = mort$temp_d0_l0 # a vector of observed temp distribution
  # degree = 4
  # temp_range = temp_range
  # mmt_quantiles = c(0.1, 0.9) # restrict mmt within these percentiles of the observed temp distribution
  # additional_grep = NULL
  # lag.include = T
  
  # Create polynomial matrix for prediction
  pred_temps_poly <- poly(temp_range, degree = degree, raw = TRUE)
  
  # Extract temperature-related coefficients and compute predicted values
  if (is.null(additional_grep)) {
    temp_coefs <- coefs[grepl("temp", names(coefs))]
  } else {
    temp_coefs <- coefs[grepl(pattern = "temp", x = names(coefs)) & grepl(pattern = additional_grep, x = names(coefs))]
  }
  
  if (lag.include){
    # add estimated coefficients to each polynomials
    poly1_coefs <- sum(temp_coefs[grepl("temp_d1", names(temp_coefs))])
    poly2_coefs <- sum(temp_coefs[grepl("temp_d2", names(temp_coefs))])
    poly3_coefs <- sum(temp_coefs[grepl("temp_d3", names(temp_coefs))])
    poly4_coefs <- sum(temp_coefs[grepl("temp_d4", names(temp_coefs))])
    
    # Append new coefficients while preserving the original named vector format
    temp_coefs <- c(temp_d1 = poly1_coefs, temp_d2 = poly2_coefs, 
                    temp_d3 = poly3_coefs, temp_d4 = poly4_coefs)
    
    pred_vector <- c(pred_temps_poly %*% temp_coefs)
  } else {
    # pred temp polynomials * estimated ceofs
    pred_vector <- c(pred_temps_poly %*% temp_coefs)
    # par(mar = c(2, 2, 2, 2))
    # plot(pred_vector)
  }
  
  # Restrict MMT to be within the observed temperature distribution
  mmt_min <- quantile(observed_temps, mmt_quantiles[1], na.rm = TRUE)
  mmt_max <- quantile(observed_temps, mmt_quantiles[2], na.rm = TRUE)
  
  # Find potential MMT options and their predicted values
  mmt_options <- which(temp_range < mmt_max & temp_range > mmt_min)
  mmt_options_values <- pred_vector[mmt_options]
  
  # Find the temperature corresponding to the minimum predicted value
  mmt_index_in_options <- which.min(mmt_options_values)
  mmt <- temp_range[mmt_options][mmt_index_in_options]
  
  # Scale predicted temperatures around MMT
  pred_temps_mmt <- poly(mmt, degree = degree, raw = TRUE)
  pred_temps_scaled <- scale(pred_temps_poly, center = pred_temps_mmt, scale = FALSE)
  
  return(list(mmt = mmt, pred_temps_scaled = pred_temps_scaled))
}

# check if this works
pooled_mmt_scaled_temp <- find_mmt_and_scale_pred_temp(
  coefs = pooled_coefs, 
  observed_temps = mort$temp_d1, 
  degree = 4, 
  temp_range = -50:50, 
  mmt_quantiles = c(0.01, 0.9), 
  additional_grep = NULL,
  lag.include = F
) # looks ok

print(pooled_mmt_scaled_temp$mmt)
print(pooled_mmt_scaled_temp$pred_temps_scaled)

pooled_mmt_scaled_temp_wlag_30days <- find_mmt_and_scale_pred_temp(
  coefs = pooled_coefs_wlag_30days,
  observed_temps = mort$temp_d1,
  degree = 4,
  temp_range = -50:50,
  mmt_quantiles = c(0.01, 0.9),
  additional_grep = NULL,
  lag.include = T) # looks ok

# ----------- 3. Then, we use mmt and mmt adjusted temp to generate prediction and calculate absolute deaths
# ------------ for each bootstrap iteration
run_bootstrap_mortality_temp_analysis <- function(boots_results, 
                                                  pred_temps_scaled,
                                                  mmt,
                                                  temp_range,
                                                  pwd, 
                                                  model_name, 
                                                  model_level,
                                                  country_name,
                                                  processed_country_dir, 
                                                  coef_pattern = "temp", 
                                                  additional_grep = NULL, 
                                                  lag.include = F, 
                                                  save_results = T) {
  
  # boots_results = pooled_bstrap_wlag
  # pred_temps_scaled = pooled_mmt_scaled_temp_wlag$pred_temps_scaled
  # mmt <- pooled_mmt_scaled_temp_wlag$mmt
  # pwd <- pooled_pwd
  # temp_range = -50:50
  # model_name = "Pooled (lag)"
  # model_level = "pooled"
  # country_name = country_name
  # processed_country_dir = processed_country_dir
  # coef_pattern = "temp"
  # log.include = T
  # additional_grep = NULL
  
  # boots_results = pooled_bstrap
  # pred_temps_scaled = pooled_mmt_scaled_temp$pred_temps_scaled
  # mmt <- pooled_mmt_scaled_temp$mmt
  # pwd <- pooled_pwd
  # temp_range = -50:50
  # model_name = "Pooled"
  # model_level = "pooled"
  # country_name = country_name
  # processed_country_dir = processed_country_dir
  # coef_pattern = "temp"
  # additional_grep = NULL
  # lag.include = F
  
  # boots_results = bstrap_at_level
  # pred_temps_scaled = mmt_scaled_temp$pred_temps_scaled
  # mmt = mmt_scaled_temp$mmt
  # pwd = pwd_at_level
  # temp_range = temp_range
  # model_name = model_name
  # model_level = level_name
  # country_name = country_name
  # processed_country_dir = processed_country_dir
  # coef_pattern = "pw_tmean"
  # additional_grep = NULL
  
  # Extract temperature-related coefficients from bootstrap results
  # if (is.null(additional_grep)) {
  #   boots_coefs <- boots_results[grepl(pattern = coef_pattern, names(boots_results))]
  # } else {
  #   boots_coefs <- boots_results[grepl(pattern = coef_pattern, names(boots_results)) & grepl(pattern = additional_grep, names(boots_results)) ]
  # }
  
  if (lag.include){
    poly_names <- c("temp_d1", "temp_d2", "temp_d3", "temp_d4")
    boots_coefs <- sapply(poly_names, function(poly) {
      # poly <- "temp_d1"
      poly_coefs <- boots_results[, grepl(poly, names(boots_results)), drop = FALSE]
      rowSums(poly_coefs)
    })
    # Compute bootstrapped predictions: each row is an iteration
    boots_pred <- boots_coefs %*% t(pred_temps_scaled)
    
  } else {
    boots_coefs <- boots_results[grepl(pattern = coef_pattern, names(boots_results))]
    boots_coefs <- as.matrix(boots_coefs)
    
    # Compute bootstrapped predictions: each row is an iteration
    boots_pred <- boots_coefs %*% t(pred_temps_scaled)
  }
  
  # Compute mean and confidence intervals for each temperature
  response <- apply(boots_pred, 2, function(coeff) {
    pred_q025 <- quantile(coeff, probs = 0.025, na.rm = TRUE)
    pred <- mean(coeff, na.rm = TRUE)
    pred_q975 <- quantile(coeff, probs = 0.975, na.rm = TRUE)
    return(c(pred_q025 = pred_q025, pred = pred, pred_q975 = pred_q975))
  })
  
  # Format result annual_metrics
  response <- as.data.table(t(response))
  names(response) <- c("pred_q025", "pred", "pred_q975")
  response$bins <- temp_range
  response$mmt <- mmt
  response$n_iterations <- nrow(boots_pred)
  # response$model_name <- model_name
  # response$model_level <- model_level
  
  setDT(pwd)
  # person-days weight for plot
  # Ensure pwd is a data.table
  pwd[, bins := ceiling(as.numeric(bins))]  # Ensure bins are numeric
  pwd_plot <- pwd[
    , .(pop_deg_days = mean(pop_deg_days, na.rm = T), n_bins = .N), by = .(bins)][
      , pop_deg_days_total := sum(pop_deg_days)][
        , pop_deg_days := pop_deg_days / pop_deg_days_total][
          , .(bins, pop_deg_days)]
  
  par(mfrow = c(4, 1))  # 2 rows, 1 column
  plot(response$bins, response$pred, type = "l")
  plot(pwd_plot$bins, pwd_plot$pop_deg_days, type = "l")
  
  # merge person-days weights with response for plot later
  response <- merge(response, pwd_plot, by = "bins")
  response[, `:=`(model_name = model_name, model_level = model_level)]
  
  # for each bootstrap iteration 
  # Set up parallel backend
  cl <- makeCluster(detectCores() - 1)  # Use available cores minus 1
  registerDoParallel(cl)
  
  # indicator whether collapse by adm units
  collapse.adm.unit <- if (model_name == "temp" | model_name == "income") T else F
  
  
  # pwd <- setDT(pwd)
  # Parallel loop
  results <- foreach(i = 1:nrow(boots_pred), 
                     .combine = "rbind", 
                     .packages = c("data.table")) %dopar% {
                       
                       # Select iteration
                       iteration_index <- i
                       boots <- boots_pred[iteration_index, ]  
                       
                       # Create iteration table
                       iteration_df <- data.table(pred = boots, bins = temp_range, mmt = mmt)
                       
                       # Merge with pwd
                       pwd_year_pred <- merge(pwd, iteration_df, by = "bins", all.x = TRUE)
                       
                       # Weight pred by pop-days exposure
                       pwd_year_pred[, deaths := pop_deg_days * pred]
                       pwd_year_pred[, iteration_index := iteration_index]
                       pwd_year_pred[, side := fifelse(bins < mmt, "cold", "hot")]
                       
                       # Collapse by administrative unit if required
                       if (collapse.adm.unit) {
                         pwd_year_pred <- pwd_year_pred[, .(
                           pop_deg_days = sum(pop_deg_days, na.rm = TRUE),
                           annual_total_pop = sum(annual_total_pop, na.rm = TRUE),
                           annual_total_death = sum(annual_total_death, na.rm = TRUE),
                           deaths = sum(deaths, na.rm = TRUE),
                           pred = first(pred),
                           mmt = first(mmt),
                           iteration_index = first(iteration_index),
                           side = first(side),
                           n_adm = .N
                         ), by = .(bins, model_level, year)]
                       }
                       
                       return(pwd_year_pred)
                     }
  
  # Stop parallel cluster
  stopCluster(cl)
  
  # deaths per bins
  deaths_with_ci <- results[, .(
    deaths_q025 = quantile(deaths, probs = 0.025, na.rm = TRUE),
    deaths      = mean(deaths, na.rm = TRUE),
    deaths_q975 = quantile(deaths, probs = 0.975, na.rm = TRUE),
    mmt         = first(mmt),
    n_iterations = .N
  ), by = .(year, bins)]
  
  # attach model name and level
  deaths_with_ci[, `:=`(model_name = model_name, model_level = model_level)]
  
  deaths_with_ci_plot <- deaths_with_ci[, .(deaths = mean(deaths, na.rm = T), n_years = .N), by = .(bins)]
  plot(deaths_with_ci_plot$bins, deaths_with_ci_plot$deaths, type = "l")
  # plot(deaths_with_ci[deaths_with_ci$year == 2000, ]$deaths, type = "l")
  
  # calculate deaths related to cold or hot temperature
  deaths_temp_related_by_bstrap <- results[, .(
    deaths = sum(deaths, na.rm = T),
    mmt = first(mmt),
    n_bins = .N
  ), by = .(year, side, iteration_index)]
  deaths_temp_related_by_bstrap[, `:=`(model_name = model_name, model_level = model_level)]
  
  # calculate deaths related to cold or hot with confidence interval
  deaths_temp_related_with_ci <- deaths_temp_related_by_bstrap[, .(
    deaths_q025 = quantile(deaths, probs = 0.025, na.rm = TRUE),
    deaths = mean(deaths, na.rm = TRUE),
    deaths_q975 = quantile(deaths, probs = 0.975, na.rm = TRUE),
    mmt = first(mmt),
    n_iterations = .N
  ), by = .(year, side)]
  
  # NOTE: fix this
  if (collapse.adm.unit) {
    annual_total_pop_death <- pwd[, .(
      annual_total_pop = sum(annual_total_pop, na.rm = TRUE),
      annual_total_death = sum(annual_total_death, na.rm = TRUE),
      n_adm = .N
    ), by = .(bins, model_level, year)]
    
    annual_total_pop_death <- unique(annual_total_pop_death[, .(year, annual_total_pop, annual_total_death)])
  } else {
    annual_total_pop_death <- unique(pwd[, .(year, annual_total_pop, annual_total_death)])
  }
  
  deaths_temp_related_with_ci <- merge(deaths_temp_related_with_ci, annual_total_pop_death, by = "year", allow.cartesian = TRUE)
  deaths_temp_related_with_ci[, `:=`(model_name = model_name, model_level = model_level)]
  
  # ADD: annual deaths
  annual_deaths_temp_related_with_ci <- deaths_temp_related_with_ci[, .(
    deaths_q025 = mean(deaths_q025, na.rm = TRUE),
    deaths = mean(deaths, na.rm = TRUE),
    deaths_q975 = mean(deaths_q975, na.rm = TRUE),
    mmt = first(mmt),
    n_iterations = first(n_iterations),
    n_years = .N), by = .(side)]
  annual_deaths_temp_related_with_ci[, `:=`(model_name = model_name, model_level = model_level)]
  
  barplot(annual_deaths_temp_related_with_ci$deaths)
  
  par(mfrow = c(1, 1))  # Reset to default 1 plot layout
  
  if (save_results) {
    # Save all data frames
    print("Saving results ... ")
    
    response_by_bins_dir <- file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_", tolower(model_name), "_response_by_bins.pq"))
    deaths_by_bins_dir <- file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_", tolower(model_name), "_deaths_by_bins.pq"))
    deaths_by_side_year_iteration <- file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_", tolower(model_name), "_deaths_by_side_year_iteration.pq"))
    deaths_by_side_year <- file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_", tolower(model_name),"_deaths_by_side_year.pq"))
    
    # Save the data frames
    write_parquet(response, response_by_bins_dir)
    write_parquet(deaths_with_ci, deaths_by_bins_dir)
    write_parquet(deaths_temp_related_by_bstrap, deaths_by_side_year_iteration)
    write_parquet(deaths_temp_related_with_ci, deaths_by_side_year)
    
    print(paste0("Results are saved successfully in the directory: ", processed_country_dir))
  }
  
  return(list(
    response_df = response, 
    deaths_by_bins = deaths_with_ci,
    deaths_by_side_iteration = deaths_temp_related_by_bstrap,
    deaths_annual_by_side = deaths_temp_related_with_ci
  ))
}


# run Pooled model
# pooled_bstrap <- as.data.frame(t(pooled_coefs))
pooled_bstrap <- read_parquet( paste(processed_country_dir, paste0(country_name, "_BootstrapResult_pooled_n100.pq"), sep = "/") )
pooled_bstrap <- as.data.frame(pooled_bstrap)

pooled_result_list <- run_bootstrap_mortality_temp_analysis(
  boots_results = pooled_bstrap,
  pred_temps_scaled = pooled_mmt_scaled_temp$pred_temps_scaled,
  mmt <- pooled_mmt_scaled_temp$mmt,
  pwd <- pooled_pwd,
  temp_range = -50:50,
  model_name = "Pooled",
  model_level = "pooled",
  country_name = country_name,
  processed_country_dir = processed_country_dir,
  coef_pattern = "temp",
  additional_grep = NULL,
  lag.include = F,
  save_results = T
)

pooled_response_df <- pooled_result_list$response_df
pooled_deaths_by_bins <- pooled_result_list$deaths_by_bins
pooled_deaths_by_side_iteration <- pooled_result_list$deaths_by_side_iteration
pooled_deaths_annual_by_side <- pooled_result_list$deaths_annual_by_side

# pooled w/ lag = 30
pooled_bstrap_wlag_30days <- as.data.frame(t(pooled_coefs_wlag_30days))
pooled_result_list_wlag_30days <- run_bootstrap_mortality_temp_analysis(
  boots_results = pooled_bstrap_wlag_30days,
  pred_temps_scaled = pooled_mmt_scaled_temp_wlag_30days$pred_temps_scaled,
  mmt <- pooled_mmt_scaled_temp_wlag_30days$mmt,
  pwd <- pooled_pwd,
  temp_range = -50:50,
  model_name = "Pooled (Lag = 30)",
  model_level = "pooled",
  country_name = country_name,
  processed_country_dir = processed_country_dir,
  coef_pattern = "temp",
  additional_grep = NULL,
  lag.include = T,
  save_results = T
)

pooled_response_df_wlag_30days <- pooled_result_list_wlag_30days$response_df
pooled_deaths_by_bins_wlag_30days <- pooled_result_list_wlag_30days$deaths_by_bins
pooled_deaths_by_side_iteration_wlag_30days <- pooled_result_list_wlag_30days$deaths_by_side_iteration
pooled_deaths_annual_by_side_wlag_30days <- pooled_result_list_wlag_30days$deaths_annual_by_side

# run Young model
# young_model_formula <- as.formula("rate_young ~ 
#                                    temp_d1_l0 + temp_d1_l.[1:30] +
#                                    temp_d2_l0 + temp_d2_l.[1:30] +
#                                    temp_d3_l0 + temp_d3_l.[1:30] +
#                                    temp_d4_l0 + temp_d4_l.[1:30] +
#                                    precip_d1_l0 + precip_d1_l.[1:30] +
#                                    precip_d2_l0 + precip_d2_l.[1:30] |  adm_id^year + adm_id^month")
young_model_formula <- as.formula("rate_young ~                                   
                                   temp_d1_l0 +
                                   temp_d2_l0 + 
                                   temp_d3_l0 + 
                                   temp_d4_l0 + 
                                   precip_d1_l0 + 
                                   precip_d2_l0 |  adm_id^year + adm_id^month")
young_coefs <- feols(young_model_formula, data = mort, weights = NULL, only.coef = T)

# young_bstrap <- run_bootstrap_feols(
#   sample_from = mort, 
#   unit_column = "adm_id", 
#   Nboot = Nboot, 
#   model_list = list(young_model_formula), 
#   weights = NULL,
#   model_name = "Young",
#   country_name = country_name,
#   processed_country_dir = processed_country_dir,
#   override = T
# )

# run Young model
young_mmt_scaled_temp <- find_mmt_and_scale_pred_temp(
  coefs = young_coefs, 
  observed_temps = mort$temp_d1, 
  degree = 4, 
  temp_range = -50:50, 
  mmt_quantiles = c(0.01, 0.9), 
  additional_grep = NULL,
  lag.include = T
)

print(young_mmt_scaled_temp$mmt)
print(young_mmt_scaled_temp$pred_temps_scaled)

# young_bstrap <- as.data.frame(t(young_coefs))
young_bstrap <- read_parquet( paste(processed_country_dir, paste0(country_name, "_BootstrapResult_young_n100.pq"), sep = "/") )
young_bstrap <- as.data.frame(young_bstrap)

young_result_list <- run_bootstrap_mortality_temp_analysis(
  boots_results = young_bstrap,
  pred_temps_scaled = young_mmt_scaled_temp$pred_temps_scaled,
  mmt <- young_mmt_scaled_temp$mmt,
  pwd <- young_pwd,
  temp_range = -50:50,
  model_name = "Young",
  model_level = "pooled",
  country_name = country_name,
  processed_country_dir = processed_country_dir,
  coef_pattern = "temp",
  additional_grep = NULL,
  save_results = T,
  lag.include = T
)

young_response_df <- young_result_list$response_df
young_deaths_by_bins <- young_result_list$deaths_by_bins
young_deaths_by_side_iteration <- young_result_list$deaths_by_side_iteration
young_deaths_annual_by_side <- young_result_list$deaths_annual_by_side

# run Adult model
# adult_model_formula <- as.formula("rate_adult ~ 
#                                    temp_d1_l0 + temp_d1_l.[1:30] +
#                                    temp_d2_l0 + temp_d2_l.[1:30] +
#                                    temp_d3_l0 + temp_d3_l.[1:30] +
#                                    temp_d4_l0 + temp_d4_l.[1:30] +
#                                    precip_d1_l0 + precip_d1_l.[1:30] +
#                                    precip_d2_l0 + precip_d2_l.[1:30] |  adm_id^year + adm_id^month")
adult_model_formula <- as.formula("rate_adult ~ 
                                   temp_d1_l0 +
                                   temp_d2_l0 + 
                                   temp_d3_l0 + 
                                   temp_d4_l0 + 
                                   precip_d1_l0 + 
                                   precip_d2_l0 |  adm_id^year + adm_id^month")
adult_coefs <- feols(adult_model_formula, data = mort, weights = NULL, only.coef = T)
# adult_bstrap <- run_bootstrap_feols(
#   sample_from = mort, 
#   unit_column = "adm_id", 
#   Nboot = Nboot, 
#   model_list = list(adult_model_formula), 
#   weights = NULL,
#   model_name = "Adult",
#   country_name = country_name,
#   processed_country_dir = processed_country_dir,
#   override = T
# )

# run Adult model
adult_mmt_scaled_temp <- find_mmt_and_scale_pred_temp(
  coefs = adult_coefs, 
  observed_temps = mort$temp_d1, 
  degree = 4, 
  temp_range = -50:50, 
  mmt_quantiles = c(0.01, 0.9), 
  additional_grep = NULL, 
  lag.include = T
)

print(adult_mmt_scaled_temp$mmt)
print(adult_mmt_scaled_temp$pred_temps_scaled)

# adult_bstrap <- as.data.frame(t(adult_coefs))
adult_bstrap <- read_parquet( paste(processed_country_dir, paste0(country_name, "_BootstrapResult_adult_n100.pq"), sep = "/") )
adult_bstrap <- as.data.frame(adult_bstrap)

adult_result_list<- run_bootstrap_mortality_temp_analysis(
  boots_results = adult_bstrap,
  pred_temps_scaled = adult_mmt_scaled_temp$pred_temps_scaled,
  mmt <- adult_mmt_scaled_temp$mmt,
  pwd <- adult_pwd,
  temp_range = -50:50,
  model_name = "Adult",
  model_level = "pooled",
  country_name = country_name,
  processed_country_dir = processed_country_dir,
  coef_pattern = "temp",
  additional_grep = NULL,
  save_results = T,
  lag.include = T
)

adult_response_df <- adult_result_list$response_df
adult_deaths_by_bins <- adult_result_list$deaths_by_bins
adult_deaths_by_side_iteration <- adult_result_list$deaths_by_side_iteration
adult_deaths_annual_by_side <- adult_result_list$deaths_annual_by_side

# run Elderly model
# elderly_model_formula <- as.formula("rate_elderly ~ 
#                                    temp_d1_l0 + temp_d1_l.[1:30] +
#                                    temp_d2_l0 + temp_d2_l.[1:30] +
#                                    temp_d3_l0 + temp_d3_l.[1:30] +
#                                    temp_d4_l0 + temp_d4_l.[1:30] +
#                                    precip_d1_l0 + precip_d1_l.[1:30] +
#                                    precip_d2_l0 + precip_d2_l.[1:30] |  adm_id^year + adm_id^month")

elderly_model_formula <- as.formula("rate_elderly ~                                    
                                   temp_d1_l0 +
                                   temp_d2_l0 + 
                                   temp_d3_l0 + 
                                   temp_d4_l0 + 
                                   precip_d1_l0 + 
                                   precip_d2_l0 |  adm_id^year + adm_id^month")
elderly_coefs <- feols(elderly_model_formula, data = mort, weights = NULL, only.coef = T)
# elderly_bstrap <- run_bootstrap_feols(
#   sample_from = mort, 
#   unit_column = "adm_id", 
#   Nboot = Nboot, 
#   model_list = list(elderly_model_formula), 
#   weights = NULL,
#   model_name = "Elderly",
#   country_name = country_name,
#   processed_country_dir = processed_country_dir,
#   override = T
# )

# run Elderly model
elderly_mmt_scaled_temp <- find_mmt_and_scale_pred_temp(
  coefs = elderly_coefs, 
  observed_temps = mort$temp_d1, 
  degree = 4, 
  temp_range = -50:50, 
  mmt_quantiles = c(0.01, 0.9), 
  additional_grep = NULL,
  lag.include = T
)

print(elderly_mmt_scaled_temp$mmt)
print(elderly_mmt_scaled_temp$pred_temps_scaled)

# elderly_bstrap <- as.data.frame(t(elderly_coefs))
elderly_bstrap <- read_parquet( paste(processed_country_dir, paste0(country_name, "_BootstrapResult_elderly_n100.pq"), sep = "/") )
elderly_bstrap <- as.data.frame(elderly_bstrap)

elderly_result_list <- run_bootstrap_mortality_temp_analysis(
  boots_results = elderly_bstrap,
  pred_temps_scaled = elderly_mmt_scaled_temp$pred_temps_scaled,
  mmt <- elderly_mmt_scaled_temp$mmt,
  pwd <- elderly_pwd,
  temp_range = -50:50,
  model_name = "Elderly",
  model_level = "pooled",
  country_name = country_name,
  processed_country_dir = processed_country_dir,
  coef_pattern = "temp",
  additional_grep = NULL,
  save_results = T,
  lag.include = T
)

elderly_response_df <- elderly_result_list$response_df
elderly_deaths_by_bins <- elderly_result_list$deaths_by_bins
elderly_deaths_by_side_iteration <- elderly_result_list$deaths_by_side_iteration
elderly_deaths_annual_by_side <- elderly_result_list$deaths_annual_by_side

# run Temp Continuous Interaction model
# for continuous interaction, we have to evaluate at mean continuous variable before generating prediction
# for temp is varying by adm unit and year
# so, we generate prediction at each unit at each percentile
# --------- 1. we run model once and get mmt and mmt adjusted temp to predict response at each percentile bins
# --------- 2. for each percentile, we evaluate at that respective percent bins value 
# --------- 3. using estimated and evaluated at continuous percentile bins coeffs, we generate predictions for each adm units for that percent bins
# --------- 4. we returns average across units results for that percent bins
# --------- 5. we repeat same sequence for each bootstrap
run_parallel_mortality_analysis <- function(mort, pwd, coefs, boots_results, temp_range, model_name, country_name, processed_country_dir, coef_pattern = "temp", save_results = F, lag.include = F) {
  library(data.table)
  library(doParallel)
  library(foreach)
  
  # mort = mort
  # pwd = adm_perc_temp_pwd
  # coefs = temp_coefs
  # boots_results = temp_bstrap
  # temp_range = -50:50
  # model_name = "temp"
  # country_name = country_name
  # coef_pattern = "temp"
  # processed_country_dir = processed_country_dir
  
  # mort = mort
  # pwd = adm_perc_income_pwd
  # coefs = income_coefs
  # boots_results = income_bstrap
  # temp_range = -50:50
  # model_name = "income"
  # country_name = country_name
  # processed_country_dir = processed_country_dir
  
  # mort = mort
  # pwd = decade_pwd
  # coefs = decade_coefs
  # boots_results = decade_bstrap
  # temp_range = -50:50
  # model_name = "decade"
  # country_name = country_name
  # processed_country_dir = processed_country_dir
  
  setDT(mort)
  setDT(pwd)
  
  # NOTE: Fix quantiles .... 
  
  if (model_name == "temp") {
    # quantiles <- sort(unique(pwd$temp_percentile))
    
    quantiles <- quantile(mort$mean_tm, probs = seq(0.01, 1, 0.01), na.rm = TRUE)
    quantiles_names <- seq(0, 99, 1)
    names(quantiles) <- quantiles_names
    
    grep_name <- "mean_tm"
    observed_temp_df <- mort[, .(temp_percentile, temp_d1)]
    setnames(observed_temp_df, old = "temp_percentile", new = "model_level")
    setnames(pwd, old = "temp_percentile", new = "model_level")
    
    # adm_year_bin <- mort[, .(adm_id, year, tperc_bins)]
    # pwd <- merge(pwd, adm_year_bin, by = c("adm_id", "year"), allow.cartesian = T)
    # setnames(pwd, old = "tperc_bins", new = "model_level")
    
  } else if (model_name == "income") {
    
    quantiles <- quantile(mort$mean_income, probs = seq(0.01, 1, 0.01), na.rm = TRUE)
    quantiles_names <- seq(0, 99, 1)
    names(quantiles) <- quantiles_names
    
    # quantiles <- quantile(mort$mean_tm, probs = seq(0.01, 1, 0.01), na.rm = TRUE)
    # quantiles_names <- paste0(seq(0, 99, 1), "%-", seq(1, 100, 1), "%")
    # names(quantiles) <- quantiles_names
    
    grep_name <- "mean_income"
    observed_temp_df <- mort[, .(income_percentile, temp_d1)]
    setnames(observed_temp_df, old = "income_percentile", new = "model_level")
    setnames(pwd, old = "income_percentile", new = "model_level")
    
    
  } else if (model_name == "decade") {
    quantiles <- as.numeric(as.character(unique(mort$decade)))
    names(quantiles) <- quantiles
    
    observed_temp_df <- mort[, .(decade, temp_d1)]
    setnames(observed_temp_df, old = "decade", new = "model_level")
    setnames(pwd, old = "decade", new = "model_level")
  }
  
  # pwd[, model_level := as.character(model_level)]
  
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  
  evaluate_at_continuous_var <- function(coefs, at_values, grep_name) {
    nth_evaluation <- coefs[grepl(grep_name, names(coefs)) & grepl(coef_pattern, names(coefs))] * at_values
    intercepts <- coefs[grepl(coef_pattern, names(coefs)) & !grepl(grep_name, names(coefs))]
    return(nth_evaluation + intercepts)
  }
  
  results_list <- foreach(i = 1:length(names(quantiles)), 
                          .combine = "rbind", 
                          .packages = c("data.table", "doParallel", "foreach"),
                          .export = c("find_mmt_and_scale_pred_temp", "evaluate_at_continuous_var", 
                                      "run_bootstrap_mortality_temp_analysis")) %dopar% {
                                        # i <- 1
                                        level_name <- names(quantiles)[i]
                                        level_value <- quantiles[[level_name]]
                                        
                                        # level_value <- log(level_value)
                                        if (model_name == "decade") {
                                          grep_name <- level_value
                                          coefs_at_level = coefs[grepl(pattern = coef_pattern, x = names(coefs)) & grepl(pattern = grep_name, x = names(coefs))]
                                          bstrap_at_level <- boots_results[grepl(pattern = coef_pattern, x = names(boots_results)) & grepl(pattern = grep_name, x = names(boots_results))]
                                        } else {
                                          coefs_at_level <- evaluate_at_continuous_var(coefs, level_value, grep_name)
                                          bstrap_at_level <- evaluate_at_continuous_var(boots_results, level_value, grep_name)
                                        }
                                        
                                        # get observed temp
                                        observed_temps <- observed_temp_df[model_level == level_name, temp_d1]
                                        
                                        mmt_scaled_temp <- find_mmt_and_scale_pred_temp(
                                          coefs = coefs_at_level, 
                                          observed_temps = observed_temps, 
                                          degree = 4, 
                                          temp_range = temp_range, 
                                          mmt_quantiles = c(0.01, 0.9), 
                                          additional_grep = NULL,
                                          lag.include = lag.include
                                        )
                                        
                                        # quantile(observed_temps, prob = c(0.01, 0.9))
                                        
                                        # hist(mort$mean_income)
                                        # hist(log(mort$mean_income))
                                        
                                        # subset person days weights to only that level
                                        pwd_at_level <- pwd[model_level == level_name, ]
                                        
                                        all_result <- run_bootstrap_mortality_temp_analysis(
                                          boots_results = bstrap_at_level,
                                          pred_temps_scaled = mmt_scaled_temp$pred_temps_scaled,
                                          mmt = mmt_scaled_temp$mmt,
                                          pwd = pwd_at_level,
                                          temp_range = temp_range,
                                          model_name = model_name,
                                          model_level = level_name,
                                          country_name = country_name,
                                          processed_country_dir = processed_country_dir,
                                          coef_pattern = coef_pattern,
                                          additional_grep = NULL,
                                          save_results = F,
                                          lag.include = lag.include
                                        )
                                        
                                        # all_result$response_df %>% plot_response(-50, 50)
                                        return(all_result)
                                      }
  
  stopCluster(cl)
  
  n_quantiles <- length(names(quantiles))  # Get the number of quantiles dynamically
  
  # Combine data frames dynamically based on the number of quantiles
  response <- do.call(rbind, results_list[1:n_quantiles])
  deaths_with_ci <- do.call(rbind, results_list[(n_quantiles + 1):(2 * n_quantiles)])
  deaths_temp_related_by_bstrap <- do.call(rbind, results_list[(2 * n_quantiles + 1):(3 * n_quantiles)])
  deaths_temp_related_with_ci <- do.call(rbind, results_list[(3 * n_quantiles + 1):(4 * n_quantiles)])
  
  if (save_results) {
    # Save all data frames
    print("Saving results ... ")
    
    response_by_bins_dir <- file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_", tolower(model_name), "_response_by_bins.pq"))
    deaths_by_bins_dir <- file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_", tolower(model_name), "_deaths_by_bins.pq"))
    deaths_by_side_year_iteration <- file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_", tolower(model_name), "_deaths_by_side_year_iteration.pq"))
    deaths_by_side_year <- file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_", tolower(model_name),"_deaths_by_side_year.pq"))
    
    # Save the data frames
    write_parquet(response, response_by_bins_dir)
    write_parquet(deaths_with_ci, deaths_by_bins_dir)
    write_parquet(deaths_temp_related_by_bstrap, deaths_by_side_year_iteration)
    write_parquet(deaths_temp_related_with_ci, deaths_by_side_year)
    
    print(paste0("Results are saved successfully in the directory: ", processed_country_dir))
  }
  
  return(list(
    response_df = response, 
    deaths_by_bins = deaths_with_ci,
    deaths_by_side_iteration = deaths_temp_related_by_bstrap,
    deaths_annual_by_side = deaths_temp_related_with_ci
  ))
}

# run Temp Continuous Interaction model
# names(mort)
# temp_lags <- names(mort)[grepl("temp_d", names(mort)) & grepl("l", names(mort))]
# precip_lags <- names(mort)[grepl("precip_d", names(mort)) & grepl("l", names(mort))]
# interaction_vars <- paste(c(temp_lags, precip_lags), collapse = " + ")
# interaction_term_temp <- paste0("(", interaction_vars, "):mean_tm")
# fixed_effects <- "adm_id^year + adm_id^month"
# temp_model_formula <- as.formula(paste("rate_all ~", interaction_vars, "+", interaction_term_temp, "|", fixed_effects))

temp_model_formula <- as.formula("rate_all ~ temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0 + precip_d1_l0 + precip_d2_l0  + 
                                 (temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0 + precip_d1_l0 + precip_d2_l0):mean_tm | adm_id^year + adm_id^month")


temp_coefs <- feols(temp_model_formula, data = mort, weights = NULL, only.coef = TRUE, lean = T)

# temp_bstrap <- run_bootstrap_feols(
#   sample_from = mort, 
#   unit_column = "adm_id", 
#   Nboot = Nboot, 
#   model_list = list(temp_model_formula), 
#   weights = NULL,
#   model_name = "temp",
#   country_name = country_name,
#   processed_country_dir = processed_country_dir,
#   override = T
# )

# temp_bstrap <- as.data.frame(t(temp_coefs))
temp_bstrap <- read_parquet( paste(processed_country_dir, paste0(country_name, "_BootstrapResult_temp_n100.pq"), sep = "/") )
temp_bstrap <- as.data.frame(temp_bstrap)

temp_results_list <- run_parallel_mortality_analysis(
  mort = mort, 
  pwd = adm_perc_temp_pwd, 
  coefs = temp_coefs, 
  boots_results = temp_bstrap, 
  temp_range = -50:50, 
  model_name = "temp", 
  country_name = country_name, 
  processed_country_dir = processed_country_dir, 
  coef_pattern = "temp",
  save_results = T,
  lag.include = T
)

temp_response_df <- temp_results_list$response_df
temp_deaths_by_bins <- temp_results_list$deaths_by_bins
temp_deaths_by_side_iteration <- temp_results_list$deaths_by_side_iteration
temp_deaths_annual_by_side <- temp_results_list$deaths_annual_by_side

# run Income Continuous Interaction model
# names(adm_perc_income_pwd) # total pop and death for unit and percent bins
# interaction_term_income <- paste0("(", interaction_vars, "):mean_income")
# income_model_formula <- as.formula(paste("rate_all ~", interaction_vars, "+", interaction_term_income, "|", fixed_effects))

income_model_formula <- as.formula("rate_all ~ temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0 + precip_d1_l0 + precip_d2_l0  + 
                                 (temp_d1_l0 + temp_d2_l0 + temp_d3_l0 + temp_d4_l0 + precip_d1_l0 + precip_d2_l0):mean_income | adm_id^year + adm_id^month")

income_coefs <- feols(income_model_formula, data = mort, weights = NULL, only.coef = TRUE)
# income_bstrap <- run_bootstrap_feols(
#   sample_from = mort, 
#   unit_column = "adm_id", 
#   Nboot = Nboot, 
#   model_list = list(income_model_formula), 
#   weights = NULL,
#   model_name = "income",
#   country_name = country_name,
#   processed_country_dir = processed_country_dir,
#   override = T
# )

# run Income Continuous Interaction model
# income_bstrap <- as.data.frame(t(income_coefs))
income_bstrap <- read_parquet( paste(processed_country_dir, paste0(country_name, "_BootstrapResult_income_n100.pq"), sep = "/") )
income_bstrap <- as.data.frame(income_bstrap)

income_results_list <- run_parallel_mortality_analysis(
  mort = mort, 
  pwd = adm_perc_income_pwd, 
  coefs = income_coefs, 
  boots_results = income_bstrap, 
  temp_range = -50:50, 
  model_name = "income", 
  country_name = country_name, 
  processed_country_dir = processed_country_dir,
  coef_pattern = "temp",
  save_results = T,
  lag.include = T
)

income_response_df <- income_results_list$response_df
income_deaths_by_bins <- income_results_list$deaths_by_bins
income_deaths_by_side_iteration <- income_results_list$deaths_by_side_iteration
income_deaths_annual_by_side <- income_results_list$deaths_annual_by_side

# run Decade model
str(mort)
names(mort)

# temp_interactions <- unlist(lapply(temp_lags, function(v) {
#   paste0("i(decade,", v, ")")
# }))
# 
# precip_interactions <- unlist(lapply(precip_lags, function(v) {
#   paste0("i(decade,", v, ")")
# }))
# 
# decade_model_formula <- as.formula(paste("rate_all ~", paste(c(temp_interactions, precip_interactions), collapse = " + "), "|", fixed_effects))
decade_model_formula <- as.formula("rate_all ~ i(decade, temp_d1_l0) + i(decade, temp_d2_l0) +i(decade, temp_d3_l0) + i(decade, temp_d4_l0) + i(decade, precip_d1_l0) + i(decade, precip_d2_l0) | adm_id^month + adm_id^year")

decade_coefs <- feols(decade_model_formula, data = mort, weights = NULL, only.coef = T)
# decade_bstrap <- run_bootstrap_feols(
#   sample_from = mort, 
#   unit_column = "adm_id", 
#   Nboot = Nboot, 
#   model_list = list(decade_model_formula), 
#   weights = NULL,
#   model_name = "decade",
#   country_name = country_name,
#   processed_country_dir = processed_country_dir,
#   override = T
# )
# decade_bstrap <- as.data.frame(t(decade_coefs))
decade_bstrap <- read_parquet( paste(processed_country_dir, paste0(country_name, "_BootstrapResult_decade_n100.pq"), sep = "/") )
decade_bstrap <- as.data.frame(decade_bstrap)

decade_results_list <- run_parallel_mortality_analysis(
  mort = mort, 
  pwd = decade_pwd, 
  coefs = decade_coefs, 
  boots_results = decade_bstrap, 
  temp_range = -50:50, 
  model_name = "decade", 
  country_name = country_name, 
  coef_pattern = "temp",
  processed_country_dir = processed_country_dir,
  save_results = T,
  lag.include = T
)

decade_response_df <- decade_results_list$response_df
decade_deaths_by_bins <- decade_results_list$deaths_by_bins
decade_deaths_by_side_iteration <- decade_results_list$deaths_by_side_iteration
decade_deaths_annual_by_side <- decade_results_list$deaths_annual_by_side

max(pooled_response_df$pred)


# ------- 4. Plot all results together
pooled_response_df <- read_parquet(file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_pooled_response_by_bins.pq")))
pooled_response_df_wlag_30days <- read_parquet(file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_pooled (lag = 30)_response_by_bins.pq")))
young_response_df <- read_parquet(file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_young_response_by_bins.pq")))
adult_response_df <- read_parquet(file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_adult_response_by_bins.pq")))
elderly_response_df <- read_parquet(file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_elderly_response_by_bins.pq")))
temp_response_df <- read_parquet(file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_temp_response_by_bins.pq")))
income_response_df <- read_parquet(file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_income_response_by_bins.pq")))
decade_response_df <- read_parquet(file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_decade_response_by_bins.pq")))

pooled_deaths_annual_by_side <- read_parquet(file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_pooled_deaths_by_side_year.pq")))
pooled_deaths_annual_by_side_wlag_30days <- read_parquet(file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_pooled (lag = 30)_deaths_by_side_year.pq")))
young_deaths_annual_by_side <- read_parquet(file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_young_deaths_by_side_year.pq")))
adult_deaths_annual_by_side <- read_parquet(file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_adult_deaths_by_side_year.pq")))
elderly_deaths_annual_by_side <- read_parquet(file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_elderly_deaths_by_side_year.pq")))
temp_deaths_annual_by_side <- read_parquet(file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_temp_deaths_by_side_year.pq")))
income_deaths_annual_by_side <- read_parquet(file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_income_deaths_by_side_year.pq")))
decade_deaths_annual_by_side <- read_parquet(file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_decade_deaths_by_side_year.pq")))


max_pop <- max(pooled_response_df$pop_deg_days,
               young_response_df$pop_deg_days, 
               adult_response_df$pop_deg_days, 
               elderly_response_df$pop_deg_days, 
               temp_response_df$pop_deg_days, 
               income_response_df$pop_deg_days, 
               decade_response_df$pop_deg_days)
max_pop <- max_pop + 0.02



plot_response <- function(response_df, min_bins, max_bins, max_pop) {
  
  # response_df <- rbind(young_response_df, 
  #                      adult_response_df, 
  #                      elderly_response_df) %>% 
  #   select(-model_level) %>% 
  #   rename(model_level = model_name) %>% 
  #   mutate(model_level = as.factor(model_level)) %>% 
  #   mutate(model_level = factor(model_level, levels = c("Young", "Adult", "Elderly"))) 
  # 
  # min_bins = 0
  # max_bins = 40
  # 
  # max_pop <- max(response_df$pop_deg_days)
  # 
  # choose a color palatte
  color_palette <- rev(met.brewer("Homer1", length(unique(response_df$model_level))))
  
  # max_pred <- max(response_df$pred)
  response_plot_df <- response_df %>% 
    group_by(model_level) %>% 
    mutate(mmt = bins[which(pred == 0)]) %>% 
    ungroup() %>% 
    mutate(model_level = paste0(model_level, " (MMT = ", mmt, "°C)"),
           model_level = as.factor(model_level)) %>% 
    filter(bins > min_bins & bins < max_bins) %>%
    mutate(pop_deg_days = pop_deg_days * max_pop)
  
  response_plot <- response_plot_df %>% 
    ggplot(aes(x = bins)) +
    
    # Plot pred values
    geom_line(aes(y = pred, color = model_level), linewidth = 1) +
    geom_ribbon(aes(ymin = pred_q025, ymax = pred_q975, fill = model_level), alpha = 0.2) +
    
    # Add hline at y=0
    geom_hline(yintercept = 0, linetype = "dashed") +
    
    # # Plot weights (scaled for visibility) with color mapping; use identity position to overlap
    # geom_col(aes(y = -pop_deg_days, fill = model_level), alpha = 0.5, position = "identity", width = 0.9) +
    
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) +
    guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) + 
    labs(color = NULL, x = "Temperature (°C)", y = "Change in deaths") +
    theme_minimal() +
    theme(
      strip.text   = element_text(size = 10, face = "bold"), 
      legend.position = c(.40, .80),
      axis.title.y = element_blank(),
      # axis.text.y  = element_blank(),
      # axis.ticks.y = element_blank(),
      axis.text.x  = element_blank(),
      axis.title.x = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      plot.margin = margin(0, 0, 0, 0)
    )
  
  exposure_plot <- response_plot_df %>% 
    ggplot(aes(x = bins, y = pop_deg_days, fill = model_level)) +
    geom_col(alpha = 0.5, position = "identity", width = 0.8) +
    scale_fill_manual(values = color_palette) +
    theme_minimal() +
    theme(
      strip.text   = element_text(size = 10, face = "bold"), 
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x  = element_text(size = 10),
      axis.title.x = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      plot.margin = margin(0, 0, 0, 0)
    )
  
  final_plot <- cowplot::plot_grid(response_plot, exposure_plot, 
                                   ncol = 1, rel_heights = c(2, 0.5),
                                   align = 'v', axis = 'lr') 
  final_plot
  return(list(response_plot, exposure_plot, final_plot))
  
}

plot_annual_deaths <- function(data) {
  
  data <- data %>% 
    mutate(
      annual_risk = deaths / annual_total_pop,                           # Annual death risk from heat or cold per person
      annual_share_total_mortality = deaths / annual_total_death         # Percent of total deaths
    ) %>% 
    group_by(side, model_level) %>% 
    summarise(
      deaths = sum(deaths, na.rm = T),
      annual_risk = mean(annual_risk, na.rm =  T),                                   # Mean across years
      annual_share_total_mortality = mean(annual_share_total_mortality, na.rm = T), # Mean across years share
      .groups = "drop"
    ) %>% 
    mutate(annual_risk_1000 = 1 - ((1 - annual_risk) ^ 1000))

  # Process data and generate the plot
  plot <- data %>% 
    ggplot(aes(x = model_level, y = annual_share_total_mortality, fill = side)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(
      values = c("cold" = "#017c88", "hot" = "#f88379"),
      labels = c(
        "Cold" = "cold",
        "Hot" = "hot"
      )
    ) +
    labs(fill = "Side", y = "% of Total Annual Deaths") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 10, face = "bold", margin = margin(r = 10)), 
      axis.title.x = element_blank(),
      legend.position = "top"
    ) +
    theme(
      strip.text   = element_text(size = 10, face = "bold"), 
      axis.title.y = element_blank(),
      # axis.text.y  = element_blank(),
      # axis.ticks.y = element_blank(),
      axis.text.x  = element_text(size = 10),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      panel.grid.major = element_blank()
    )
  return(list(plot = plot, df = data))
}

# -------------------------------
# 1. Plot Pooled Model Results
# -------------------------------
temp_min = -20
temp_max = 40
pooled_plot <- rbind(
  pooled_response_df, 
  pooled_response_df_wlag_30days) %>% 
  select(-model_level) %>% 
  rename(model_level = model_name) %>% 
  mutate(model_level = case_when(
    model_level == "Pooled" ~ "Pooled",
    model_level == "Pooled (Lag = 30)" ~ "Pooled (Lag = 30d)"
  )) %>%
  plot_response(temp_min, temp_max, max_pop)

annual_deaths_pooled_result <- rbind(pooled_deaths_annual_by_side, 
                                     pooled_deaths_annual_by_side_wlag_30days) %>%
  select(-model_level) %>% 
  rename(model_level = model_name) %>% 
  mutate(model_level = case_when(
    model_level == "Pooled" ~ "Pooled",
    model_level == "Pooled (Lag = 30)" ~ "Pooled (Lag = 30d)"
  )) %>% 
  plot_annual_deaths()

annual_deaths_pooled <- annual_deaths_pooled_result$plot

# -------------------------------
# 2. Plot Age-Related Results
# -------------------------------
# Prepare response plot for age groups
age_response_plot <- rbind(young_response_df, 
                           adult_response_df, 
                           elderly_response_df) %>% 
  select(-model_level) %>% 
  rename(model_level = model_name) %>% 
  mutate(model_level = as.factor(model_level)) %>% 
  mutate(model_level = factor(model_level, levels = c("Young", "Adult", "Elderly"))) %>% 
  plot_response(temp_min, temp_max, max_pop)

age_response_plot[[1]] +facet_wrap(~ model_level, scales = "free_y")


# Prepare annual deaths plot for age groups

age_deaths_result <- rbind(young_deaths_annual_by_side, 
                           adult_deaths_annual_by_side, 
                           elderly_deaths_annual_by_side) %>% 
  select(-model_level) %>% 
  rename(model_level = model_name) %>% 
  mutate(model_level = as.factor(model_level)) %>% 
  mutate(model_level = factor(model_level, levels = c("Young", "Adult", "Elderly"))) %>% 
  plot_annual_deaths()

age_deaths_plot <- age_deaths_result$plot
# -------------------------------
# 3. Plot Temperature Results
# -------------------------------
# Filter desired levels for temperature response
temp_response_plot <- temp_response_df %>% 
  filter(model_level %in% c('10', '50', '90')) %>%
  mutate(model_level = case_when(
    model_level == '10' ~ "10th",
    model_level == '50' ~ "50th",
    model_level == '90' ~ "90th",
  )) %>% 
  plot_response(temp_min, temp_max, max_pop)

# Prepare temperature annual deaths plot, parsing model level as numeric
temp_deaths_result <- temp_deaths_annual_by_side %>% 
  mutate(percentile = as.numeric(model_level)) %>% 
  mutate(percentile = ntile(percentile, n = 3)) %>% 
  mutate(model_level = case_when(
    is.na(percentile) ~ model_level,  # Preserve original if no percentile found
    percentile == 1 ~ "Low Temp", 
    percentile == 2 ~ "Medium Temp", 
    percentile == 3 ~ "High Temp"
  )) %>% 
  group_by(year, side, model_level) %>% 
  summarise(across(starts_with("deaths"), sum, na.rm = TRUE),  # Avoid NA sum errors
            annual_total_pop = sum(annual_total_pop, na.rm = TRUE),
            annual_total_death = sum(annual_total_death, na.rm = TRUE),
            .groups = "drop") %>% 
  mutate(model_level = factor(model_level, levels = c("Low Temp", "Medium Temp", "High Temp"))) %>% 
  plot_annual_deaths() 

temp_deaths_plot <- temp_deaths_result$plot 

temp_deaths_annual_by_side %>% 
  mutate(model_level = as.numeric(model_level)) %>% 
  arrange(side, model_level) %>% 
  plot_annual_deaths()
# -------------------------------
# 4. Plot Income Results
# -------------------------------
# Filter desired levels for income response
income_response_plot <- income_response_df %>% 
  filter(model_level %in% c('10', '50', '90')) %>%
  mutate(model_level = case_when(
    model_level == '10' ~ "10th",
    model_level == '50' ~ "50th",
    model_level == '90' ~ "90th",
  ))%>% 
  plot_response(temp_min, temp_max, max_pop)

# Prepare income annual deaths plot, parsing model level as numeric
income_deaths_result <- income_deaths_annual_by_side %>% 
  filter(annual_total_death != 0) %>% 
  mutate(percentile = as.numeric(model_level)) %>% 
  mutate(percentile = ntile(percentile, n = 3)) %>% 
  mutate(model_level = case_when(
    is.na(percentile) ~ model_level,  # Preserve original if no percentile found
    percentile == 1 ~ "Low Income", 
    percentile == 2 ~ "Medium Income", 
    percentile == 3 ~ "High Income"
  )) %>% 
  group_by(year, side, model_level) %>% 
  summarise(across(starts_with("deaths"), sum, na.rm = TRUE),  # Avoid NA sum errors
            annual_total_pop = sum(annual_total_pop, na.rm = TRUE),
            annual_total_death = sum(annual_total_death, na.rm = TRUE),
            .groups = "drop") %>% 
  mutate(model_level = factor(model_level, levels = c("Low Income", "Medium Income", "High Income"))) %>% 
  plot_annual_deaths() 

income_deaths_plot <- income_deaths_result$plot

income_deaths_annual_by_side %>% 
  filter(annual_total_death != 0) %>% 
  mutate(model_level = as.numeric(model_level)) %>% 
  arrange(side, model_level) %>% 
  plot_annual_deaths() 

# -------------------------------
# 5. Plot Decade Results
# -------------------------------
# Prepare decade response plot with recoded model levels
decade_response_plot <- decade_response_df %>% 
  mutate(model_level = case_when(
    model_level == 1970 ~ "1970s",
    model_level == 1980 ~ "1980s",
    model_level == 1990 ~ "1990s",
    model_level == 2000 ~ "2000s",
    model_level == 2010 ~ "2010s",
    model_level == 2020 ~ "2020s",
    TRUE ~ as.character(model_level)
  )) %>% 
  plot_response(temp_min, temp_max, max_pop)

# Prepare decade annual deaths plot with recoded model levels
decade_deaths_result <- decade_deaths_annual_by_side %>% 
  mutate(model_level = case_when(
    model_level == 1970 ~ "1970s",
    model_level == 1980 ~ "1980s",
    model_level == 1990 ~ "1990s",
    model_level == 2000 ~ "2000s",
    model_level == 2010 ~ "2010s",
    model_level == 2020 ~ "2020s",
    TRUE ~ as.character(model_level)
  )) %>% 
  plot_annual_deaths()

decade_deaths_plot <-  decade_deaths_result$plot

# -------------------------------
# Final Layout (example)
# -------------------------------
# (pooled_plot | age_response_plot | temp_response_plot | income_response_plot | decade_response_plot) /
# (annual_deaths_pooled | age_deaths_plot | temp_deaths_plot | income_deaths_plot | decade_deaths_plot)
library(cowplot)
# age_response_plot | temp_response_plot + coord_cartesian(ylim = c(-0.00000338, 0.000113))

# First row: response plots (5 columns)
responses_row_1 <- plot_grid(
  pooled_plot[[3]], 
  age_response_plot[[3]], 
  temp_response_plot[[3]], 
  income_response_plot[[3]], 
  decade_response_plot[[3]], 
  ncol = 5, 
  labels = "auto"
)

# Second row: annual deaths plots (5 columns)
# common_y_scale <- c(0, 0.13)
deaths_row_1 <- plot_grid(
  annual_deaths_pooled, 
  age_deaths_plot, 
  temp_deaths_plot, 
  income_deaths_plot, 
  decade_deaths_plot, 
  ncol = 5, 
  labels = "auto",
  align = "h"
)

# Combine the two rows into a final grid with 2 rows and 5 columns in total
final_plot_1 <- plot_grid(responses_row_1, deaths_row_1, ncol = 1)
final_plot_1
plot_out_dir <- file.path(fig_country_dir, paste0(country_name, "_MortalityTempAnalysis_combined_plot.png"))
ggsave(plot = final_plot_1, filename = plot_out_dir, width = 20, height = 8, units = "in", dpi = 300)
# 
# # List of mortality response plots for different models
# mortality_response_plots <- list(
#   pooled_plot, 
#   age_response_plot, 
#   temp_response_plot, 
#   income_response_plot, 
#   decade_response_plot
# )
# 
# # Calculate shared y-axis limits
# max_y <- max(sapply(mortality_response_plots, function(p) {
#   y_range <- ggplot_build(p)$layout$panel_params[[1]]$y.range
#   if (is.null(y_range)) stop("One or more plots have no y-axis range.")
#   return(y_range[2])
# }))
# min_y <- min(sapply(mortality_response_plots, function(p) {
#   y_range <- ggplot_build(p)$layout$panel_params[[1]]$y.range
#   if (is.null(y_range)) stop("One or more plots have no y-axis range.")
#   return(y_range[1])
# }))
# 
# y_axis_limits <- c(min_y, max_y)
# 
# # Apply shared y-axis limits and theme adjustments
# mortality_response_plots <- lapply(seq_along(mortality_response_plots), function(i) {
#   plot <- mortality_response_plots[[i]] + 
#     ylim(y_axis_limits) +
#     theme_minimal() + 
#     theme(
#       strip.text = element_text(size = 10, face = "bold"), 
#       axis.title.y = if (i == 1) element_text(size = 10, face = "bold", margin = margin(r = 10)) else element_blank(),
#       axis.text.y = if (i == 1) element_text(size = 10) else element_blank(),
#       axis.ticks.y = if (i == 1) element_line() else element_blank(),
#       axis.text.x = element_text(size = 10),
#       axis.title.x = element_blank(),
#       legend.position = c(.20, .90),
#       panel.grid.major = element_blank()
#     )
#   return(plot)
# })
# 
# # Combine plots into a single grid
# response_row_2 <- wrap_plots(mortality_response_plots, ncol = length(mortality_response_plots))
# 
# # Combine the two rows into a final grid with 2 rows and 5 columns in total
# final_plot_2 <- plot_grid(responses_row_2, deaths_row, ncol = 1, align = "vh")
# final_plot_2
# plot_out_dir <- file.path(fig_country_dir, paste0(country_name, "_MortalityTempAnalysis_combined_plot_2.png"))
# ggsave(plot = final_plot_2, filename = plot_out_dir, width = 20, height = 8, units = "in", dpi = 300)
