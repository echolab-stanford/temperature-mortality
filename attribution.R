################################################################################
# ENVIRONMENT SETUP
################################################################################

# Clear environment
rm(list = ls())
gc()

# Load required libraries
pacman::p_load(fastverse, tidyverse, arrow, fixest, matrixStats, ISOweek, here)

# Set configuration
setFixest_nthreads(6, save = TRUE)

################################################################################
# DIRECTORY AND PATH SETUP
################################################################################

# here::here() resolves paths relative to the .Rproj root — no setwd() needed.
# Open the project via temperature-mortality.Rproj to ensure here() works correctly.

# Load country configuration
source(here::here("config", "countries.R"))

# ── USER SETTINGS ─────────────────────────────────────────────────────────────
# country_name: one of the single-dataset keys in config/countries.R
#   e.g. "US", "EU", "MEX", "BRA", "LA", "PHL", "MYS"
country_name <- "BRA"

# dataset_tag: set ONLY for countries with multiple cleaned datasets.
#   NULL  → single-dataset country; cleaned_dir = data/{country}/cleaned/
#   "THL_1997_2025" → data/THL/cleaned/THL_1997_2025/  (total deaths, 913 AMCs)
#   "THL_2005_2025" → data/THL/cleaned/THL_2005_2025/  (age-specific, 920 AMCs)
dataset_tag <- NULL

# si_folder: sensitivity analysis variant to run on the selected dataset.
#   NULL (default/winsor), "log_rate", "same_year", "age_standardized"
si_folder <- "winsor"
# ──────────────────────────────────────────────────────────────────────────────

# Resolve config key: dataset_tag when set, else country_name
config_key <- if (!is.null(dataset_tag)) dataset_tag else country_name

# Validate config key
if (!config_key %in% names(country_config)) {
  stop("'", config_key, "' not found in config/countries.R. ",
       "Available options: ", paste(names(country_config), collapse = ", "))
}

# Pull country-specific config
cfg <- country_config[[config_key]]

# For dataset_tag entries, cfg$country_name holds the actual country code used
# in file names (e.g. "THL"). Override country_name so all downstream paste0()
# calls that build file names remain unchanged.
if (!is.null(cfg$country_name)) country_name <- cfg$country_name

# Set up directory structure
data_dir      <- here::here("data")
processed_dir <- here::here("processed")

if (!is.null(dataset_tag)) {
  # Multi-dataset country: nest both cleaned input and processed output under dataset_tag
  cleaned_dir           <- file.path(data_dir,      country_name, "cleaned",   dataset_tag)
  processed_country_dir <- file.path(processed_dir, country_name, dataset_tag)
} else {
  cleaned_dir           <- file.path(data_dir,      country_name, "cleaned")
  processed_country_dir <- file.path(processed_dir, country_name)
}

# Ensure base processed directory exists before si_folder logic appends to it
if (!dir.exists(processed_country_dir)) {
  dir.create(processed_country_dir, recursive = TRUE, showWarnings = FALSE)
}

# Load helper functions and data
# source(paste(script_dir, "helper_functions.R", sep = "/"))
# mort <- read_parquet( paste(cleaned_dir, paste0(country_name, "FullPanel_complete.pq"), sep = "/") )
# temp_bins <- read_parquet( paste(cleaned_dir, paste0(country_name, "TempBins_monthly_complete.pq"), sep = "/") )
# daily_temp <- read_parquet( paste(cleaned_dir, paste0(country_name, "Daily_temp_covariates.pq"), sep = "/") )

# mort <- read_parquet( paste(
#   "/Users/khusel/BurkeLab Dropbox/projects/temperature-mortality-descriptive/data/", country_name, "/cleaned", 
#   paste0(country_name, "FullPanel_complete.pq"), sep = "/") )
# temp_bins <- read_parquet( paste(
#   "/Users/khusel/BurkeLab Dropbox/projects/temperature-mortality-descriptive/data/", country_name, "/cleaned", 
#   paste0(country_name, "TempBins_monthly_complete.pq"), sep = "/") )
# daily_temp_cols <- if (length(cfg$models_to_run) == 1 && cfg$models_to_run == "Pooled") {
#   c("adm_id", "year", "month", "day", "temp1", "temp2", "temp3", "temp4")
# } else {
#   NULL  # load all columns for heterogeneity models
# }
# daily_temp <- read_parquet(
#   paste("/Users/khusel/BurkeLab Dropbox/projects/temperature-mortality-descriptive/data/",
#         country_name, "/cleaned",
#         paste0(country_name, "Daily_temp_covariates.pq"), sep = "/"),
#   col_select = if (!is.null(daily_temp_cols)) all_of(daily_temp_cols) else everything()
# )

mort <- read_parquet(file.path(cleaned_dir, paste0(country_name, "FullPanel_complete.pq")))
temp_bins <- read_parquet(file.path(cleaned_dir, paste0(country_name, "TempBins_monthly_complete.pq")))
daily_temp_cols <- if (length(cfg$models_to_run) == 1 && cfg$models_to_run == "Pooled") {
  c("adm_id", "year", "month", "day", "temp1", "temp2", "temp3", "temp4")
} else {
  NULL  # load all columns for heterogeneity models
}
daily_temp <- read_parquet(
  file.path(cleaned_dir, paste0(country_name, "Daily_temp_covariates.pq")),
  col_select = if (!is.null(daily_temp_cols)) all_of(daily_temp_cols) else everything()
)

setDT(mort)
setDT(temp_bins)
setDT(daily_temp)

################################################################################
# COUNTRY-SPECIFIC DATA PREPARATION
################################################################################

# Pull variables from config
time_agg_var       <- cfg$time_agg_var
right_hand_var_list <- cfg$rhs_vars
right_hand_vars    <- paste(right_hand_var_list, collapse = " + ")
fixed_effects      <- cfg$fixed_effects

# # Year subset (from config; NULL = no filter)
# if (!is.null(cfg$year_range)) {
#   message("Subsetting to years ", cfg$year_range[1], "–", cfg$year_range[2])
#   mort <- mort %>% filter(year >= cfg$year_range[1] & year <= cfg$year_range[2])
#   daily_temp <- daily_temp %>% filter(year >= cfg$year_range[1] & year <= cfg$year_range[2])
# }

# Set data.table keys based on temporal resolution
if (time_agg_var == "week") {

  setkeyv(mort,       c("adm_id", "year", "week"))
  setkeyv(temp_bins,  c("adm_id", "year", "week"))
  setkeyv(daily_temp, c("adm_id", "year", "month", "day"))

  # Create month variable (needed for some downstream steps)
  mort[, date  := ISOweek2date(sprintf("%d-W%02d-1", year, week))]
  mort[, month := month(date)]

  # Generate weekly lag variables (EU uses lags t1–t4)
  new_lag_names            <- c("temp_d1", "temp_d2", "temp_d3", "temp_d4",
                                "precip_d1", "precip_d2")
  weather_polynomial_names <- c("temp_d1_l0", "temp_d2_l0", "temp_d3_l0", "temp_d4_l0",
                                "precip_d1_l0", "precip_d2_l0")
  for (lag_value in 1:4) {
    lag_cols <- paste0(new_lag_names, "_t", lag_value)
    mort[, (lag_cols) := lapply(.SD, shift, n = lag_value),
         by = adm_id, .SDcols = weather_polynomial_names]
  }

} else {

  setkeyv(mort,       c("adm_id", "year", "month"))
  setkeyv(temp_bins,  c("adm_id", "year", "month"))
  setkeyv(daily_temp, c("adm_id", "year", "month", "day"))

}

################################################################################
# ANALYSIS SETUP
################################################################################

if (si_folder == "log_rate") {
  
  # Set up log rate analysis
  processed_country_dir <- file.path(processed_country_dir, si_folder)
  
  if (!dir.exists(processed_country_dir)) {
    dir.create(processed_country_dir, recursive = FALSE)
  }
  
  # mort <- mort %>% mutate(
  #   rate_all = log(rate_all), 
  #   rate_young = log(rate_young), 
  #   rate_adult = log(rate_adult), 
  #   rate_elderly = log(rate_elderly)
  #   ) # no longer need, we use poisson regression
  
  log_rate <- TRUE
  
} else if (si_folder == "same_year") {
  
  # Set up same year analysis
  processed_country_dir <- file.path(processed_country_dir, si_folder)
  
  if (!dir.exists(processed_country_dir)) {
    dir.create(processed_country_dir, recursive = FALSE)
  }
  
  mort <- mort %>% filter(year > 1999 & year < 2020)
  log_rate <- F
  
} else if (si_folder == "age_standardized") {

  # Age-standardized rate (full sample)
  processed_country_dir <- file.path(processed_country_dir, si_folder)
  if (!dir.exists(processed_country_dir)) dir.create(processed_country_dir, recursive = FALSE)
  if (!"age_std_rate" %in% names(mort)) stop("age_std_rate not found in data for country: ", country_name)
  mort <- mort %>%
    mutate(rate_all = age_std_rate) %>%
    filter(!is.na(age_std_rate), is.finite(age_std_rate))
  log_rate <- FALSE

} else if (si_folder == "age_std_no_agg") {

  # Age-standardized rate; exclude ALL AMCs that have any has_agg_only_deaths == TRUE
  processed_country_dir <- file.path(processed_country_dir, si_folder)
  if (!dir.exists(processed_country_dir)) dir.create(processed_country_dir, recursive = FALSE)
  if ("has_agg_only_deaths" %in% names(mort)) {
    agg_only_amcs <- mort[has_agg_only_deaths == TRUE, unique(adm_id)]
    message("Excluding ", length(agg_only_amcs), " AMCs (has_agg_only_deaths == TRUE); ",
            uniqueN(mort$adm_id) - length(agg_only_amcs), " AMCs remain.")
    mort       <- mort[!adm_id %in% agg_only_amcs]
    daily_temp <- daily_temp[!adm_id %in% agg_only_amcs]
    temp_bins  <- temp_bins[!adm_id %in% agg_only_amcs]
  } else {
    message("has_agg_only_deaths not found in data; no AMC exclusion applied.")
  }
  if (!"age_std_rate" %in% names(mort)) stop("age_std_rate not found in data for country: ", country_name)
  mort <- mort %>%
    mutate(rate_all = age_std_rate) %>%
    filter(!is.na(age_std_rate), is.finite(age_std_rate))
  log_rate <- FALSE

} else if (si_folder == "age_std_2014_2025") {

  # Age-standardized rate; 2014-2025 year subset
  processed_country_dir <- file.path(processed_country_dir, si_folder)
  if (!dir.exists(processed_country_dir)) dir.create(processed_country_dir, recursive = FALSE)
  mort       <- mort[year >= 2014]
  daily_temp <- daily_temp[year >= 2014]
  temp_bins  <- temp_bins[year >= 2014]
  if (!"age_std_rate" %in% names(mort)) stop("age_std_rate not found in data for country: ", country_name)
  mort <- mort %>%
    mutate(rate_all = age_std_rate) %>%
    filter(!is.na(age_std_rate), is.finite(age_std_rate))
  log_rate <- FALSE

} else if (si_folder == "age_std_2005_2012") {

  # Age-standardized rate; 2005-2012 year subset
  processed_country_dir <- file.path(processed_country_dir, si_folder)
  if (!dir.exists(processed_country_dir)) dir.create(processed_country_dir, recursive = FALSE)
  mort       <- mort[year <= 2012]
  daily_temp <- daily_temp[year <= 2012]
  temp_bins  <- temp_bins[year <= 2012]
  if (!"age_std_rate" %in% names(mort)) stop("age_std_rate not found in data for country: ", country_name)
  mort <- mort %>%
    mutate(rate_all = age_std_rate) %>%
    filter(!is.na(age_std_rate), is.finite(age_std_rate))
  log_rate <- FALSE

} else if (si_folder == "crude_str_fe") {

  # Crude mortality rate; FEs: adm_id^year_str_break + adm_id^pre_post_break^month
  # year_str_break splits 2013 into 20130 / 20131 (structural break).
  # pre_post_break: "pre" = through 20130, "post" = 20131 onwards.
  # Both variables are now in the FullPanel for both THL_1997_2025 and THL_2005_2025.
  processed_country_dir <- file.path(processed_country_dir, si_folder)
  if (!dir.exists(processed_country_dir)) dir.create(processed_country_dir, recursive = FALSE)

  mort <- mort %>%
    mutate(across(starts_with("rate"), ~ fifelse(.x > collapse::fnth(.x, 0.99), collapse::fnth(.x, 0.99), .x))) %>%
    mutate(mean_income = fifelse(mean_income > collapse::fnth(mean_income, 0.99), collapse::fnth(mean_income, 0.99), mean_income))

  fixed_effects <- "adm_id^year_str_break + adm_id^pre_post_break^month"
  log_rate <- FALSE

} else if (si_folder == "age_std_2015_str_fe") {

  # Age-standardized rate, full 2005-2025 sample; FEs: adm_id^year_str_break + adm_id^pre_post_break^month
  processed_country_dir <- file.path(processed_country_dir, si_folder)
  if (!dir.exists(processed_country_dir)) dir.create(processed_country_dir, recursive = FALSE)

  if (!"age_std_rate" %in% names(mort)) stop("age_std_rate not found in data for country: ", country_name)
  mort <- mort %>%
    mutate(rate_all = age_std_rate) %>%
    filter(!is.na(age_std_rate), is.finite(age_std_rate))

  fixed_effects <- "adm_id^year_str_break + adm_id^pre_post_break^month"
  log_rate <- FALSE

} else if (si_folder == "newdata_agerate_fixed" | si_folder == "newdata_original" | si_folder == "newdata_1990_2023") {
  
  # Set up age standardized analysis
  processed_country_dir <- file.path(processed_country_dir, si_folder)
  
  if (!dir.exists(processed_country_dir)) {
    dir.create(processed_country_dir, recursive = FALSE)
  }
  
  # mort <- mort %>% mutate(rate_all = age_std_rate) %>% 
  log_rate <- F
  
} else {
  # Default: winsorizing rates 
  mort <- mort %>% mutate(across(starts_with("rate"), ~ fifelse(.x >  collapse::fnth(.x, 0.99),  collapse::fnth(.x, 0.99), .x)))
  # winsor income
  mort <- mort %>% mutate(mean_income = fifelse(mean_income>  collapse::fnth(mean_income, 0.99),  collapse::fnth(mean_income, 0.99), mean_income))

  log_rate <- F
}

################################################################################
# MODEL FORMULA DEFINITIONS
################################################################################

# Create model formulas
pooled_model_formula <- as.formula(paste("rate_all ~", right_hand_vars, "|", fixed_effects))
young_model_formula <- as.formula(paste("rate_young ~", right_hand_vars, "|", fixed_effects))
adult_model_formula <- as.formula(paste("rate_adult ~", right_hand_vars, "|", fixed_effects))
elderly_model_formula <- as.formula(paste("rate_elderly ~", right_hand_vars, "|", fixed_effects))
temp_model_formula <- as.formula(paste("rate_all ~", right_hand_vars, " + ", paste0("(", right_hand_vars, "):mean_tm"),  "|", fixed_effects))
income_model_formula <- as.formula(paste("rate_all ~", right_hand_vars, " + ", paste0("(", right_hand_vars, "):mean_income"),  "|", fixed_effects))
baseline_model_formula <- as.formula(paste("rate_all ~", right_hand_vars, " + ", paste0("(", right_hand_vars, "):baseline_mean"),  "|", fixed_effects))

decade_right_hand_vars <- unlist(lapply(right_hand_var_list, function(v) {paste0("i(decade,", v, ")")}))
decade_model_formula <- as.formula(paste("rate_all ~", paste(decade_right_hand_vars, collapse = " + "), "|", fixed_effects))

# Create model lists
model_formula_list <- list(
  "Pooled" = pooled_model_formula,
  "Young" = young_model_formula,
  "Adult" = adult_model_formula,
  "Elderly" = elderly_model_formula,
  "temp" = temp_model_formula,
  "income" = income_model_formula, 
  "decade" = decade_model_formula,
  "baseline" = baseline_model_formula
)

################################################################################
# MAIN ANALYSIS LOOP
################################################################################

# Filter model list to only models configured for this country
model_formula_list <- model_formula_list[names(model_formula_list) %in% cfg$models_to_run]

# Further restrict for same_year sensitivity (pooled only)
if (si_folder == "same_year") {
  model_formula_list <- model_formula_list["Pooled"]
}

# Main loop through models
for (model_name in names(model_formula_list)) {

  # model_name <- "Pooled"
  message("Processing model: ", model_name)
  
  # Set analysis parameters
  Nboot <- 100
  formula <- model_formula_list[[model_name]]
  weights <- NULL
  override <- F
  
  ############################################################################
  # SECTION 1: BOOTSTRAPPING
  ############################################################################
  
  out_dir <- file.path(processed_country_dir, paste0(country_name, "_BootstrapResult_", tolower(model_name), "_n", Nboot, ".rds"))
  if (!dir.exists(dirname(out_dir))) {
  dir.create(dirname(out_dir), recursive = TRUE)
  }

  if (!override & file.exists(out_dir)){
    
    message("Loading existing bootstrap results from: ", out_dir)
    bstrap <- readRDS(out_dir)
    
  } else {
    
    options(warn = -1)
    set.seed(12)
    
    setkey(mort, adm_id)
    unique_ids <- unique(mort$adm_id)
    n_units <- length(unique_ids)
    
    results <- lapply(1:Nboot, function(i) {
      message("Bootstrapping: ", i)
      
      # Resample spatial units with replacement
      sampled_units <- unique_ids[sample.int(n_units, n_units, replace = TRUE)]
      dtboot <- mort[.(sampled_units), allow.cartesian = TRUE]
      
      # dtboot <- mort[adm_id %in% sampled_units]
      
      if (is.null(weights)) {
        if (log_rate == TRUE) {
          message("Running quasipoisson model ...")
          fit <- fixest::feglm(formula, family = "quasipoisson", data = dtboot, nthreads = 6)
        } else {
          fit <- fixest::feols(formula, data = dtboot, nthreads = 6)
        }
 
      } else {
        if (log_rate == TRUE) {
          message("Running quasipoisson model ...")
          fit <- fixest::feglm(formula, family = "quasipoisson" , data = dtboot, weights = dtboot[[weights]], nthreads = 6)
        } else {
          fit <- fixest::feols(formula, data = dtboot, weights = dtboot[[weights]], nthreads = 6)
        }
      }
      
      model_results <- coef(fit)
      model_results <- as.data.table(t(model_results), stringsAsFactors = FALSE)
      model_results$iteration <- i
      return(model_results)
      
    })
    
    bstrap <- do.call(rbind, results)
    
    out_dir <- file.path(processed_country_dir, paste0(country_name, "_BootstrapResult_", tolower(model_name), "_n", Nboot, ".rds"))
    
    message("Saving bootstrap results to: ", out_dir)
    write_rds(bstrap, out_dir)
    options(warn = 0)
  }
  
  ############################################################################
  # SECTION 2: MODEL ESTIMATION
  ############################################################################
  setkey(mort, adm_id)
  setDT(bstrap)
  if (log_rate == TRUE) {
    coefs <- feglm(formula, family = "quasipoisson", data = mort, weights = NULL, only.coef = TRUE)
  } else {
    coefs <- feols(formula, data = mort, weights = NULL, only.coef = TRUE)
  }
  
  # Set analysis parameters
  # pwd <- mapping_pooled_pwd[[model_name]]
  coef_pattern <- "temp"
  pred_temp_range <- -50:50
  pred_temps_poly <- poly(pred_temp_range, degree = 4, raw = TRUE)
  
  ############################################################################
  # SECTION 3: MODEL-SPECIFIC ANALYSIS
  ############################################################################

  if (model_name == "Pooled" | model_name == "Pooled (lag)" | model_name == "Young" | model_name == "Adult" | model_name == "Elderly") {
    
    ########################################################################
    # POOLED ANALYSIS SETUP
    ########################################################################
    
    # 0. Population-Specific Variables
    # Assign population and death columns based on model_name
    if (model_name == "Young") {
      pop_var <- "pop_young"
      death_var <- "n_deaths_young"
      outcome_var <- "rate_young"
    } else if (model_name == "Adult") {
      pop_var <- "pop_adult"
      death_var <- "n_deaths_adult"
      outcome_var <- "rate_adult"
    } else if (model_name == "Elderly") {
      pop_var <- "pop_elderly"
      death_var <- "n_deaths_elderly"
      outcome_var <- "rate_elderly"
    } else {
      pop_var <- "pop_all"              # Default: All populations
      death_var <- "n_deaths_all"
      outcome_var <- "rate_all"
    }
    
    if (country_name == "EU") { 
      time_agg_var <- "week"
    } else {
      time_agg_var <- "month"
    }
    
    # 1. Model Coefficient Processing
    # Transpose coefficient matrix for downstream operations
    coefs_before_sum <- t(as.matrix(coefs))  # coefs: Original regression coefficients
    
    # 2. Cumulative Lag Effects Calculation
    # Sum coefficients for each polynomial term (temp_d1 to temp_d4), excluding precipitation terms
    poly_groups <- c("temp_d1", "temp_d2", "temp_d3", "temp_d4")
    coefs_sum <- sapply(poly_groups, function(poly) {
      # Identify columns matching the polynomial group (e.g., temp_d1_lag1, temp_d1_lag2)
      poly_names = colnames(coefs_before_sum)[grepl(poly, colnames(coefs_before_sum)) & !grepl("prec", colnames(coefs_before_sum))]
      rowSums(coefs_before_sum[, poly_names, drop = FALSE])  # Sum across lags
    })
    
    coefs_sum <- t(as.matrix(coefs_sum))  # Transpose back to original orientation
    
    # 3. Temperature Response
    # Predict mortality response using polynomial-transformed temperatures
    pred_vector <- pred_temps_poly %*% t(coefs_sum)  # pred_temps_poly: Temp values transformed via polynomial basis
    
    # 4. Temperature Range Trimming
    # Restrict predictions to 1%-99% percentile of observed temperatures to exclude extremes
    observed_temps <- daily_temp$temp1  # Observed temperature data
    t1_index <- collapse::fquantile(observed_temps, 0.01, na.rm = TRUE)  # 1st percentile
    t99_index <- collapse::fquantile(observed_temps, 0.998, na.rm = TRUE)  # 99th percentile
    mmt_options <- which(pred_temp_range > t1_index & pred_temp_range < t99_index)  # Indices within range
    new_temp_range <- pred_temp_range[mmt_options]  # Subset temperature range
    
    # 5. Minimum Mortality Temperature (MMT)
    # Find temperature with lowest predicted mortality (MMT)
    mmt <- new_temp_range[which.min(pred_vector[mmt_options])]
    mmt_matrix <- poly(mmt, degree = 4, raw = T)  # Create polynomial terms for MMT
    pred_temps_scaled <- scale(pred_temps_poly, center = as.numeric(mmt_matrix), scale = FALSE)  # Center temps at MMT
    
    # Predict mortality response using centered temperatures
    pred_vector <- pred_temps_scaled %*% t(coefs_sum)
    
    # 6. Bootstrap Coefficient Processing
    # Extract temperature-related coefficients from bootstrap results
    temp_cols <- names(bstrap)[grep("temp", names(bstrap))]  # Columns with "temp" in name
    bst_coefs_before_sum <- bstrap[, ..temp_cols]  # Subset bootstrap coefficients
    bst_coefs_before_sum <- as.matrix(bst_coefs_before_sum)  # Convert to matrix
    
    # 7. Bootstrap Response Calculation
    # Sum bootstrap coefficients by polynomial group (same as step 2 but for bootstrapped coefs)
    bst_coefs_sum <- sapply(poly_groups, function(poly) {
      poly_names = colnames(bst_coefs_before_sum)[grepl(poly, colnames(bst_coefs_before_sum)) & !grepl("prec", colnames(bst_coefs_before_sum))]
      rowSums(bst_coefs_before_sum[, poly_names, drop = FALSE])
    })
    
    # Predict mortality response for each bootstrap iteration
    bstrap_result <- pred_temps_scaled %*% t(bst_coefs_sum)  # Rows: temps, cols: bootstrap samples
    
    # Calculate response statistics (2.5%, median, 97.5%, mean)
    response_stats <- cbind(
      matrixStats::rowQuantiles(bstrap_result, probs = c(0.025, 0.5, 0.975), na.rm = TRUE),
      mean = matrixStats::rowMeans2(bstrap_result, na.rm = TRUE)
    )
    
    # 8. Response Result 
    response <- data.table(
      pred_q025 = response_stats[, 1],  # 2.5% quantile
      pred = response_stats[, 4],       # Mean response
      pred_q975 = response_stats[, 3],  # 97.5% quantile
      bins = pred_temp_range            # Temperature bins
    )[, `:=` (
      key = "pooled",                   # Identifier for pooled analysis
      model = model_name,               # Population group (e.g., "Young")
      t1 = t1_index,                    # 1st percentile temp
      t99 = t99_index                   # 99th percentile temp
    )]
    
    # 9. Add annual average person days exposure
    midpoints <- (pred_temp_range[-1] + pred_temp_range[-length(pred_temp_range)]) / 2
    midpoints <- as.character(midpoints)
    
    cols <- c("adm_id", "year", time_agg_var, pop_var)
    
    # Population and degree days bin exposures
    pwd <- merge(mort[, ..cols], temp_bins, by = c("adm_id", "year", time_agg_var))
    setnames(pwd, pop_var, "pop")
    
    # Annual population person-day exposure
    annual_pwd_bins <- pwd[!is.na(pop)]
    rm(pwd)
    gc(full = TRUE)
    
    annual_pwd_bins[, (midpoints) := .SD * pop, .SDcols = midpoints]
    annual_pwd_bins <- annual_pwd_bins[, lapply(.SD, sum), .SDcols = midpoints, by = .(year)]
    annual_pwd_bins <- melt(annual_pwd_bins, measure.vars = midpoints, variable.name = "bins", value.name = "pop_deg_days", na.rm = TRUE)
    
    # Get overall deaths and merge with annual person-day exposure
    annual_pwd_bins <- merge(annual_pwd_bins, mort[, .(total_deaths = sum(get(death_var), na.rm = TRUE)), by = .(year)], by = "year")
    annual_pwd_bins[, bins := ceiling(as.numeric(as.character(bins)))]
    
    # Take average annual person-day exposure and merge with response
    response <- merge(response, annual_pwd_bins[, .(pop_deg_days = mean(pop_deg_days)), by = .(bins)], by = "bins", all.x = TRUE) 
    
    #   10. Bin-specific deaths for Fig 3 (temperature attribution histogram)
    #
    # AGGREGATION METHOD: mean-of-ratios (Method A)
    #
    #   We compute the annual attributable fraction for each temperature bin:
    #     share_{b,t} = bins_deaths_{b,t} / total_deaths_t
    #   and then average across years:
    #     share_b = (1/T) * sum_t( bins_deaths_{b,t} / total_deaths_t )
    #
    #   Interpretation: "In a typical year, X% of [age-group] deaths are
    #   attributable to temperature bin b."
    #
    #   NOTE: total_deaths uses the age-specific death variable (death_var), so
    #   the fraction is relative to the same-age denominator — e.g. adult deaths
    #   attributable / total adult deaths — not all-cause mortality.
    bstrap_long <- melt(data.table(bins = -50:50, bstrap_result), id.vars = "bins", variable.name = "iteration", value.name = "pred")
    bstrap_year_bins_full <- bstrap_long[annual_pwd_bins, on = "bins", allow.cartesian = TRUE][, bins_deaths := pred * pop_deg_days]

    # deaths_by_year_bins: average bootstrap iterations but KEEP year structure.
    # The per-year (bins_deaths, total_deaths) pairs are passed to plot_deaths_at_bins(),
    # which computes the annual ratio and then averages — implementing Method A above.
    deaths_by_year_bins <- bstrap_year_bins_full[, .(bins_deaths = mean(bins_deaths, na.rm = TRUE), total_deaths = mean(total_deaths, na.rm = TRUE)), by = .(year, bins)]
    deaths_by_year_bins[, `:=` (key = "pooled", model = model_name)]

    # response_by_year_bins: collapse year AND iteration — used for inline diagnostic plot only.
    response_by_year_bins <- bstrap_year_bins_full[, .(bins_deaths = mean(bins_deaths), total_deaths = mean(total_deaths)), by = .(bins)]
    response_by_year_bins[, `:=` (key = "pooled", model = model_name, t1 = t1_index, t99 = t99_index)]

    rm(annual_pwd_bins, bstrap_long, bstrap_year_bins_full)
    gc(full = TRUE)
    
    # 11. Hot/Cold Day Classification
    # Classify in-place (no copy) to avoid duplicating the 53M-row daily_temp.
    daily_temp[, side := fifelse(temp1 > mmt, "hot", "cold")]

    # 12. Prepare aggregated panel — aggregate immediately, then drop the side column.
    narrow_panel <- daily_temp[!is.na(side), .(
      temp1 = sum(temp1),
      temp2 = sum(temp2),
      temp3 = sum(temp3),
      temp4 = sum(temp4),
      ndays = .N
    ), by = c("adm_id", "year", time_agg_var, "side")]
    daily_temp[, side := NULL]
    
    #   12. Merge Population/Death Data  
    narrow_panel <- merge(narrow_panel, mort %>% select(adm_id, year, all_of(time_agg_var), pop = all_of(pop_var), n_deaths = all_of(death_var)), by = c("adm_id", "year", time_agg_var))
    
    #   13. Death Prediction via Bootstrap  
    temp_cols <- paste0("temp", 1:4)
    pred_mat <- as.matrix(narrow_panel[, ..temp_cols])
    pred_mat <- pred_mat - matrix(narrow_panel$ndays %o% t(mmt_matrix), ncol = 4)  # Scale by days
    bstrap_deaths <- pred_mat %*% t(bst_coefs_sum)
    
    #   14. Number of estimated deaths   
    if (log_rate == TRUE) {
      # transform back to level
      typical_rate <- mean(mort[[outcome_var]], na.rm = TRUE)
      # bstrap_deaths <- sweep(expm1(bstrap_deaths), 1, narrow_panel$pop, "*") * typical_rate
      bstrap_deaths <- expm1(bstrap_deaths)*narrow_panel$pop*typical_rate
    } else (
      bstrap_deaths <- sweep(bstrap_deaths, 1, narrow_panel$pop, "*")
    )
    
    # aggregate annual deaths
    annual_total_deaths <- mort[, .(total_deaths = sum(get(death_var), na.rm = TRUE)), by = .(year)]

    # deaths by year and side: for each (year, side) group, colSums the bootstrap matrix
    # then summarise across bootstrap iterations — identical to the old wide-DT approach
    groups <- narrow_panel[, .(row_idx = list(.I)), by = .(year, side)]
    deaths_by_year_side <- rbindlist(lapply(seq_len(nrow(groups)), function(i) {
      boot_sums <- colSums(bstrap_deaths[groups$row_idx[[i]], , drop = FALSE], na.rm = TRUE)
      data.table(
        year        = groups$year[i],
        side        = groups$side[i],
        deaths      = mean(boot_sums, na.rm = TRUE),
        deaths_q025 = quantile(boot_sums, 0.025, na.rm = TRUE, names = FALSE),
        deaths_q975 = quantile(boot_sums, 0.975, na.rm = TRUE, names = FALSE)
      )
    }))
    rm(bstrap_deaths, groups)
    deaths_by_year_side[, `:=`(key = "pooled", model = model_name)]
    deaths_by_year_side <- merge(deaths_by_year_side, annual_total_deaths, by = "year", all.x = TRUE)
    deaths_by_year_side <- deaths_by_year_side[, .(year, side, deaths, deaths_q025, deaths_q975, key, model, total_deaths)]

    # deaths_by_year_side <- narrow_panel %>%
    #   group_by(year, side) %>%
    #   summarise(across(starts_with("b"), ~sum(.x, na.rm = T))) %>%
    #   ungroup()  %>% group_by(year, side) %>%
    #   summarise(
    #     deaths = mean(c_across(starts_with("b")), na.rm = TRUE),
    #     deaths_q025 = quantile(c_across(starts_with("b")), 0.025, na.rm = TRUE),
    #     deaths_q975 = quantile(c_across(starts_with("b")), 0.975, na.rm = TRUE),
    #     .groups = "drop"
    #   ) %>%
    #   mutate(key = "pooled", model = model_name) %>% left_join(annual_total_deaths) %>% as.data.table()
    
    # deaths by side 
    deaths <- deaths_by_year_side[, .(                  # Aggregate by year and hot/cold
      total_deaths = mean(total_deaths, na.rm = TRUE),  # Average across years
      deaths_q025 = mean(deaths_q025, na.rm = TRUE),
      deaths = mean(deaths, na.rm = TRUE),
      deaths_q975 = mean(deaths_q975, na.rm = TRUE)
    ), by = .(side)][, key := "pooled"][, model := model_name] 
    
    #   15. Final Output  
    out <- list(
      response_by_year_bins = response_by_year_bins, # Temperature-mort response by year and bins
      deaths_by_year_bins   = deaths_by_year_bins,   # Per-year per-bin deaths (for Fig 3)
      response = response,  # Temperature-mortality curve
      deaths = deaths,       # Hot/cold death estimates
      deaths_by_year_side = deaths_by_year_side # Annual hot/cold death estimates
    )
    
  } else if (model_name %in% c("temp", "income", "baseline")) {
    
    ########################################################################
    # HETEROGENEITY ANALYSIS SETUP
    ########################################################################
    # model_name <- 'temp'
    if (model_name == 'temp') {
      interaction_var <- "mean_tm"
      prc_var <- "temp_percentile"
      pop_var <- "pop_all"              # Default: All populations
      death_var <- "n_deaths_all"
      outcome_var <- "rate_all"
    } else if (model_name == "income") {
      interaction_var <- "mean_income"
      prc_var <- "income_percentile"
      pop_var <- "pop_all"              # Default: All populations
      death_var <- "n_deaths_all"
      outcome_var <- "rate_all"
    } else if (model_name == "baseline") {
      interaction_var <- "baseline_mean"
      prc_var <- "baseline_percentile"
      pop_var <- "pop_all"              # Default: All populations
      death_var <- "n_deaths_all"
      outcome_var <- "rate_all"
    }
    
    # run model one time, and predict temp x each percentile evaluated 
    # find percentile values at 0-100
    quantiles_vec <- collapse::fquantile(mort[[interaction_var]], probs = seq(0, 1, by = 0.01), na.rm = TRUE)
    quantiles_vec <- quantiles_vec[!duplicated(names(quantiles_vec))] 
    
    # Get column names
    cols <- names(coefs)

    # Extract values
    interacted_coefs <- coefs[cols[grepl(interaction_var, cols) & grepl("temp", cols)]]
    intercept_coefs <- coefs[cols[grepl("temp", cols) & !grepl(interaction_var, cols)]]
    
    # outer matrix multification
    quantiles_treated <- quantiles_vec %o% interacted_coefs
    coefs_before_sum <- sweep(quantiles_treated, 2, intercept_coefs, "+") 

    # cumulative effect of lags amd keep matrix format
    poly_groups <- c("temp_d1", "temp_d2", "temp_d3", "temp_d4")
    coefs_sum <- sapply(poly_groups, function(poly) {
      poly_names = colnames(coefs_before_sum)[grepl(poly, colnames(coefs_before_sum)) & !grepl("prec", colnames(coefs_before_sum))]
      rowSums(coefs_before_sum[, poly_names, drop = FALSE])
    })
    
    # for each percentile
    prc_pred_vector <- pred_temps_poly %*% t(coefs_sum)
    
    # Prepare bootstrapped coefficients
    cols <- names(bstrap)
    
    # Find matching columns with better error handling
    eval_cols <- cols[grepl(interaction_var, cols) & grepl("temp", cols)]
    int_cols <- cols[grepl("temp", cols) & !grepl(interaction_var, cols)]
    
    # Convert to matrixec once
    nth_evaluation <- as.matrix(bstrap[, ..eval_cols])
    intercepts <- as.matrix(bstrap[, ..int_cols])
    
    # Process each quantile
    # quantiles_to_loop <- names(quantiles_vec)[names(quantiles_vec) %in% c("10%", "50%", "90%")]
    quantiles_to_loop <- names(quantiles_vec)
    final <- lapply(quantiles_to_loop, function(current_quantile) {
        
      # extract current quantile value, mmt, 1% and 99% percentile observed temp
      # current_quantile <- "7%"
      print(current_quantile)
      quantile_num <- as.numeric(sub("%", "", current_quantile))
      current_quantile_value <- quantiles_vec[[current_quantile]]
      
      current_mort <- mort[as.integer(get(prc_var)) == quantile_num]
      
      if (country_name == "EU") {
        time_agg_var <- "week"
        adm_year_pop_death <- current_mort[, .(adm_id, year, week, pop = get(pop_var), n_deaths = get(death_var))]
      } else {
        time_agg_var <- "month"
        adm_year_pop_death <- current_mort[, .(adm_id, year, month, pop = get(pop_var), n_deaths = get(death_var))]
      }
      
      current_narrow_panel <- daily_temp[as.integer(get(prc_var)) == quantile_num]
      
      current_pred <- as.data.table(prc_pred_vector)[[current_quantile]]
      current_observed_temp <- daily_temp[as.integer(get(prc_var)) == quantile_num]$temp1
      
      t1_index = collapse::fquantile(current_observed_temp, 0.01, na.rm = TRUE)
      t99_index = collapse::fquantile(current_observed_temp, 0.998, na.rm = TRUE)
      
      prc_mmt_options <- which(pred_temp_range > t1_index & pred_temp_range < t99_index)
      prc_new_temp_range <- pred_temp_range[prc_mmt_options]
      mmt <- prc_new_temp_range[which.min(current_pred[prc_mmt_options])]
      
      # Guard: skip quantiles with no valid temp range (empty mmt)
      if (length(mmt) == 0 || length(current_observed_temp) == 0) return(NULL)
      
      mmt_matrix <- poly(mmt, degree = 4, raw = T)
      bst_pred_temps_scaled <- scale(pred_temps_poly, center = as.numeric(mmt_matrix), scale = FALSE)
      
      # Calculate coefficients
      bst_quantile_treated <- nth_evaluation*current_quantile_value
      bst_coefs_before_sum <- bst_quantile_treated + intercepts
      
      # Sum coefficients by poly group (vectorized)
      poly_groups <- c("temp_d1", "temp_d2", "temp_d3", "temp_d4")
      bst_coefs_sum <- sapply(poly_groups, function(poly) {
        poly_names = colnames(bst_coefs_before_sum)[grepl(poly, colnames(bst_coefs_before_sum)) & !grepl("prec", colnames(bst_coefs_before_sum))]
        rowSums(bst_coefs_before_sum[, poly_names, drop = FALSE])
      })
      
      # Predict dose response
      bstrap_result <- bst_pred_temps_scaled %*% t(bst_coefs_sum)
      
      # Calculate response statistics (vectorized)
      response_stats <- cbind(
        matrixStats::rowQuantiles(bstrap_result, probs = c(0.025, 0.5, 0.975), na.rm = TRUE),
        mean = matrixStats::rowMeans2(bstrap_result, na.rm = TRUE)
      )
      
      response <- data.table(
        pred_q025 = response_stats[, 1],
        pred = response_stats[, 4],  # Using mean instead of separate median
        pred_q975 = response_stats[, 3],
        bins = pred_temp_range
      )[, `:=` (
        key = current_quantile,
        model = model_name,
        t1 = t1_index,
        t99 = t99_index
      )]
      
      #   9. Add annual average person days exposure
      midpoints <- (pred_temp_range[-1] + pred_temp_range[-length(pred_temp_range)]) / 2
      midpoints <- as.character(midpoints)
      
      cols <- c("adm_id", "year", time_agg_var, pop_var)
      pwd <- merge( current_mort[, ..cols], temp_bins, by = c("adm_id", "year", time_agg_var) )
      setnames(pwd, pop_var, "pop")
      pwd <- pwd[!is.na(pop)]
      
      pwd[, (midpoints) := lapply(.SD, function(x) x * pop), .SDcols = midpoints]
      pwd <- pwd[, lapply(.SD, sum), .SDcols = midpoints, by = .(year)]
      annual_pwd <- pwd[, lapply(.SD, mean), .SDcols = midpoints]
      annual_pwd <- melt(annual_pwd, measure.vars = midpoints, variable.name = "bins", value.name = "pop_deg_days", na.rm = TRUE)
      annual_pwd[, bins := ceiling(as.numeric(as.character(bins)))]
      response <- merge(response, annual_pwd, by = "bins", all.x = T)
      
      rm(pwd, annual_pwd)
      gc()
      
      # create side indicator in relation to mmt
      current_narrow_panel[, side := fifelse(temp1 > mmt, "hot", "cold")]
      
      #   12. Prepare aggregated panel  
      current_narrow_panel <- current_narrow_panel[!is.na(side), .(
        temp1 = sum(temp1), 
        temp2 = sum(temp2), 
        temp3 = sum(temp3), 
        temp4 = sum(temp4), 
        ndays = .N
      ), by = c("adm_id", "year", time_agg_var, "side")]
      
      #   12. Merge Population/Death Data  
      current_narrow_panel <- merge(current_narrow_panel, adm_year_pop_death, by = c("adm_id", "year", time_agg_var))
      
      #   13. Death Prediction via Bootstrap  
      temp_cols <- paste0("temp", 1:4)
      pred_mat <- as.matrix(current_narrow_panel[, ..temp_cols])
      pred_mat <- pred_mat - matrix(current_narrow_panel$ndays %o% t(mmt_matrix), ncol = 4)  # Scale by days
      bstrap_deaths <- pred_mat %*% t(bst_coefs_sum)
      
      #   14. Number of estimated deaths   
      if (log_rate == TRUE) {
        # transform back to level
        typical_rate <- mean(current_mort %>% pull(outcome_var), na.rm = TRUE)
        bstrap_deaths <- sweep(expm1(bstrap_deaths), 1, current_narrow_panel$pop, "*") * typical_rate
      } else (
        bstrap_deaths <- sweep(bstrap_deaths, 1, current_narrow_panel$pop, "*")
      )
      
      # aggregate annual deaths
      annual_total_deaths <- adm_year_pop_death[, .(total_deaths = sum(n_deaths, na.rm = TRUE)), by = .(year)]

      # deaths by year and side: group row indices, colSums bootstrap matrix per group
      groups <- current_narrow_panel[, .(row_idx = list(.I)), by = .(year, side)]
      deaths_by_year_side <- rbindlist(lapply(seq_len(nrow(groups)), function(i) {
        boot_sums <- colSums(bstrap_deaths[groups$row_idx[[i]], , drop = FALSE], na.rm = TRUE)
        data.table(
          year        = groups$year[i],
          side        = groups$side[i],
          deaths      = mean(boot_sums, na.rm = TRUE),
          deaths_q025 = quantile(boot_sums, 0.025, na.rm = TRUE, names = FALSE),
          deaths_q975 = quantile(boot_sums, 0.975, na.rm = TRUE, names = FALSE)
        )
      }))
      rm(bstrap_deaths, groups)
      deaths_by_year_side[, `:=`(key = current_quantile, model = model_name)]
      deaths_by_year_side <- merge(deaths_by_year_side, annual_total_deaths, by = "year", all.x = TRUE)
      
      # deaths by side 
      deaths <- deaths_by_year_side[, .(          # Aggregate by year and hot/cold
        total_deaths = mean(total_deaths, na.rm = TRUE),  # Average across years
        deaths_q025 = mean(deaths_q025, na.rm = TRUE),
        deaths = mean(deaths, na.rm = TRUE),
        deaths_q975 = mean(deaths_q975, na.rm = TRUE)
      ), by = .(side)][, key := current_quantile][, model := model_name] 
      
      
      list(response = response, deaths = deaths, deaths_by_year_side = deaths_by_year_side)
        
      })
    
    final <- Filter(Negate(is.null), final)
    out <- list(
      response = rbindlist(lapply(final, `[[`, "response")), 
      deaths = rbindlist(lapply(final, `[[`, "deaths")), 
      deaths_by_year_side = rbindlist(lapply(final, `[[`, "deaths_by_year_side"))
      )
    
    } else if (model_name == "decade") {
      
      ########################################################################
      # DECADE ANALYSIS
      ########################################################################
      
      decades <- as.numeric(as.character(unique(mort$decade)))
      pop_var <- "pop_all"              # Default: All populations
      death_var <- "n_deaths_all"
      outcome_var <- "rate_all"
      
      decade_results <- lapply(decades, function(decade_key) {
        
        # decade_key <- 2000
        print(paste("Processing decade:", decade_key))
      
        current_mort <- mort[decade == decade_key]  
        
        if (country_name == "EU") {
          time_agg_var <- "week"
          adm_year_pop_death <- current_mort[, .(adm_id, year, week, pop = get(pop_var), n_deaths = get(death_var))]
        } else {
          time_agg_var <- "month"
          adm_year_pop_death <- current_mort[, .(adm_id, year, month, pop = get(pop_var), n_deaths = get(death_var))]
        }
        
        current_narrow_panel <- daily_temp[decade == decade_key]
        current_observed_temp <- current_narrow_panel$temp1
        
        decade_cols <- names(coefs)[grepl("temp", names(coefs)) & grepl(as.character(decade_key), names(coefs))]
        
        coefs_before_sum <- t(as.matrix(coefs[decade_cols]))
        # Cumulative effect of lags and keep matrix format
        poly_groups <- c("temp_d1", "temp_d2", "temp_d3", "temp_d4")
        coefs_sum <- sapply(poly_groups, function(poly) {
          poly_names <- colnames(coefs_before_sum)[grepl(poly, colnames(coefs_before_sum)) & !grepl("prec", colnames(coefs_before_sum))]
          rowSums(coefs_before_sum[, poly_names, drop = FALSE])
        })
        
        coefs_sum <- t(as.matrix(coefs_sum))
        pred_vector <- pred_temps_poly %*% t(coefs_sum)
        
        # Crop predictions by 1% and 99% percentile of observed temp distribution
        t1_index <- quantile(current_observed_temp, 0.01, na.rm = TRUE) 
        t99_index <- quantile(current_observed_temp, 0.998, na.rm = TRUE)
        mmt_options <- which(pred_temp_range > t1_index & pred_temp_range < t99_index)
        new_temp_range <- pred_temp_range[mmt_options]
        
        mmt <- new_temp_range[which.min(pred_vector[mmt_options])]
        mmt_matrix <- poly(mmt, degree = 4, raw = TRUE)
        pred_temps_scaled <- scale(pred_temps_poly, center = as.numeric(mmt_matrix), scale = FALSE)
        
        # Bootstrap coefficients
        temp_cols <- names(bstrap)[grepl("temp", names(bstrap)) & grepl(as.character(decade_key), names(bstrap))]
        bst_coefs_before_sum <- bstrap[, ..temp_cols]
        bst_coefs_before_sum <- as.matrix(bst_coefs_before_sum)
        
        # Sum coefficients by poly group (vectorized)
        poly_groups <- c("temp_d1", "temp_d2", "temp_d3", "temp_d4")
        bst_coefs_sum <- sapply(poly_groups, function(poly) {
          poly_names = colnames(bst_coefs_before_sum)[grepl(poly, colnames(bst_coefs_before_sum)) & !grepl("prec", colnames(bst_coefs_before_sum))]
          rowSums(bst_coefs_before_sum[, poly_names, drop = FALSE])
        })
        
        bstrap_result <- pred_temps_scaled %*% t(bst_coefs_sum)
        
        # Calculate response statistics (vectorized)
        response_stats <- cbind(
          matrixStats::rowQuantiles(bstrap_result, probs = c(0.025, 0.5, 0.975), na.rm = TRUE),
          mean = matrixStats::rowMeans2(bstrap_result, na.rm = TRUE)
        )
        
        response <- data.table(
          pred_q025 = response_stats[, 1],
          pred = response_stats[, 4],  # Using mean
          pred_q975 = response_stats[, 3],
          bins = pred_temp_range
        )[, `:=` (
          key = as.character(decade_key),
          model = model_name,
          t1 = t1_index,
          t99 = t99_index
        )]
        
        #   9. Add annual average person days exposure
        midpoints <- (pred_temp_range[-1] + pred_temp_range[-length(pred_temp_range)]) / 2
        midpoints <- as.character(midpoints)
        cols <- c("adm_id", "year", time_agg_var, pop_var)
        pwd <- merge( current_mort[, ..cols], temp_bins, by = c("adm_id", "year", time_agg_var) )
        setnames(pwd, pop_var, "pop")
        pwd <- pwd[!is.na(pop)]
        
        pwd[, (midpoints) := lapply(.SD, function(x) x * pop), .SDcols = midpoints]
        pwd <- pwd[, lapply(.SD, sum), .SDcols = midpoints, by = .(year)]
        annual_pwd <- pwd[, lapply(.SD, mean), .SDcols = midpoints]
        annual_pwd <- melt(annual_pwd, measure.vars = midpoints, variable.name = "bins", value.name = "pop_deg_days", na.rm = TRUE)
        annual_pwd[, bins := ceiling(as.numeric(as.character(bins)))]
        response <- merge(response, annual_pwd, by = "bins", all.x = T)
        
        rm(pwd, annual_pwd)
        gc()
        
        # create side indicator in relation to mmt
        current_narrow_panel[, side := fifelse(temp1 > mmt, "hot", "cold")]
        
        #   12. Prepare aggregated panel  
        current_narrow_panel <- current_narrow_panel[!is.na(side), .(
          temp1 = sum(temp1), 
          temp2 = sum(temp2), 
          temp3 = sum(temp3), 
          temp4 = sum(temp4), 
          ndays = .N
        ), by = c("adm_id", "year", time_agg_var, "side")]
        
        #   12. Merge Population/Death Data  
        current_narrow_panel <- merge(current_narrow_panel, adm_year_pop_death, by = c("adm_id", "year", time_agg_var))
        
        #   13. Death Prediction via Bootstrap  
        temp_cols <- paste0("temp", 1:4)
        pred_mat <- as.matrix(current_narrow_panel[, ..temp_cols])
        pred_mat <- pred_mat - matrix(current_narrow_panel$ndays %o% t(mmt_matrix), ncol = 4)  # Scale by days
        bstrap_deaths <- pred_mat %*% t(bst_coefs_sum)
        
        #   14. Number of estimated deaths   
        if (log_rate == TRUE) {
          # transform back to level
          typical_rate <- mean(current_mort %>% pull(outcome_var), na.rm = TRUE)
          bstrap_deaths <- sweep(expm1(bstrap_deaths), 1, current_narrow_panel$pop, "*") * typical_rate
        } else (
          bstrap_deaths <- sweep(bstrap_deaths, 1, current_narrow_panel$pop, "*")
        )
        
        # aggregate annual deaths
        annual_total_deaths <- adm_year_pop_death[, .(total_deaths = sum(n_deaths, na.rm = TRUE)), by = .(year)]

        # deaths by year and side: group row indices, colSums bootstrap matrix per group
        groups <- current_narrow_panel[, .(row_idx = list(.I)), by = .(year, side)]
        deaths_by_year_side <- rbindlist(lapply(seq_len(nrow(groups)), function(i) {
          boot_sums <- colSums(bstrap_deaths[groups$row_idx[[i]], , drop = FALSE], na.rm = TRUE)
          data.table(
            year        = groups$year[i],
            side        = groups$side[i],
            deaths      = mean(boot_sums, na.rm = TRUE),
            deaths_q025 = quantile(boot_sums, 0.025, na.rm = TRUE, names = FALSE),
            deaths_q975 = quantile(boot_sums, 0.975, na.rm = TRUE, names = FALSE)
          )
        }))
        rm(bstrap_deaths, groups)
        deaths_by_year_side[, `:=`(key = as.character(decade_key), model = model_name)]
        deaths_by_year_side <- merge(deaths_by_year_side, annual_total_deaths, by = "year", all.x = TRUE)
        
        # deaths by side 
        deaths <- deaths_by_year_side[, .(          # Aggregate by year and hot/cold
          total_deaths = mean(total_deaths, na.rm = TRUE),  # Average across years
          deaths_q025 = mean(deaths_q025, na.rm = TRUE),
          deaths = mean(deaths, na.rm = TRUE),
          deaths_q975 = mean(deaths_q975, na.rm = TRUE)
        ), by = .(side)][, key := as.character(decade_key)][, model := model_name] 
        
        
        # Fixed-1970 counterfactual: reuse pred_mat + mmt_matrix from this iteration;
        # only swap bootstrap coefs to 1970. US-only, non-log-rate.
        deaths_fixed1970              <- NULL
        deaths_by_year_side_fixed1970 <- NULL
        if (country_name == "US" & !log_rate) {
          temp_cols_1970 <- names(bstrap)[grepl("temp", names(bstrap)) & grepl("1970", names(bstrap))]
          bst_coefs_before_sum_1970 <- as.matrix(bstrap[, ..temp_cols_1970])
          bst_coefs_sum_1970 <- sapply(poly_groups, function(poly) {
            poly_names <- colnames(bst_coefs_before_sum_1970)[grepl(poly, colnames(bst_coefs_before_sum_1970)) & !grepl("prec", colnames(bst_coefs_before_sum_1970))]
            rowSums(bst_coefs_before_sum_1970[, poly_names, drop = FALSE])
          })

          bstrap_deaths_1970 <- sweep(pred_mat %*% t(bst_coefs_sum_1970), 1, current_narrow_panel$pop, "*")

          groups_1970 <- current_narrow_panel[, .(row_idx = list(.I)), by = .(year, side)]
          deaths_by_year_side_fixed1970 <- rbindlist(lapply(seq_len(nrow(groups_1970)), function(i) {
            boot_sums <- colSums(bstrap_deaths_1970[groups_1970$row_idx[[i]], , drop = FALSE], na.rm = TRUE)
            data.table(
              year        = groups_1970$year[i],
              side        = groups_1970$side[i],
              deaths      = mean(boot_sums, na.rm = TRUE),
              deaths_q025 = quantile(boot_sums, 0.025, na.rm = TRUE, names = FALSE),
              deaths_q975 = quantile(boot_sums, 0.975, na.rm = TRUE, names = FALSE)
            )
          }))
          rm(bstrap_deaths_1970, groups_1970)
          deaths_by_year_side_fixed1970[, `:=`(key = as.character(decade_key), model = model_name)]
          deaths_by_year_side_fixed1970 <- merge(deaths_by_year_side_fixed1970, annual_total_deaths, by = "year", all.x = TRUE)

          deaths_fixed1970 <- deaths_by_year_side_fixed1970[, .(
            total_deaths = mean(total_deaths, na.rm = TRUE),
            deaths_q025  = mean(deaths_q025,  na.rm = TRUE),
            deaths       = mean(deaths,        na.rm = TRUE),
            deaths_q975  = mean(deaths_q975,  na.rm = TRUE)
          ), by = .(side)][, key := as.character(decade_key)][, model := model_name]
        }

        list(response = response, deaths = deaths, deaths_by_year_side = deaths_by_year_side,
             deaths_fixed1970 = deaths_fixed1970, deaths_by_year_side_fixed1970 = deaths_by_year_side_fixed1970)
      })

      # Combine results from all decades
      out <- list(
        response            = rbindlist(lapply(decade_results, `[[`, "response")),
        deaths              = rbindlist(lapply(decade_results, `[[`, "deaths")),
        deaths_by_year_side = rbindlist(lapply(decade_results, `[[`, "deaths_by_year_side"))
      )

      # Save fixed-1970 counterfactual (US only, non-log-rate)
      if (country_name == "US" & !log_rate) {
        decade_fixed1970_out <- list(
          deaths              = rbindlist(Filter(Negate(is.null), lapply(decade_results, `[[`, "deaths_fixed1970"))),
          deaths_by_year_side = rbindlist(Filter(Negate(is.null), lapply(decade_results, `[[`, "deaths_by_year_side_fixed1970")))
        )
        decade_fixed1970_out_file <- file.path(processed_country_dir, paste0(country_name, "_", model_name, "_fixed1970_results.rds"))
        saveRDS(decade_fixed1970_out, decade_fixed1970_out_file)
      }
  }
  
  # Save output for each model
  output_file <- file.path(processed_country_dir, paste0(country_name, "_", model_name, "_results.rds"))
  message("Saving results to: ", output_file)
  saveRDS(out, output_file)
}

################################################################################
# END OF SCRIPT
################################################################################