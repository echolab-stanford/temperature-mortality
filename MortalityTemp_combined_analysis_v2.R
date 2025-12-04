################################################################################
# ENVIRONMENT SETUP
################################################################################

# Clear environment
rm(list = ls())
gc()

# Load required libraries
pacman::p_boot()
pacman::p_load(fastverse, tidyverse, arrow, fixest, matrixStats, ISOweek)

# Set configuration
setFixest_nthreads(6, save = TRUE)

################################################################################
# DIRECTORY AND PATH SETUP
################################################################################

# setwd("~/path/to/temperature-mortality")
setwd("~/BurkeLab Dropbox/projects/temperature-mortality")

# Define main parameters
# To replicate analysis, you can set the following parameters
# outputs from this script will be saved in "processed" folder for respective countries
# for combined figures, we will run different script
country_name <- "US" # Other options: EU, MEX
# if one wants to replicate supplmentary figure inputs
# options include: log_rate, same_year, age_standardized, data_1990_2023 (only when country_name == MEX) or winsor 
si_folder <- "winsor" 

# Set up directory structure
data_dir <- 'data'
cleaned_dir <- paste(data_dir, country_name, "cleaned", sep = "/")
processed_dir <- 'processed'
processed_country_dir <- paste(processed_dir, country_name, sep = "/")

# Load helper functions and data
# source(paste(script_dir, "helper_functions.R", sep = "/"))
mort <- read_parquet( paste(cleaned_dir, paste0(country_name, "FullPanel_complete.pq"), sep = "/") )
temp_bins <- read_parquet( paste(cleaned_dir, paste0(country_name, "TempBins_monthly_complete.pq"), sep = "/") )
daily_temp <- read_parquet( paste(cleaned_dir, paste0(country_name, "Daily_temp_covariates.pq"), sep = "/") )
setDT(mort)
setDT(temp_bins)
setDT(daily_temp)

################################################################################
# COUNTRY-SPECIFIC DATA PREPARATION
################################################################################

if (country_name == "EU") {
  
  setkeyv(mort, c("adm_id", "year", "week"))
  setkeyv(temp_bins, c("adm_id", "year", "week"))
  setkeyv(daily_temp, c("adm_id", "year", "month", "day"))
  
  # create month variable
  mort[, date := ISOweek2date(sprintf("%d-W%02d-1", year, week))]
  mort[, month := month(date)]
  
  # Create lag variables
  new_lag_names <- c('temp_d1', 'temp_d2', 'temp_d3', 'temp_d4', 'precip_d1', 'precip_d2')
  weather_polynomial_names <- c('temp_d1_l0','temp_d2_l0','temp_d3_l0','temp_d4_l0', 'precip_d1_l0', 'precip_d2_l0')
  lags <- 1:4
  
  # Generate lag variables using data.table
  for (lag_value in lags){
    lag_cols <- paste0(new_lag_names, '_t', lag_value)
    mort[, (lag_cols) := lapply(.SD, shift, n = lag_value), by = adm_id, .SDcols = weather_polynomial_names]
  }
  
  # Define model variables for EU
  right_hand_var_list <- c("temp_d1_l0", "temp_d1_t1", "temp_d1_t2", "temp_d1_t3",
                           "temp_d2_l0", "temp_d2_t1", "temp_d2_t2", "temp_d2_t3",  
                           "temp_d3_l0", "temp_d3_t1", "temp_d3_t2", "temp_d3_t3",  
                           "temp_d4_l0", "temp_d4_t1", "temp_d4_t2", "temp_d4_t3",  
                           "precip_d1_l0", "precip_d1_t1", "precip_d1_t2", "precip_d1_t3",  
                           "precip_d2_l0", "precip_d2_t1", "precip_d2_t2", "precip_d2_t3")
  right_hand_vars <- paste(right_hand_var_list, collapse = " + ")
  fixed_effects <- "adm_id^year + adm_id^week"
  
} else if (country_name %in% c("US", "MEX")) {
  
  setkeyv(mort, c("adm_id", "year", "month"))
  setkeyv(temp_bins, c("adm_id", "year", "month"))
  setkeyv(daily_temp, c("adm_id", "year", "month", "day"))
  
  # # Precompute days_in_month outside the loop
  # mort[, ndays := days_in_month(paste0(year, '-', month, '-01'))]
  
  # Define model variables for US
  right_hand_var_list <- c("temp_d1_l0", 
                           "temp_d2_l0", 
                           "temp_d3_l0", 
                           "temp_d4_l0", 
                           "precip_d1_l0",  
                           "precip_d2_l0")
  right_hand_vars <- paste(right_hand_var_list, collapse = " + ")
  fixed_effects <- "adm_id^year + adm_id^month"
  
} else {
  message("Country must be: US, EU, or MEX")
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
  
  # Check if country_name equals "US" or "MEX"
  if (country_name %in% c("US", "MEX")) {
    # Set up age-standardized analysis
    processed_country_dir <- file.path(processed_country_dir, si_folder)
    
    if (!dir.exists(processed_country_dir)) {
      dir.create(processed_country_dir, recursive = FALSE)
    }
    
    mort <- mort %>%
      mutate(rate_all = age_std_rate) %>%
      filter(!is.na(age_std_rate)) %>%
      filter(is.finite(age_std_rate)) # dropping NA and infinity
    log_rate <- FALSE
  } else {
    message("Age-standardized analysis can only be run for US or MEX.")
  }

  
}  else if (si_folder == "newdata_agerate_fixed" | si_folder == "newdata_original" | si_folder == "newdata_1990_2023") {
  
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

# Filter models for same_year analysis
if (si_folder == "same_year") {
  model_formula_list <- model_formula_list[1]
}

# Main loop through models
for (model_name in names(model_formula_list)[8:8]) { 

  # model_name <- "baseline"
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
    
    # --- 0. Population-Specific Variables ---
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
    
    # --- 1. Model Coefficient Processing ---
    # Transpose coefficient matrix for downstream operations
    coefs_before_sum <- t(as.matrix(coefs))  # coefs: Original regression coefficients
    
    # --- 2. Cumulative Lag Effects Calculation ---
    # Sum coefficients for each polynomial term (temp_d1 to temp_d4), excluding precipitation terms
    poly_groups <- c("temp_d1", "temp_d2", "temp_d3", "temp_d4")
    coefs_sum <- sapply(poly_groups, function(poly) {
      # Identify columns matching the polynomial group (e.g., temp_d1_lag1, temp_d1_lag2)
      poly_names = colnames(coefs_before_sum)[grepl(poly, colnames(coefs_before_sum)) & !grepl("prec", colnames(coefs_before_sum))]
      rowSums(coefs_before_sum[, poly_names, drop = FALSE])  # Sum across lags
    })
    
    coefs_sum <- t(as.matrix(coefs_sum))  # Transpose back to original orientation
    
    # --- 3. Temperature Response Curve ---
    # Predict mortality response using polynomial-transformed temperatures
    pred_vector <- pred_temps_poly %*% t(coefs_sum)  # pred_temps_poly: Temp values transformed via polynomial basis
    
    # --- 4. Temperature Range Trimming ---
    # Restrict predictions to 1%-99% percentile of observed temperatures to exclude extremes
    observed_temps <- daily_temp$temp1  # Observed temperature data
    t1_index <- collapse::fquantile(observed_temps, 0.01, na.rm = TRUE)  # 1st percentile
    t99_index <- collapse::fquantile(observed_temps, 0.998, na.rm = TRUE)  # 99th percentile
    mmt_options <- which(pred_temp_range > t1_index & pred_temp_range < t99_index)  # Indices within range
    new_temp_range <- pred_temp_range[mmt_options]  # Subset temperature range
    
    # --- 5. Minimum Mortality Temperature (MMT) ---
    # Find temperature with lowest predicted mortality (MMT)
    mmt <- new_temp_range[which.min(pred_vector[mmt_options])]
    mmt_matrix <- poly(mmt, degree = 4, raw = T)  # Create polynomial terms for MMT
    pred_temps_scaled <- scale(pred_temps_poly, center = mmt_matrix, scale = FALSE)  # Center temps at MMT
    
    pred_vector <- pred_temps_scaled %*% t(coefs_sum)
    
    # --- 6. Bootstrap Coefficient Processing ---
    # Extract temperature-related coefficients from bootstrap results
    temp_cols <- names(bstrap)[grep("temp", names(bstrap))]  # Columns with "temp" in name
    bst_coefs_before_sum <- bstrap[, ..temp_cols]  # Subset bootstrap coefficients
    bst_coefs_before_sum <- as.matrix(bst_coefs_before_sum)  # Convert to matrix
    
    # --- 7. Bootstrap Response Calculation ---
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
    
    # --- 8. Response Result ---
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
    
    # --- 9. Add annual average person days exposure
    midpoints <- (pred_temp_range[-1] + pred_temp_range[-length(pred_temp_range)]) / 2
    midpoints <- as.character(midpoints)
    
    cols <- c("adm_id", "year", time_agg_var, pop_var)
    pwd <- merge( mort[, ..cols], temp_bins, by = c("adm_id", "year", time_agg_var) )
    setnames(pwd, pop_var, "pop")
    annual_pwd <- pwd[!is.na(pop)][, 
                                   (midpoints) := lapply(.SD, function(x) x * pop), .SDcols = midpoints
                                   ][, 
                                     lapply(.SD, sum), .SDcols = midpoints, by = .(year)
                                     ][, 
                                       lapply(.SD, mean), .SDcols = midpoints
                                       ]
    annual_pwd <- melt(annual_pwd, measure.vars = midpoints, variable.name = "bins", value.name = "pop_deg_days", na.rm = TRUE)
    annual_pwd[, bins := ceiling(as.numeric(as.character(bins)))]
    response <- merge(response, annual_pwd, by = "bins", all.x = T)
    
    rm(pwd, annual_pwd)
    gc()

    
    # --- 11. Hot/Cold Day Classification ---
    # Classify days as hot (>MMT) or cold (<=MMT)
    narrow_panel <- daily_temp[, side := fifelse(temp1 > mmt, "hot", "cold")]
    
    # --- 12. Prepare aggregated panel ---
    narrow_panel <- narrow_panel[!is.na(side), .(
      temp1 = sum(temp1), 
      temp2 = sum(temp2), 
      temp3 = sum(temp3), 
      temp4 = sum(temp4), 
      ndays = .N
    ), by = c("adm_id", "year", time_agg_var, "side")]
    
    # --- 12. Merge Population/Death Data ---
    narrow_panel <- merge(narrow_panel, mort %>% select(adm_id, year, all_of(time_agg_var), pop = all_of(pop_var), n_deaths = all_of(death_var)), by = c("adm_id", "year", time_agg_var))
    
    # --- 13. Death Prediction via Bootstrap ---
    temp_cols <- paste0("temp", 1:4)
    pred_mat <- as.matrix(narrow_panel[, ..temp_cols])
    pred_mat <- pred_mat - matrix(narrow_panel$ndays %o% t(mmt_matrix), ncol = 4)  # Scale by days
    bstrap_deaths <- pred_mat %*% t(bst_coefs_sum)
    
    # --- 14. Number of estimated deaths --- 
    if (log_rate == TRUE) {
      # transform back to level
      typical_rate <- mean(mort %>% group_by(get(time_agg_var)) %>% pull(outcome_var), na.rm = TRUE)
      # bstrap_deaths <- sweep(expm1(bstrap_deaths), 1, narrow_panel$pop, "*") * typical_rate
      bstrap_deaths <- expm1(bstrap_deaths)*narrow_panel$pop*typical_rate
    } else (
      bstrap_deaths <- sweep(bstrap_deaths, 1, narrow_panel$pop, "*")
    )
    
    # create narrow panel and bootstrap panel
    bstrap_deaths <- as.data.table(bstrap_deaths)
    colnames(bstrap_deaths) <- paste0("b", 1:100)
    narrow_panel <- cbind(narrow_panel, bstrap_deaths)
    
    # aggregate annual deaths 
    annual_total_deaths <- mort[, .(total_deaths = sum(get(death_var), na.rm = TRUE)), by = .(year)]
    
    # deaths by year and side
    deaths_by_year_side <- narrow_panel[, lapply(.SD, sum, na.rm=TRUE), by=.(year, side), .SDcols=colnames(bstrap_deaths)]

    deaths_by_year_side[, deaths := rowMeans(.SD, na.rm=TRUE), .SDcols = colnames(bstrap_deaths)]
    deaths_by_year_side[, deaths_q025 := apply(.SD, 1, quantile, probs=0.025, na.rm=TRUE), .SDcols = colnames(bstrap_deaths)]
    deaths_by_year_side[, deaths_q975 := apply(.SD, 1, quantile, probs=0.975, na.rm=TRUE), .SDcols = colnames(bstrap_deaths)]

    deaths_by_year_side[, key := "pooled"]
    deaths_by_year_side[, model := model_name]
    deaths_by_year_side <- merge(deaths_by_year_side, annual_total_deaths, by=c("year"), all.x=TRUE)
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
    
    # --- 15. Final Output ---
    out <- list(
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
      mmt_matrix <- poly(mmt, degree = 4, raw = T)
      bst_pred_temps_scaled <- scale(pred_temps_poly, center = mmt_matrix, scale = FALSE)
      
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
      
      # --- 9. Add annual average person days exposure
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
      
      # --- 12. Prepare aggregated panel ---
      current_narrow_panel <- current_narrow_panel[!is.na(side), .(
        temp1 = sum(temp1), 
        temp2 = sum(temp2), 
        temp3 = sum(temp3), 
        temp4 = sum(temp4), 
        ndays = .N
      ), by = c("adm_id", "year", time_agg_var, "side")]
      
      # --- 12. Merge Population/Death Data ---
      current_narrow_panel <- merge(current_narrow_panel, adm_year_pop_death, by = c("adm_id", "year", time_agg_var))
      
      # --- 13. Death Prediction via Bootstrap ---
      temp_cols <- paste0("temp", 1:4)
      pred_mat <- as.matrix(current_narrow_panel[, ..temp_cols])
      pred_mat <- pred_mat - matrix(current_narrow_panel$ndays %o% t(mmt_matrix), ncol = 4)  # Scale by days
      bstrap_deaths <- pred_mat %*% t(bst_coefs_sum)
      
      # --- 14. Number of estimated deaths --- 
      if (log_rate == TRUE) {
        # transform back to level
        typical_rate <- mean(current_mort %>% pull(outcome_var), na.rm = TRUE)
        bstrap_deaths <- sweep(expm1(bstrap_deaths), 1, current_narrow_panel$pop, "*") * typical_rate
      } else (
        bstrap_deaths <- sweep(bstrap_deaths, 1, current_narrow_panel$pop, "*")
      )
      
      bstrap_deaths <- as.data.table(bstrap_deaths)
      colnames(bstrap_deaths) <- paste0("b", 1:100)
      current_narrow_panel <- cbind(current_narrow_panel, bstrap_deaths)
      
      # aggregate annual deaths 
      annual_total_deaths <- adm_year_pop_death[, .(total_deaths = sum(n_deaths, na.rm = TRUE)), by = .(year)]
      
      # deaths by year and side
      deaths_by_year_side <- current_narrow_panel %>% 
        group_by(year, side) %>% 
        summarise(across(starts_with("b"), ~sum(.x, na.rm = T)), .groups = "drop") %>% 
        ungroup() %>% group_by(year, side) %>%  
        summarise(
          deaths = mean(c_across(starts_with("b")), na.rm = TRUE),
          deaths_q025 = quantile(c_across(starts_with("b")), 0.025, na.rm = TRUE),
          deaths_q975 = quantile(c_across(starts_with("b")), 0.975, na.rm = TRUE), 
          .groups = "drop"
        ) %>% 
        mutate(key = current_quantile, model = model_name) %>% left_join(annual_total_deaths) %>% as.data.table()
      
      # deaths by side 
      deaths <- deaths_by_year_side[, .(          # Aggregate by year and hot/cold
        total_deaths = mean(total_deaths, na.rm = TRUE),  # Average across years
        deaths_q025 = mean(deaths_q025, na.rm = TRUE),
        deaths = mean(deaths, na.rm = TRUE),
        deaths_q975 = mean(deaths_q975, na.rm = TRUE)
      ), by = .(side)][, key := current_quantile][, model := model_name] 
      
      
      list(response = response, deaths = deaths, deaths_by_year_side = deaths_by_year_side)
        
      })
    
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
        pred_temps_scaled <- scale(pred_temps_poly, center = mmt_matrix, scale = FALSE)
        
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
        
        # --- 9. Add annual average person days exposure
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
        
        # --- 12. Prepare aggregated panel ---
        current_narrow_panel <- current_narrow_panel[!is.na(side), .(
          temp1 = sum(temp1), 
          temp2 = sum(temp2), 
          temp3 = sum(temp3), 
          temp4 = sum(temp4), 
          ndays = .N
        ), by = c("adm_id", "year", time_agg_var, "side")]
        
        # --- 12. Merge Population/Death Data ---
        current_narrow_panel <- merge(current_narrow_panel, adm_year_pop_death, by = c("adm_id", "year", time_agg_var))
        
        # --- 13. Death Prediction via Bootstrap ---
        temp_cols <- paste0("temp", 1:4)
        pred_mat <- as.matrix(current_narrow_panel[, ..temp_cols])
        pred_mat <- pred_mat - matrix(current_narrow_panel$ndays %o% t(mmt_matrix), ncol = 4)  # Scale by days
        bstrap_deaths <- pred_mat %*% t(bst_coefs_sum)
        
        # --- 14. Number of estimated deaths --- 
        if (log_rate == TRUE) {
          # transform back to level
          typical_rate <- mean(current_mort %>% pull(outcome_var), na.rm = TRUE)
          bstrap_deaths <- sweep(expm1(bstrap_deaths), 1, current_narrow_panel$pop, "*") * typical_rate
        } else (
          bstrap_deaths <- sweep(bstrap_deaths, 1, current_narrow_panel$pop, "*")
        )
        
        bstrap_deaths <- as.data.table(bstrap_deaths)
        colnames(bstrap_deaths) <- paste0("b", 1:100)
        current_narrow_panel <- cbind(current_narrow_panel, bstrap_deaths)
        
        # aggregate annual deaths 
        annual_total_deaths <- adm_year_pop_death[, .(total_deaths = sum(n_deaths, na.rm = TRUE)), by = .(year)]
        
        # deaths by year and side
        deaths_by_year_side <- current_narrow_panel %>% 
          group_by(year, side) %>% 
          summarise(across(starts_with("b"), ~sum(.x, na.rm = T))) %>% 
          ungroup() %>% group_by(year, side) %>%  
          summarise(
            deaths = mean(c_across(starts_with("b")), na.rm = TRUE),
            deaths_q025 = quantile(c_across(starts_with("b")), 0.025, na.rm = TRUE),
            deaths_q975 = quantile(c_across(starts_with("b")), 0.975, na.rm = TRUE), 
            .groups = "drop"
          ) %>% 
          mutate(key = as.character(decade_key), model = model_name) %>% left_join(annual_total_deaths) %>% as.data.table()
        
        # deaths by side 
        deaths <- deaths_by_year_side[, .(          # Aggregate by year and hot/cold
          total_deaths = mean(total_deaths, na.rm = TRUE),  # Average across years
          deaths_q025 = mean(deaths_q025, na.rm = TRUE),
          deaths = mean(deaths, na.rm = TRUE),
          deaths_q975 = mean(deaths_q975, na.rm = TRUE)
        ), by = .(side)][, key := as.character(decade_key)][, model := model_name] 
        
        
        list(response = response, deaths = deaths, deaths_by_year_side = deaths_by_year_side)
      })
      
      # Combine results from all decades
      out <- list(
        response = rbindlist(lapply(decade_results, `[[`, "response")), 
        deaths = rbindlist(lapply(decade_results, `[[`, "deaths")),
        deaths_by_year_side = rbindlist(lapply(decade_results, `[[`, "deaths_by_year_side"))
      )
      
      # Decade attributable deaths if response function was fixed at 1970
      # with decade specific mmt
      if (country_name == "US" & !log_rate) {
        pop_var <- "pop_all"              # Default: All populations
        death_var <- "n_deaths_all"
        
        decade_fixed1970 <- lapply(decades, function(decade_key) {
          
          fixed_decade <- 1970
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
          pred_temps_scaled <- scale(pred_temps_poly, center = mmt_matrix, scale = FALSE)
          
          # Bootstrap coefficients
          temp_cols <- names(bstrap)[grepl("temp", names(bstrap)) & grepl(as.character(fixed_decade), names(bstrap))]
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
          
          # --- 9. Add annual average person days exposure
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
          
          # --- 12. Prepare aggregated panel ---
          current_narrow_panel <- current_narrow_panel[!is.na(side), .(
            temp1 = sum(temp1), 
            temp2 = sum(temp2), 
            temp3 = sum(temp3), 
            temp4 = sum(temp4), 
            ndays = .N
          ), by = c("adm_id", "year", time_agg_var, "side")]
          
          # --- 12. Merge Population/Death Data ---
          current_narrow_panel <- merge(current_narrow_panel, adm_year_pop_death, by = c("adm_id", "year", time_agg_var))
          
          # --- 13. Death Prediction via Bootstrap ---
          temp_cols <- paste0("temp", 1:4)
          pred_mat <- as.matrix(current_narrow_panel[, ..temp_cols])
          pred_mat <- pred_mat - matrix(current_narrow_panel$ndays %o% t(mmt_matrix), ncol = 4)  # Scale by days
          bstrap_deaths <- pred_mat %*% t(bst_coefs_sum)
          
          # --- 14. Number of estimated deaths --- 
          if (log_rate == TRUE) {
            # transform back to level
            typical_rate <- mean(current_mort %>% pull(outcome_var), na.rm = TRUE)
            bstrap_deaths <- sweep(expm1(bstrap_deaths), 1, current_narrow_panel$pop, "*") * typical_rate
          } else (
            bstrap_deaths <- sweep(bstrap_deaths, 1, current_narrow_panel$pop, "*")
          )
          
          bstrap_deaths <- as.data.table(bstrap_deaths)
          colnames(bstrap_deaths) <- paste0("b", 1:100)
          current_narrow_panel <- cbind(current_narrow_panel, bstrap_deaths)
          
          # aggregate annual deaths 
          annual_total_deaths <- adm_year_pop_death[, .(total_deaths = sum(n_deaths, na.rm = TRUE)), by = .(year)]
          
          # deaths by year and side
          deaths_by_year_side <- current_narrow_panel %>% 
            group_by(year, side) %>% 
            summarise(across(starts_with("b"), ~sum(.x, na.rm = T))) %>% 
            ungroup() %>% group_by(year, side) %>%  
            summarise(
              deaths = mean(c_across(starts_with("b")), na.rm = TRUE),
              deaths_q025 = quantile(c_across(starts_with("b")), 0.025, na.rm = TRUE),
              deaths_q975 = quantile(c_across(starts_with("b")), 0.975, na.rm = TRUE), 
              .groups = "drop"
            ) %>% 
            mutate(key = as.character(decade_key), model = model_name) %>% left_join(annual_total_deaths) %>% as.data.table()
          
          # deaths by side 
          deaths <- deaths_by_year_side[, .(          # Aggregate by year and hot/cold
            total_deaths = mean(total_deaths, na.rm = TRUE),  # Average across years
            deaths_q025 = mean(deaths_q025, na.rm = TRUE),
            deaths = mean(deaths, na.rm = TRUE),
            deaths_q975 = mean(deaths_q975, na.rm = TRUE)
          ), by = .(side)][, key := as.character(decade_key)][, model := model_name] 
          
          
          list(response = response, deaths = deaths, deaths_by_year_side = deaths_by_year_side)
        })
        
        decade_fixed1970_out <- list(
          response = rbindlist(lapply(decade_fixed1970, `[[`, "response")), 
          deaths = rbindlist(lapply(decade_fixed1970, `[[`, "deaths")),
          deaths_by_year_side = rbindlist(lapply(decade_fixed1970, `[[`, "deaths_by_year_side"))
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