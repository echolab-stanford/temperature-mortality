find_mmt_and_scale_pred_temp <- function(coefs, 
                                         observed_temps, 
                                         degree = 4, 
                                         temp_range, 
                                         mmt_quantiles, 
                                         additional_grep = NULL) {
  
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
  
  # add estimated coefficients to each polynomials
  poly1_coefs <- sum(temp_coefs[grepl("temp_d1", names(temp_coefs))])
  poly2_coefs <- sum(temp_coefs[grepl("temp_d2", names(temp_coefs))])
  poly3_coefs <- sum(temp_coefs[grepl("temp_d3", names(temp_coefs))])
  poly4_coefs <- sum(temp_coefs[grepl("temp_d4", names(temp_coefs))])
  
  # append new coefficients while preserving the original named vector format
  temp_coefs <- c(temp_d1 = poly1_coefs, temp_d2 = poly2_coefs, 
                  temp_d3 = poly3_coefs, temp_d4 = poly4_coefs)
  
  pred_vector <- c(pred_temps_poly %*% temp_coefs)
  
  
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

run_bootstrap_mortality_temp_analysis <- function(boots_results, 
                                                  pred_temps_scaled,
                                                  mmt,
                                                  temp_range,
                                                  pwd, 
                                                  log.rate = FALSE, 
                                                  model_name, 
                                                  model_level,
                                                  country_name,
                                                  processed_country_dir, 
                                                  coef_pattern = "temp", 
                                                  additional_grep = NULL, 
                                                  save_results = T) {
  
  # boots_results = bstrap
  # pred_temps_scaled = mmt_scaled_pred_temp$pred_temps_scaled
  # mmt = mmt_scaled_pred_temp$mmt
  # temp_range = -50:50
  # pwd = pwd
  # model_name = model_name
  # model_level = model_level
  # country_name = country_name
  # processed_country_dir = processed_country_dir
  # coef_pattern = coef_pattern
  # additional_grep = NULL
  # save_results = TRUE
  # log.rate = TRUE
  
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
  # temp_range = -50:50
  # model_name = model_name
  # model_level = level_name
  # country_name = country_name
  # processed_country_dir = processed_country_dir
  # coef_pattern = "pw_tmean"
  # additional_grep = NULL
  # log.rate = TRUE
  
  boots_results <- as.data.frame(boots_results)
  
  # Define polynomial term names
  poly_names <- c("temp_d1", "temp_d2", "temp_d3", "temp_d4")
  
  # Extract coefficients & compute row sums for each polynomial term
  boots_coefs <- sapply(poly_names, function(poly) {
    
    poly_coefs <- boots_results[grepl(poly, names(boots_results))]
    
    if (ncol(poly_coefs) == 1) {
      return(as.vector(poly_coefs[[1]]))
    } else {
      return(rowSums(poly_coefs))
    }
  })
  
  # Ensure boots_coefs is a matrix
  boots_coefs <- if (is.vector(boots_coefs)) matrix(boots_coefs, ncol = length(poly_names)) else boots_coefs
  
  # Compute bootstrapped predictions: each row is an iteration
  boots_pred <- boots_coefs %*% t(pred_temps_scaled)
  
  # Compute mean and confidence intervals for each temperature
  response <- apply(boots_pred, 2, function(coeff) {
    pred_q025 <- quantile(coeff, probs = 0.025, na.rm = TRUE)
    pred <- mean(coeff, na.rm = TRUE)
    pred_q975 <- quantile(coeff, probs = 0.975, na.rm = TRUE)
    return(c(pred_q025 = pred_q025, pred = pred, pred_q975 = pred_q975))
  })

  # Format response function
  response <- as.data.table(t(response))
  names(response) <- c("pred_q025", "pred", "pred_q975")
  response$bins <- temp_range
  response$mmt <- mmt
  response$n_iterations <- nrow(boots_pred)
  
  # clean predicted values at temps 
  boots_pred <- as.data.table(boots_pred)
  colnames(boots_pred) <- as.character(temp_range)
  boots_pred[, iteration_index := 1:.N]
  boots_pred <- melt(boots_pred, id.vars = "iteration_index",
       variable.name = "bins",
       value.name = "pred")
  boots_pred[, mmt := mmt]
  boots_pred[, bins := ceiling(as.numeric(as.character(bins)))] # Ensure bins are numeric

  # we assume predicted rate is same for each year so apply merge like this
  setDT(pwd)
  pwd <- pwd[, bins := ceiling(as.numeric(as.character(bins)))]  # Ensure bins are numeric
  results <- merge(pwd, boots_pred, by = "bins", all.x = TRUE, allow.cartesian = T)
  results[, side := fifelse(bins < mmt, "cold", "hot")]
  
  # indicator whether collapse by adm units
  collapse.adm.unit <- if (model_name == "temp" | model_name == "income" | model_name == "baseline") T else F
  
  if (collapse.adm.unit) {
    results <- results[, .(
      pop_deg_days = sum(pop_deg_days, na.rm = TRUE),
      pw_deg_days = mean(pw_deg_days, na.rm = TRUE), # pw_deg_days is population weighted degrees
      annual_total_pop = sum(annual_total_pop, na.rm = TRUE),
      annual_total_death = sum(annual_total_death, na.rm = TRUE),
      pred = mean(pred),
      mmt = mean(mmt),
      side = first(side),
      n_adm = .N
    ), by = .(bins, year, model_level, iteration_index)]

    # annual_total_pop_death <- pwd[, .(
    #   annual_total_pop = sum(annual_total_pop, na.rm = TRUE),
    #   annual_total_death = sum(annual_total_death, na.rm = TRUE),
    #   n_adm = .N
    # ), by = .(bins, model_level, year)]
    
    annual_total_pop_death <- unique(results[, .(year, annual_total_pop, annual_total_death)])
  } else {
    annual_total_pop_death <- unique(results[, .(year, annual_total_pop, annual_total_death)])
  }
  
  # response by bins and year and pop deg days
  pop_deg_days <- results[, .(pop_deg_days = mean(pop_deg_days, na.rm = TRUE)), by = .(bins)]
  response <- merge(response, pop_deg_days, by = "bins")
  response <- response[pop_deg_days != 0]
  response[, `:=`(model_name = model_name, model_level = model_level)]
  
  # calculate deaths if log_rate is true
  if (log.rate){
    results[, deaths := pw_deg_days * (exp(pred) - 1)] # pw_deg_days is population weighted temperature degree days; sum(pw_deg_days) per year is equal 30 days
  } else {
    results[, deaths := pop_deg_days * pred] # pop_deg_days is person degree days exposure
  }
  
  unique(results$iteration_index)
  
  # deaths per bins (i.e, collapsing by iterations)
  deaths_with_ci <- results[, .(
    deaths_q025 = quantile(deaths, probs = 0.025, na.rm = TRUE),
    deaths      = mean(deaths, na.rm = TRUE),
    deaths_q975 = quantile(deaths, probs = 0.975, na.rm = TRUE),
    mmt         = first(mmt),
    n_iterations = .N
  ), by = .(year, bins)]
  
  # attach model name and level
  deaths_with_ci[, `:=`(model_name = model_name, model_level = model_level)]
  
  # plot(deaths_with_ci[, .(deaths = mean(deaths, na.rm = T), n_years = .N), by = .(bins)]$bins, deaths_with_ci[, .(deaths = mean(deaths, na.rm = T), n_years = .N), by = .(bins)]$deaths, type = "l")
  # plot(deaths_with_ci[deaths_with_ci$year == 2000, ]$deaths, type = "l")
  
  # calculate deaths attributed to cold or hot temperature
  deaths_temp_related_by_bstrap <- results[, .(
    deaths = sum(deaths, na.rm = T),
    mmt = first(mmt),
    n_bins = .N
  ), by = .(year, side, iteration_index)]
  deaths_temp_related_by_bstrap[, `:=`(model_name = model_name, model_level = model_level)]
  
  # calculate deaths attributed to cold or hot with confidence interval
  deaths_temp_related_with_ci <- deaths_temp_related_by_bstrap[, .(
    deaths_q025 = quantile(deaths, probs = 0.025, na.rm = TRUE),
    deaths = mean(deaths, na.rm = TRUE),
    deaths_q975 = quantile(deaths, probs = 0.975, na.rm = TRUE),
    mmt = first(mmt),
    n_iterations = .N
  ), by = .(year, side)]
  
  deaths_temp_related_with_ci <- merge(deaths_temp_related_with_ci, annual_total_pop_death, by = "year", allow.cartesian = TRUE)
  deaths_temp_related_with_ci[, `:=`(model_name = model_name, model_level = model_level)]
  
  if (save_results) {
    # Save all data frames
    print("Saving results ... ")
    
    response_by_bins_dir <- file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_", tolower(model_name), "_response_by_bins.rds"))
    deaths_by_bins_dir <- file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_", tolower(model_name), "_deaths_by_bins.rds"))
    deaths_by_side_year_iteration <- file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_", tolower(model_name), "_deaths_by_side_year_iteration.rds"))
    deaths_by_side_year <- file.path(processed_country_dir, paste0(country_name, "_MortalityTempAnalysis_", tolower(model_name),"_deaths_by_side_year.rds"))
    
    # Save the data frames
    write_rds(response, response_by_bins_dir)
    write_rds(deaths_with_ci, deaths_by_bins_dir)
    write_rds(deaths_temp_related_by_bstrap, deaths_by_side_year_iteration)
    write_rds(deaths_temp_related_with_ci, deaths_by_side_year)
    
    print(paste0("Results are saved successfully in the directory: ", processed_country_dir))
  }
  
  return(list(
    response_df = response, 
    deaths_by_bins = deaths_with_ci,
    deaths_by_side_iteration = deaths_temp_related_by_bstrap,
    deaths_annual_by_side = deaths_temp_related_with_ci
  ))
}

generate_response_plot <- function(df_subset) {
  
  # Ensure color palette matches unique model levels
  unique_levels <- unique(df_subset$model_level)
  color_palette <- rev(met.brewer("Homer1", length(unique_levels)))
  
  response <- ggplot(df_subset, aes(x = bins)) +
    geom_line(aes(y = pred, color = model_level), linewidth = 1) +
    geom_ribbon(aes(ymin = pred_q025, ymax = pred_q975, fill = model_level), alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    
    scale_color_manual(values = setNames(color_palette, unique_levels)) +
    scale_fill_manual(values = setNames(color_palette, unique_levels)) +
    
    guides(color = guide_legend(title = df_subset$level_new[1]), 
           fill = guide_legend(title = df_subset$level_new[1])) + 
    labs(x = "Temperature (°C)", y = "Change in deaths") +
    theme_minimal() +
    theme(
      strip.text   = element_text(size = 10, face = "bold"), 
      legend.position = c(.40, .80),
      legend.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x  = element_blank(),
      axis.title.x = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      plot.margin = margin(0, 0, 0, 0)
    )
  
  exposure <- df_subset %>% 
    ggplot(aes(x = bins, y = pop_deg_days, fill = model_level)) +
    geom_col(alpha = 0.5, position = "identity", width = 0.8) +
    scale_fill_manual(values = setNames(color_palette, unique_levels)) +
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
  
  return(list(response = response, exposure = exposure))
}