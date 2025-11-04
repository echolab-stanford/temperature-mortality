#' Run Bootstrap for Multiple FEOLS Models
#'
#' This function performs a bootstrap procedure to estimate fixed effects models (`feols`) using resampled data. 
#' The bootstrap process involves resampling spatial units with replacement, estimating coefficients for the provided model formulas, 
#' and storing the results in a tidy format. If a saved bootstrap file exists, it can either load the file or override it based on user input.
#'
#' @param sample_from A data frame or tibble from which the resampling is performed.
#' @param unit_column A character string indicating the column name of the spatial units to resample by. 
#' These units will be resampled with replacement for each bootstrap iteration.
#' @param Nboot Integer. The number of bootstrap iterations to perform.
#' @param model_list A list of model formulas to be passed to `fixest::feols()`. Each formula specifies a model to be estimated in each iteration.
#' @param weights Optional. A character string specifying the column name of the weights variable in `sample_from`. If `NULL`, no weights are applied.
#' @param model_name A character string specifying the name of the model (used for naming the output file).
#' @param country_name A character string specifying the country name (used for naming the output file).
#' @param processed_country_dir A character string specifying the directory where the bootstrap results will be saved.
#' @param override Logical. If `TRUE`, the function will always run the bootstrap and override any existing file. If `FALSE`, it will load the existing file if it exists.
#'
#' @return A data frame where each row corresponds to one bootstrap iteration, and the columns represent the estimated coefficients for all models across the iterations. 
#' @examples
#' # Define some models
#' model_list <- list(
#'   y ~ x1 + x2,
#'   y ~ x1 + x2 + x3
#' )
#' 
#' # Run bootstrap with 100 iterations
#' results <- run_bootstrap_feols(
#'   sample_from = data, 
#'   unit_column = "region", 
#'   Nboot = 100, 
#'   model_list = model_list, 
#'   weights = "pop_weights",
#'   model_name = "Young",
#'   country_name = "USA",
#'   processed_country_dir = "path/to/dir",
#'   override = FALSE
#' )
#'
#' @import dplyr
#' @import fixest
#' @import foreach
#' @import parallel
#' @import doParallel
#' @import tidyr
#' @import data.table
#' @import arrow
#' @export

run_bootstrap_feols <- function(sample_from, Nboot, formula, weights = NULL, model_name, country_name, processed_country_dir, override = FALSE) {
  
  # Define the output file path
  out_dir <- file.path(processed_country_dir, paste0(country_name, "_BootstrapResult_", tolower(model_name), "_n", Nboot, ".pq"))
  
  # Check if the file exists and override is FALSE
  if (!override && file.exists(out_dir)) {
    message("Loading existing bootstrap results from: ", out_dir)
    return(read_parquet(out_dir))
  }
  
  # If file does not exist or override is TRUE, run the bootstrap
  message("Running bootstrap for model: ", model_name)
  
  # Parallel backend for foreach
  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  doParallel::registerDoParallel(cl)
  
  # Nboot <- 2
  # formula <- pooled_model_formula
  # weights <- NULL
  
  # Bootstrap loop
  results <- foreach(i = 1:Nboot, .combine = 'rbind', .packages = c("fixest")) %dopar% {
    
    # Resample spatial units with replacement
    sampled_units <- sample(unique(mort[['adm_id']]), size = length(unique(mort[['adm_id']])), replace = TRUE)
    dtboot <- mort[adm_id %in% sampled_units]
    # dtboot <- sample_from %>% filter(!!sym(unit_column) %in% sampled_units)
    
    if (is.null(weights)) {
      fit <- fixest::feols(formula, data = dtboot)
    } else {
      fit <- fixest::feols(formula, data = dtboot, weights = dtboot[[weights]])
    }
    model_results <- coef(fit)
    
    # Combine results into a single row
    coefs_vector <- unlist(model_results, recursive = TRUE, use.names = TRUE)
    coefs_df <- as.data.table(t(coefs_vector), stringsAsFactors = FALSE)
    coefs_df$iteration <- i
    coefs_df
  }
  
  # Stop parallel backend
  parallel::stopCluster(cl)
  
  # results <- as.data.table(results)
  # Save the results
  message("Saving bootstrap results to: ", out_dir)
  write_parquet(results, out_dir)
  
  return(results)
}