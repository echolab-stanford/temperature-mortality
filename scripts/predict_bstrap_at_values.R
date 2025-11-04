#' Predict Bootstrap at Values
#'
#' A function to calculate predictions or derivatives using bootstrap coefficients
#' for specific values of a predictor variable (`at_values`) across interaction levels.
#'
#' @param coefs_btstrap A matrix or data frame of bootstrap coefficients.
#' @param at_values A numeric vector of values at which predictions or derivatives are to be calculated.
#' @param add_identity Logical. If TRUE, adds an intercept (identity) term to the predictions. Default is FALSE.
#' @param derivative Logical. If TRUE, calculates derivatives instead of raw predictions. Default is FALSE.
#' @param center Optional. A single numeric value. If provided, centers the predictions at this value.
#' @param interacted_levels A vector of interaction levels to calculate marginal effects for.
#' @param poly_degree Optional. An integer specifying the degree of the polynomial terms for the predictor.
#'
#' @return A matrix containing marginal effects or predictions at specified `at_values` for each interaction level.
#' @examples
#' # Example usage
#' coefs <- matrix(runif(100), nrow = 10) # Example bootstrap coefficients
#' at_vals <- seq(1, 10, by = 1)
#' levels <- c("level1", "level2")
#' result <- predict_bstrap_at_values(coefs, at_vals, FALSE, FALSE, NULL, levels, 3)
#'
#' @export
predict_bstrap_at_values <- function(coefs_btstrap, 
                                     at_values,
                                     add_identity = FALSE, 
                                     derivative = FALSE, 
                                     center = NULL, 
                                     interacted_levels, 
                                     poly_degree = NULL) {
  
  # Step 1: Create the matrix for the predictor variable, either as polynomial terms or as derivatives
  if (!is.null(poly_degree)) {
    if (derivative) {
      # Create a derivative matrix
      xxr <- t(sapply(1:poly_degree, function(deg) deg * at_values^(deg - 1)))
    } else {
      # Create a polynomial matrix
      xxr <- t(sapply(1:poly_degree, function(deg) at_values^deg))
    }
  } else {
    # Treat as a single variable if no polynomial terms are specified
    xxr <- matrix(at_values, nrow = 1)
  }
  
  # Step 2: Calculate marginal effects across specified interaction levels
  marginal_values <- lapply(interacted_levels, function(x) {
    # Prepare coefficients for the current interaction level
    if (add_identity) {
      xxr <- rbind(1, xxr)  # Add intercept term
      coefs_matrix <- as.data.frame(coefs_btstrap) %>%
        dplyr::select(dplyr::contains(as.character(x))) %>%
        as.matrix()
    } else {
      coefs_matrix <- as.data.frame(coefs_btstrap) %>%
        dplyr::select(dplyr::contains(as.character(x)) & !contains("pw_prec")) %>%
        as.matrix()
    }
    
    # Perform matrix multiplication to calculate predictions
    yy <- coefs_matrix %*% xxr
    
    # Center predictions if `center` is provided
    if (!is.null(center)) {
      yy <- t(apply(yy, 1, function(row) row - row[at_values == center]))
    }
    
    # Append interaction level and assign column names
    yy <- cbind(yy, x)
    colnames(yy) <- c(as.character(at_values), "interacted_levels")
    return(yy)
  })
  
  # Step 3: Combine all marginal effect matrices into one final matrix
  all_marginals <- do.call(rbind, marginal_values)
  
  # Step 4: Return the final matrix
  return(all_marginals)
}
