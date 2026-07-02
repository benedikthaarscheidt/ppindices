#' Coefficient of variation of trait values (CV_t)
#'
#' Computes the coefficient of variation across a genotype's trait values measured
#' over multiple environments.
#'
#' @param trait_values Numeric vector of trait measurements across environments,
#'   with no missing values.
#' @return A single numeric value: the coefficient of variation
#'   (`sd(trait_values) / mean(trait_values)`).
#' @examples
#' calculate_CVt(c(2, 4, 6, 8))
#' @export
calculate_CVt = function(trait_values) {
  if (length(trait_values) < 2) {
    return(NA)  # Avoid division by zero for single values
  }

  if (any(is.na(trait_values))) {
    stop("The input vector contains missing values. Please handle missing data before calling this function.")
  }

  # Calculate and return the CVt
  return(sd(trait_values) / mean(trait_values))
}

################################

#' Coefficient of variation of means (CVm)
#'
#' Computes the coefficient of variation across group (e.g. environment) means of
#' a trait, where groups are supplied either as a list of vectors or via
#' `group_labels`.
#'
#' @param trait_values Either a list of numeric vectors (one per group) or a
#'   numeric vector of trait measurements to be split by `env_values`.
#' @param env_values Optional grouping vector used to split `trait_values` when
#'   `trait_values` is not already a list.
#' @return A single numeric value: the coefficient of variation of the group means.
#' @note The example below is illustrative only and is wrapped in `\dontrun{}`
#'   because the non-list branch of this function references an undefined
#'   `group_labels` object rather than its `env_values` argument.
#' @examples
#' \dontrun{
#' calculate_CVm(list(c(2, 3, 4), c(5, 6, 7), c(1, 2, 3)))
#' }
#' @export
calculate_CVm = function(trait_values, env_values = NULL) {
  # Handle case where trait_values is a list of vectors
  if (is.list(trait_values)) {
    # Calculate group means directly
    group_means = sapply(trait_values, mean, na.rm = TRUE)
  } else if (!is.null(group_labels)) {
    # Calculate group means based on group labels
    group_means = tapply(trait_values, group_labels, mean, na.rm = TRUE)
  } else {
    stop("If trait_values is not a list, group_labels must be provided.")
  }

  # Calculate standard deviation and mean of group means
  sd_of_means = sd(group_means, na.rm = TRUE)
  mean_of_means = mean(group_means, na.rm = TRUE)

  # Calculate and return CVm
  CVm = sd_of_means / mean_of_means
  return(CVm)
}

################################

#' Coefficient of variation of medians (CVmd)
#'
#' Computes the coefficient of variation across group (e.g. environment) medians
#' of a trait, where groups are supplied either as a list of vectors or via
#' `group_labels`.
#'
#' @param trait_values Either a list of numeric vectors (one per group) or a
#'   numeric vector of trait measurements to be split by `group_labels`.
#' @param group_labels Optional grouping vector used to split `trait_values` when
#'   `trait_values` is not already a list.
#' @return A single numeric value: the coefficient of variation of the group medians.
#' @examples
#' calculate_CVmd(c(2, 4, 6, 8, 10, 12), group_labels = c(1, 1, 2, 2, 3, 3))
#' @export
calculate_CVmd = function(trait_values, group_labels = NULL) {
  # Handle case where trait_values is a list of vectors
  if (is.list(trait_values)) {
    # Calculate group medians directly
    group_medians = sapply(trait_values, median, na.rm = TRUE)
  } else if (!is.null(group_labels)) {
    # Calculate group medians based on group labels
    group_medians = tapply(trait_values, group_labels, median, na.rm = TRUE)
  } else {
    stop("If trait_values is not a list, group_labels must be provided.")
  }

  # Calculate standard deviation and mean of group medians
  sd_of_medians = sd(group_medians, na.rm = TRUE)
  mean_of_medians = mean(group_medians, na.rm = TRUE)

  # Calculate and return CVmd
  CVmd = sd_of_medians / mean_of_medians
  return(CVmd)
}

################################

#' Coefficient of environmental variation (CEV)
#'
#' Computes the coefficient of variation of a trait across environments, expressed
#' as a percentage.
#'
#' @param trait_values Numeric vector of trait measurements across environments.
#' @return A single numeric value: `100 * sd(trait_values) / mean(trait_values)`.
#'   Returns a list with `CEV = NA` (plus `Mean`, `SD`, `Valid`) when fewer than
#'   two values are supplied or the mean is zero.
#' @examples
#' calculate_CEV(c(2, 4, 6, 8))
#' @export
calculate_CEV = function(trait_values) {
  # Ensure input is a numeric vector
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }

  # Check that there are enough data points to compute a standard deviation
  if (length(trait_values) < 2) {
    return(list(CEV = NA, Mean = NA, SD = NA, Valid = FALSE))
  }

  # Calculate mean and standard deviation (ignoring NA values)
  mean_val = mean(trait_values, na.rm = TRUE)
  sd_val = sd(trait_values, na.rm = TRUE)

  # Avoid division by zero if the mean is zero
  if (mean_val == 0) {
    return(list(CEV = NA, Mean = mean_val, SD = sd_val, Valid = FALSE))
  }

  # Calculate the Coefficient of Environmental Variation (CEV)
  cev = (sd_val / mean_val) * 100

  return(cev)
}

################################

#' Environmental Variance Sensitivity (EVS)
#'
#' Computes the ratio of trait variance to environmental variance, a measure of
#' how much trait variability is expressed relative to the spread of the
#' environmental gradient.
#'
#' @param trait_values Numeric vector of trait measurements across environments.
#' @param env Optional numeric or factor vector of environment values, the same
#'   length as `trait_values`. Defaults to equidistant indices `1, 2, ..., n`.
#' @return A single numeric value: `var(trait_values) / var(env)`. Returns `NA`
#'   with a warning when the environmental variance is zero.
#' @examples
#' calculate_EVS(c(2, 4, 6, 8), env = c(1, 2, 3, 4))
#' @export
calculate_EVS <- function(trait_values, env = NULL) {
  # Validate trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }

  # If no environment vector is provided, assume equidistant environments.
  if (is.null(env)) {
    env <- seq_along(trait_values)
  } else {
    if (length(trait_values) != length(env)) {
      stop("trait_values and env must have the same length.")
    }
  }

  # Compute the variance in trait values.
  trait_variance <- var(trait_values, na.rm = TRUE)

  # If env is not numeric, convert it to numeric via as.factor.
  if (!is.numeric(env)) {
    env_numeric <- as.numeric(as.factor(env))
  } else {
    env_numeric <- env
  }

  # Compute the variance in the environment.
  env_variance <- var(env_numeric, na.rm = TRUE)

  if (env_variance == 0) {
    warning("Environmental variance is zero; EVS is undefined. Returning NA.")
    return(NA)
  }

  # Calculate EVS as the ratio of the two variances.
  EVS_value <- trait_variance / env_variance

  return(EVS_value)
}

################################

#' Plasticity stability Index (PSI)
#'
#' Fits a linear regression of trait values on environmental values and converts
#' the slope into a bounded stability score.
#'
#' @param trait_values Numeric vector of trait measurements across environments.
#' @param env_values Optional numeric vector of environment values, the same
#'   length as `trait_values`. Defaults to equidistant indices `1, 2, ..., n`.
#' @return A single numeric value in `(0, 1]`: `1 / (1 + abs(slope))`, where
#'   `slope` is the coefficient from regressing `trait_values` on `env_values`.
#' @examples
#' calculate_PSI(c(2, 4, 6, 8), env_values = c(1, 2, 3, 4))
#' @export
calculate_PSI = function(trait_values, env_values = NULL) {
  # Input validation
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector")
  }

  n_values = length(trait_values)

  # If no environmental values are provided, create equidistant values
  if (is.null(env_values)) {
    env_values = seq_len(n_values)
  }

  # Ensure environmental values are numeric
  if (!is.numeric(env_values)) {
    stop("env_values must be numeric")
  }

  # Check for length mismatch
  if (length(trait_values) != length(env_values)) {
    stop("trait_values and env_values must have the same length")
  }

  # Check for near-constant environmental values
  if (sd(env_values) < .Machine$double.eps) {
    warning("Environmental factor has too little variation, making the PSI calculation unreliable.")
  }

  # Perform a linear regression of trait values on environmental values
  model = lm(trait_values ~ env_values)

  # Extract the regression coefficient (slope)
  beta = unname(coef(model)[2])

  # Calculate the stability score
  stability_score = 1 / (1 + abs(beta))

  return(stability_score)
}

################################

#' Relative stability index (RSI)
#'
#' Computes a stability index based on the coefficient of variation of
#' environment-level trait means, expressed as `1 - CV`.
#'
#' @param trait_values Numeric vector of trait measurements across environments.
#' @param env Optional grouping vector of environment labels, the same length as
#'   `trait_values`. Defaults to equidistant indices `1, 2, ..., n` (each
#'   observation its own environment).
#' @return A single numeric value: `1 - sd(env_means) / mean(env_means)`. Returns
#'   `NA` with a warning when the overall mean of environment means is zero.
#' @examples
#' calculate_RSI(c(2, 4, 6, 8), env = c(1, 1, 2, 2))
#' @export
calculate_RSI <- function(trait_values, env = NULL) {
  # Validate trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }

  # If no environment vector is provided, assume equidistant environments.
  if (is.null(env)) {
    env <- seq_along(trait_values)
  } else {
    if (length(trait_values) != length(env)) {
      stop("trait_values and env must have the same length.")
    }
  }

  # Compute the mean trait value for each unique environment.
  env_means <- tapply(trait_values, env, mean, na.rm = TRUE)

  # Calculate the standard deviation and mean of the environment means.
  overall_sd <- sd(env_means, na.rm = TRUE)
  overall_mean <- mean(env_means, na.rm = TRUE)

  if (overall_mean == 0) {
    warning("The overall mean of environment means is 0; RSI is undefined. Returning NA.")
    return(NA)
  }

  # Compute the Relative Stability Index (RSI)
  RSI_value <- 1 - (overall_sd / overall_mean)

  return(RSI_value)
}

################################

#' Stability index (SI)
#'
#' Computes a stability index as the variance of environment-level trait means
#' relative to their overall mean.
#'
#' @param trait_values Numeric vector of trait measurements across environments.
#' @param env Optional grouping vector of environment labels, the same length as
#'   `trait_values`. Defaults to equidistant indices `1, 2, ..., n` (each
#'   observation its own environment).
#' @return A single numeric value: `var(env_means) / mean(env_means)`. Returns
#'   `NA` with a warning when the overall mean of environment means is zero.
#' @examples
#' calculate_SI(c(2, 4, 6, 8), env = c(1, 1, 2, 2))
#' @export
calculate_SI <- function(trait_values, env = NULL) {
  # Validate trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }

  # If no environment vector is provided, assume equidistant environments.
  if (is.null(env)) {
    env <- seq_along(trait_values)
  } else {
    if (length(trait_values) != length(env)) {
      stop("trait_values and env must have the same length.")
    }
  }

  # Compute the mean trait value for each unique environment.
  env_means <- tapply(trait_values, env, mean, na.rm = TRUE)

  # Calculate the variance among environment means and their overall mean.
  variance_env <- var(env_means, na.rm = TRUE)
  overall_mean <- mean(env_means, na.rm = TRUE)

  if (overall_mean == 0) {
    warning("The overall mean of environment means is 0; SI is undefined. Returning NA.")
    return(NA)
  }

  # Compute the Stability Index (SI)
  SI_value <- variance_env / overall_mean

  return(SI_value)
}

################################

#' Cross-environment covariance and correlation
#'
#' Computes the covariance (and optionally the correlation) between trait values
#' and environmental values.
#'
#' @param trait_values Numeric vector of trait measurements across environments.
#' @param env_values Optional numeric vector of environment values, the same
#'   length as `trait_values`. Defaults to equidistant indices `1, 2, ..., n`.
#' @param return_correlation Logical; if `TRUE`, also compute and return the
#'   Pearson correlation between `trait_values` and `env_values`.
#' @return A list with element `covariance` (and `correlation` when
#'   `return_correlation = TRUE`).
#' @examples
#' cross_env_cov(c(2, 4, 6, 8), env_values = c(1, 2, 3, 4))
#' cross_env_cov(c(2, 4, 6, 8), env_values = c(1, 2, 3, 4), return_correlation = TRUE)
#' @export
cross_env_cov = function(trait_values, env_values = NULL, return_correlation = FALSE) {
  # Input validation
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector")
  }

  n = length(trait_values)

  # If no env_values are provided, assume equidistant environments
  if (is.null(env_values)) {
    env_values = seq_len(n)  # Equidistant environments
  }

  # Ensure correct length
  if (length(env_values) != n) {
    stop("trait_values and env_values must have the same length")
  }

  # Calculate covariance
  covariance = cov(trait_values, env_values)

  # Calculate correlation if requested
  if (return_correlation) {
    correlation = cor(trait_values, env_values)
    return(list(covariance = covariance, correlation = correlation))
  } else {
    return(list(covariance = covariance))
  }
}
