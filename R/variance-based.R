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
