#' Relative distance plasticity index (RDPI)
#'
#' Computes, for every pair of environments, the relative distance between
#' trait values (`abs(x2 - x1) / (x1 + x2)`), and averages across all pairs.
#'
#' @param trait_values Numeric vector of trait measurements across environments.
#' @param env_values Optional vector of environment values, the same length as
#'   `trait_values`. Defaults to equidistant indices `1, 2, ..., n`.
#' @return A single numeric value: the mean relative distance across all
#'   pairwise environment comparisons.
#' @examples
#' calculate_rdpi(c(2, 4, 6, 8), env_values = c(1, 2, 3, 4))
#' @export
calculate_rdpi = function(trait_values, env_values = NULL) {
  # Input validation
  if (!is.numeric(trait_values)) {
    stop("trait_values must be numeric")
  }

  n_envs = length(trait_values)

  # If no environment values provided, create sequential environments
  if (is.null(env_values)) {
    env_values = seq_len(n_envs)
  }

  # Ensure same length
  if (length(trait_values) != length(env_values)) {
    stop("trait_values and env_values must have the same length")
  }

  # Get all pairs of environment indices
  env_pairs = combn(n_envs, 2)
  n_pairs = ncol(env_pairs)

  # Calculate RDPIs for each pair
  rdpis = numeric(n_pairs)

  for (i in seq_len(n_pairs)) {
    idx1 = env_pairs[1, i]
    idx2 = env_pairs[2, i]

    # Calculate relative distance for this pair
    abs_diff = abs(trait_values[idx2] - trait_values[idx1])
    sum_vals = trait_values[idx1] + trait_values[idx2]

    # Handle potential division by zero
    if (sum_vals == 0) {
      rdpis[i] = 0
    } else {
      rdpis[i] = abs_diff / sum_vals
    }
  }

  # Calculate final RDPI as mean of all RDPIs
  rdpi = mean(rdpis, na.rm = TRUE)

  return(rdpi)
}

################################

#' Environment-standardized plasticity index (ESPI)
#'
#' Computes the range of environment-level trait means, standardized by the
#' range of the environment values. Accepts a single trait (numeric vector) or
#' multiple traits (columns of a data frame).
#'
#' @param trait_values A numeric vector of trait measurements, or a data frame
#'   where each column is a trait.
#' @param env_values Optional vector of environment values. For the
#'   single-trait (vector) form it must have the same length as `trait_values`.
#'   Defaults to equidistant indices `1, 2, ..., n`.
#' @return A single numeric value for a vector input, or a named numeric vector
#'   (one value per trait column) for a data frame input.
#' @examples
#' calculate_ESPI(c(2, 3, 6, 7), env_values = c(1, 1, 2, 2))
#' @export
calculate_ESPI = function(trait_values, env_values = NULL) {

  if (is.null(env_values)) {
    env_values = seq_along(trait_values)

  }

  if (!is.vector(trait_values) && !is.data.frame(trait_values)) {
    stop("trait_values must be a numeric vector or a data frame where each column is a trait.")
  }

  if (length(env_values) != length(trait_values)) {
    stop("env_values must be the same length as trait_values.")
  }

  # Function to calculate ESPI for a single trait
  calculate_single_espi = function(single_trait) {
    means = tapply(single_trait, env_values, mean, na.rm = TRUE)
    max_mean = max(means, na.rm = TRUE)
    min_mean = min(means, na.rm = TRUE)

    abs_env_distance = abs(max(as.numeric(env_values), na.rm = TRUE) - min(as.numeric(env_values), na.rm = TRUE))
    if (abs_env_distance > 0) {
      return((max_mean - min_mean) / abs_env_distance)
    } else {
      return(NA)
    }
  }

  # Handle single trait
  if (is.vector(trait_values)) {
    return(calculate_single_espi(trait_values))
  } else if (is.data.frame(trait_values)) {
    # Handle multiple traits
    espi_results = sapply(trait_values, calculate_single_espi)
    names(espi_results) = colnames(trait_values)
    return(espi_results)
  }
}

################################

#' Environment-standardized plasticity index with individual distances (ESPIID)
#'
#' For every pair of (equidistant) environments, computes the mean (or median)
#' absolute pairwise difference between individual trait values in the two
#' environments, standardized by the distance between environment indices, and
#' averages across all pairs.
#'
#' @param trait_values Numeric vector of trait measurements across equidistant
#'   environments (each element is one observation, in environment order).
#' @param use_median Logical; if `TRUE`, use the median absolute pairwise
#'   difference instead of the mean. Defaults to `FALSE`.
#' @return A single numeric value: the mean ESPIID across all environment
#'   pairs.
#' @examples
#' calculate_espiid(c(2, 4, 6, 8))
#' @export
calculate_espiid = function(trait_values, use_median = FALSE) {

  env_values = seq_along(trait_values) # Create equidistant environment indices

  # Input validation
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }

  if (length(trait_values) != length(env_values)) {
    stop("trait_values and env_values must have the same length.")
  }

  env_values = as.factor(env_values) # Convert to factor

  # Get all unique pairs of environments
  env_levels = levels(env_values)
  n_envs = length(env_levels)

  if (n_envs < 2) {
    stop("At least two environments are required to calculate ESPIID.")
  }

  # Initialize storage for ESPIID values
  espiid_values = numeric()

  # Loop over all environment pairs
  for (i in 1:(n_envs - 1)) {
    for (j in (i + 1):n_envs) {
      # Extract trait values for the two environments
      trait_i = trait_values[env_values == env_levels[i]]
      trait_j = trait_values[env_values == env_levels[j]]

      # Calculate absolute differences between all pairs of individuals
      abs_diff = abs(outer(trait_i, trait_j, "-"))

      # Calculate mean or median absolute phenotypic distance
      mean_abs_diff = if (use_median) median(abs_diff, na.rm = TRUE) else mean(abs_diff, na.rm = TRUE)

      # Calculate absolute distance between environmental values
      abs_env_distance = abs(as.numeric(env_levels[i]) - as.numeric(env_levels[j]))

      # Handle potential division by zero
      espiid_value = if (abs_env_distance > 0) mean_abs_diff / abs_env_distance else NA

      # Store ESPIID value for this environment pair
      espiid_values = c(espiid_values, espiid_value)
    }
  }

  # Return the mean ESPIID across all environment pairs
  return(mean(espiid_values, na.rm = TRUE))
}

################################

#' Median-based phenotypic plasticity index (PImd)
#'
#' Computes environment-level medians of a trait and expresses their range
#' relative to the maximum median.
#'
#' @param trait_values Numeric vector of trait measurements.
#' @param env Vector of environment labels, the same length as `trait_values`;
#'   coerced to a factor.
#' @return A single numeric value: `(max(medians) - min(medians)) / max(medians)`.
#'   Returns `NA` with a warning when the maximum median is zero.
#' @examples
#' calculate_PImd(c(2, 4, 6, 8), env = factor(c(1, 1, 2, 2)))
#' @export
calculate_PImd = function(trait_values, env) {
  # Ensure env is treated as a factor
  if (!is.factor(env)) {
    env = factor(env)
  }

  # Get the levels of the environment
  levels = levels(env)
  medians = numeric(length(levels))

  # Calculate medians for each environment level
  for (i in seq_along(levels)) {
    level = levels[i]
    medians[i] = median(trait_values[env == level], na.rm = TRUE)
  }

  # Check if maximum median is zero to avoid division by zero
  if (max(medians) == 0) {
    warning("Maximum median is zero; PImd cannot be calculated.")
    return(NA)
  }

  # Calculate PImd as (Max - Min) / Max
  PImd = (max(medians) - min(medians)) / max(medians)

  return(PImd)
}

################################

#' Least-squares-mean phenotypic plasticity index (PILSM)
#'
#' Fits a linear model of trait values on environment (optionally adjusting for
#' covariates), extracts least-squares means per environment via `emmeans`, and
#' expresses their range relative to the maximum least-squares mean.
#'
#' @param trait_values Numeric vector of trait measurements.
#' @param env Vector of environment labels, the same length as `trait_values`;
#'   coerced to a factor.
#' @param covariates Optional numeric vector or data frame of covariates to
#'   include additively in the model.
#' @return A single numeric value: `(max_lsm - min_lsm) / max_lsm`.
#' @examples
#' trait_values <- c(2, 2.5, 3, 4, 4.5, 5)
#' env <- factor(rep(c("E1", "E2"), each = 3))
#' calculate_PILSM(trait_values, env)
#' @export
calculate_PILSM = function(trait_values, env, covariates = NULL) {
  # Ensure env is treated as a factor
  if (!is.factor(env)) {
    env = factor(env)
  }

  # Fit the linear model
  if (is.null(covariates)) {
    model = lm(trait_values ~ env)
  } else {
    # If covariates are a vector, convert to a data frame
    if (is.vector(covariates)) {
      covariates = data.frame(Covariate = covariates)
    }
    # Combine env and covariates into a single data frame
    data = data.frame(trait_values = trait_values, env = env, covariates)
    formula = as.formula("trait_values ~ env + .")
    model = lm(formula, data = data)
  }

  # Calculate least square means (LSMs) for each environment
  lsm = as.data.frame(emmeans::emmeans(model, ~ env))

  # Calculate PILSM
  max_lsm = max(lsm$emmean, na.rm = TRUE)
  min_lsm = min(lsm$emmean, na.rm = TRUE)
  PILSM = (max_lsm - min_lsm) / max_lsm

  return(PILSM)
}

################################

#' General phenotypic distance (PD)
#'
#' Computes a distance-based plasticity measure between environments, either
#' from an explicit control/stress grouping vector, or via one of three
#' methods: pairwise mean absolute differences across all environment pairs,
#' differences relative to the lowest-mean reference environment, or the raw
#' trait range.
#'
#' @param trait_values Numeric vector of trait measurements.
#' @param env_values Optional vector of environment labels, the same length as
#'   `trait_values`. Defaults to equidistant indices `1, 2, ..., n`.
#' @param control_stress_vector Optional vector, the same length as
#'   `trait_values`, labelling each observation as `"Control"`/`"Stress"` or
#'   `0`/`1`. When supplied, `method` is ignored.
#' @param method One of `"pairwise"` (default), `"reference"`, or
#'   `"variability"`; see Description.
#' @return A single numeric value summarizing phenotypic distance.
#' @examples
#' calculate_general_PD(c(2, 4, 6, 8), env_values = c(1, 1, 2, 2))
#' @export
calculate_general_PD = function(trait_values, env_values = NULL, control_stress_vector = NULL, method = "pairwise") {
  if (!is.numeric(trait_values)) stop("trait_values must be numeric.")

  num_values = length(trait_values)
  if (num_values < 2) stop("At least two trait values are required.")

  # If no env_values are provided, assume equidistant environments
  if (is.null(env_values)) {
    env_values = seq_len(num_values)
  }

  if (length(env_values) != num_values) stop("trait_values and env_values must have the same length.")

  unique_envs = unique(env_values)
  num_envs = length(unique_envs)

  if (num_envs < 2) stop("At least two distinct environments are required.")

  # If a control vs. stress vector is provided, use it directly
  if (!is.null(control_stress_vector)) {
    if (length(control_stress_vector) != num_values) {
      stop("control_stress_vector must have the same length as trait_values.")
    }

    # Automatically detect whether the vector is categorical or numeric
    if (all(control_stress_vector %in% c("Control", "Stress"))) {
      control_values = trait_values[control_stress_vector == "Control"]
      stress_values = trait_values[control_stress_vector == "Stress"]
    } else if (all(control_stress_vector %in% c(0, 1))) {
      control_values = trait_values[control_stress_vector == 0]
      stress_values = trait_values[control_stress_vector == 1]
    } else {
      stop("control_stress_vector must contain 'Control'/'Stress' labels or numeric 0/1 values.")
    }

    if (length(control_values) != length(stress_values)) {
      stop("Control and stress groups must have the same number of values.")
    }

    # Compute PD as mean absolute difference
    return(mean(abs(stress_values - control_values), na.rm = TRUE))
  }

  # If no control vs. stress is provided, use the selected method
  if (method == "pairwise") {
    pd_values = c()
    env_combinations = combn(unique_envs, 2)
    for (i in 1:ncol(env_combinations)) {
      env1 = env_combinations[1, i]
      env2 = env_combinations[2, i]

      values1 = trait_values[env_values == env1]
      values2 = trait_values[env_values == env2]

      if (length(values1) != length(values2)) next  # Skip if unequal sample sizes

      pd_values = c(pd_values, mean(abs(values1 - values2), na.rm = TRUE))
    }
    return(mean(pd_values, na.rm = TRUE))

  } else if (method == "reference") {
    env_means = tapply(trait_values, env_values, mean, na.rm = TRUE)
    reference_env = names(which.min(env_means))

    pd_values = c()
    for (env in unique_envs) {
      if (env == reference_env) next
      values1 = trait_values[env_values == reference_env]
      values2 = trait_values[env_values == env]
      if (length(values1) != length(values2)) next  # Skip if unequal sample sizes
      pd_values = c(pd_values, mean(abs(values1 - values2), na.rm = TRUE))
    }
    return(mean(pd_values, na.rm = TRUE))

  } else if (method == "variability") {
    return(max(trait_values, na.rm = TRUE) - min(trait_values, na.rm = TRUE))

  } else {
    stop("Invalid method. Choose 'pairwise', 'reference', or 'variability'.")
  }
}

################################

#' Fitness-relative plasticity index (FPI)
#'
#' Computes the relative change in trait value between environments (or
#' between control and stress conditions), expressed as a fraction of the
#' baseline (control / first) value, averaged across environment pairs.
#'
#' @param trait_values Numeric vector of trait measurements.
#' @param env_values Optional vector of environment labels, the same length as
#'   `trait_values`. Defaults to equidistant indices `1, 2, ..., n`.
#' @param control_stress Optional vector, the same length as `trait_values`,
#'   labelling each observation as `0`/`1` or `"Control"`/`"Stress"`. When
#'   supplied, `env_values` pairing is bypassed in favor of a direct
#'   control-vs-stress comparison.
#' @return A single numeric value: the mean fractional change in trait value.
#'   Returns `NA` with a warning when only one unique environment is present.
#' @examples
#' calculate_FPI(c(2, 4, 6, 8), env_values = c(1, 1, 2, 2))
#' @export
calculate_FPI = function(trait_values, env_values = NULL, control_stress = NULL) {
  # Ensure trait_values is numeric
  if (!is.numeric(trait_values)) stop("trait_values must be a numeric vector.")

  n = length(trait_values)

  # Assume equidistant environments if not provided
  if (is.null(env_values)) {
    env_values = seq_len(n)
  }

  # Ensure env_values has correct length
  if (length(env_values) != n) stop("env_values must have the same length as trait_values.")

  unique_envs = unique(env_values)

  # If only one environment is present, return NA
  if (length(unique_envs) < 2) {
    warning("Only one unique environment present. FPI cannot be calculated.")
    return(NA)
  }

  # If control-stress mapping is provided
  if (!is.null(control_stress)) {
    if (length(control_stress) != n) stop("control_stress must have the same length as trait_values.")

    # Ensure binary format (0/1 or Control/Stress)
    if (!all(control_stress %in% c(0, 1, "Control", "Stress"))) {
      stop("control_stress must be a vector of 0/1 or 'Control'/'Stress'.")
    }

    # Convert to 0/1 if needed
    if (is.character(control_stress)) {
      control_stress = ifelse(control_stress == "Control", 0, 1)
    }

    control_values = trait_values[control_stress == 0]
    stress_values = trait_values[control_stress == 1]

    # Ensure equal sample size
    min_len = min(length(control_values), length(stress_values))
    control_values = control_values[seq_len(min_len)]
    stress_values = stress_values[seq_len(min_len)]

    # Compute FPI
    FPI_values = (stress_values - control_values) / control_values
    return(mean(FPI_values, na.rm = TRUE))
  }

  # Otherwise, compute FPI for all environment pairs
  env_combinations = combn(unique_envs, 2, simplify = TRUE)
  FPI_results = numeric(ncol(env_combinations))

  for (i in seq_len(ncol(env_combinations))) {
    env1 = env_combinations[1, i]
    env2 = env_combinations[2, i]

    idx1 = which(env_values == env1)
    idx2 = which(env_values == env2)

    min_len = min(length(idx1), length(idx2))
    idx1 = idx1[seq_len(min_len)]
    idx2 = idx2[seq_len(min_len)]

    FPI_values = (trait_values[idx2] - trait_values[idx1]) / trait_values[idx1]
    FPI_results[i] = mean(FPI_values, na.rm = TRUE)
  }

  FPI_grand=mean(FPI_results)
  names(FPI_results) = apply(env_combinations, 2, paste, collapse = "-")
  return(FPI_grand)
}

################################

#' Transplant phenotypic shift (TPS)
#'
#' Computes the relative change in trait value between a native environment and
#' a transplanted environment, expressed as a fraction of the native value.
#'
#' @param trait_values Numeric vector of trait measurements.
#' @param env_values Vector of environment labels, the same length as
#'   `trait_values`.
#' @param native_env Label (in `env_values`) of the native environment.
#' @param transplanted_env Label (in `env_values`) of the transplanted
#'   environment.
#' @return A single numeric value: the mean of
#'   `(transplanted - native) / native`.
#' @examples
#' calculate_TPS(c(2, 4, 6, 8), env_values = c(1, 1, 2, 2), native_env = 1, transplanted_env = 2)
#' @export
calculate_TPS = function(trait_values, env_values, native_env, transplanted_env) {
  # Ensure trait_values is numeric
  if (!is.numeric(trait_values)) stop("trait_values must be a numeric vector.")

  # Ensure env_values is provided and has correct length
  if (is.null(env_values)) stop("env_values must be provided.")
  if (length(env_values) != length(trait_values)) stop("env_values must have the same length as trait_values.")

  # Check if specified native and transplanted environments exist in env_values
  unique_envs = unique(env_values)
  if (!(native_env %in% unique_envs) || !(transplanted_env %in% unique_envs)) {
    stop("native_env and transplanted_env must be in env_values.")
  }

  # Extract trait values for native and transplanted environments
  native_data = trait_values[env_values == native_env]
  transplanted_data = trait_values[env_values == transplanted_env]

  # Ensure equal sample size
  min_len = min(length(native_data), length(transplanted_data))
  native_data = native_data[seq_len(min_len)]
  transplanted_data = transplanted_data[seq_len(min_len)]

  # Compute TPS
  TPS_values = (transplanted_data - native_data) / native_data
  return(mean(TPS_values, na.rm = TRUE))
}

################################

#' Daily plasticity index (DPI)
#'
#' Computes the rate of change in trait value between two time points, divided
#' by the elapsed time interval.
#'
#' @param trait_values_time1 Numeric vector of trait measurements at time 1.
#' @param trait_values_time2 Numeric vector of trait measurements at time 2,
#'   the same length as `trait_values_time1`.
#' @param time_interval Optional numeric vector of positive time intervals, the
#'   same length as `trait_values_time1`. Defaults to a vector of ones.
#' @return A numeric vector: `(trait_values_time2 - trait_values_time1) / time_interval`.
#' @examples
#' calculate_DPI(c(2, 4, 6, 8), c(3, 5, 7, 9))
#' @export
calculate_DPI = function(trait_values_time1, trait_values_time2, time_interval = NULL) {
  # Ensure trait values are numeric vectors
  if (!is.numeric(trait_values_time1) || !is.numeric(trait_values_time2)) {
    stop("trait_values_time1 and trait_values_time2 must be numeric vectors.")
  }

  # If no time_interval vector is provided, assume a vector of ones
  if (is.null(time_interval)) {
    time_interval = rep(1, length(trait_values_time1))
  } else if (!is.numeric(time_interval)) {
    stop("time_interval must be a numeric vector.")
  }

  # Ensure all vectors have the same length
  if (length(trait_values_time1) != length(trait_values_time2) || length(trait_values_time1) != length(time_interval)) {
    stop("trait_values_time1, trait_values_time2, and time_interval must all have the same length.")
  }

  # Ensure all elements of time_interval are positive numeric values
  if (any(time_interval <= 0)) {
    stop("All elements of time_interval must be positive numeric values.")
  }

  # Calculate DPI for each sample
  dpi = (trait_values_time2 - trait_values_time1) / time_interval

  return(dpi)
}
