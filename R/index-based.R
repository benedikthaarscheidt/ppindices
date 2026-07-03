#' Grand plasticity index from adjusted environment means
#'
#' Fits a linear model of trait values on environment (optionally adjusting for
#' a covariate), extracts adjusted (least-squares) means per environment via
#' `emmeans`, and computes the coefficient of variation of those means, with an
#' optional exclusion of a designated control environment.
#'
#' @param trait_values Numeric vector of trait measurements.
#' @param env_data Vector of environment labels, the same length as
#'   `trait_values`; coerced to a factor.
#' @param covariate_data Optional numeric covariate, the same length as
#'   `trait_values`, added to the model as an additive term.
#' @param control_env Optional environment label to exclude from the treatment
#'   means before computing the index.
#' @return A single numeric value: `sd(treatment_means) / mean(treatment_means)`.
#' @examples
#' trait_values <- c(2, 2.5, 3, 4, 4.5, 5, 6, 6.5, 7)
#' env_data <- factor(rep(c("E1", "E2", "E3"), each = 3))
#' calculate_grand_plasticity(trait_values, env_data)
#' @export
calculate_grand_plasticity = function(trait_values, env_data, covariate_data = NULL, control_env = NULL) {

  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("The 'emmeans' package is required but not installed.")
  }

  env_data = factor(env_data)

  if (!is.null(covariate_data)) {
    model = lm(trait_values ~ covariate_data + env_data)
  } else {
    model = lm(trait_values ~ env_data)
  }
  adjusted_means = emmeans::emmeans(model, ~ env_data)
  adjusted_means_summary = summary(adjusted_means)

  if (!is.null(control_env)) {
    treatment_means = adjusted_means_summary$emmean[adjusted_means_summary$env_data != control_env]
  } else {
    treatment_means = adjusted_means_summary$emmean
  }

  if (length(treatment_means) == 0) {
    stop("Could not extract treatment means. Check your input data and control environment.")
  }

  sd_means = sd(treatment_means, na.rm = TRUE)
  grand_mean = mean(treatment_means, na.rm = TRUE)

  grand_plasticity = sd_means / grand_mean

  return(grand_plasticity)
}

################################

#' Phenotypic Plasticity Factor (PPF) 
#'
#' For each pair of environment groups (all pairs by default), fits a linear
#' model and computes the percentage difference between the least-squares means
#' of the pair, then averages across pairs.
#'
#' @param trait_values Numeric vector of trait measurements.
#' @param env_groups Vector of environment group labels, the same length as
#'   `trait_values`; coerced to a factor.
#' @param covariate_values Optional numeric vector or data frame of covariates
#'   to include additively in the model.
#' @param env_pairs Optional list of environment-label pairs to evaluate.
#'   Defaults to all pairwise combinations of `levels(env_groups)`.
#' @return A single numeric value: the mean, across environment pairs, of
#'   `100 * abs(LSM1 - LSM2) / LSM1`.
#' @note The example below is illustrative only and is wrapped in `\dontrun{}`.
#'   The function fits each pairwise model on locally renamed variables
#'   (`subset_trait_values`, `subset_env_groups`) but then calls
#'   `emmeans::emmeans(model, ~ env_groups)`, which looks for a term named
#'   `env_groups` that the model does not contain, so it errors with
#'   "No variable named env_groups in the reference grid".
#' @examples
#' \dontrun{
#' trait_values <- c(2, 2.5, 3, 4, 4.5, 5)
#' env_groups <- factor(rep(c("A", "B"), each = 3))
#' calculate_PPF(trait_values, env_groups)
#' }
#' @export
calculate_PPF = function(trait_values, env_groups, covariate_values = NULL, env_pairs = NULL) {
  # Ensure the environment groups are treated as factors

  env_groups = as.factor(env_groups)

  # Default: Calculate PPF for all possible combinations of environment pairs
  if (is.null(env_pairs)) {
    env_pairs = combn(levels(env_groups), 2, simplify = FALSE)
  } else if (is.vector(env_pairs) && length(env_pairs) == 2) {
    env_pairs = list(env_pairs)  # Convert a single pair into a list of one pair
  }

  # Initialize a vector to store PPF values for each environment pair
  PPF_values = numeric(length(env_pairs))
  names(PPF_values) = sapply(env_pairs, function(pair) paste(pair, collapse = "-"))

  # Create a progress bar
  pb <- txtProgressBar(min = 0, max = length(env_pairs), style = 3)

  # Loop through each environment pair
  for (i in seq_along(env_pairs)) {
    setTxtProgressBar(pb, i)
    env_pair = env_pairs[[i]]

    # Subset the data to include only the current environment pair
    subset_idx = which(env_groups %in% env_pair)
    subset_trait_values = trait_values[subset_idx]
    subset_env_groups = env_groups[subset_idx]

    if (!is.null(covariate_values)) {
      # Handle single or multiple covariates
      if (is.vector(covariate_values)) {
        subset_covariates = covariate_values[subset_idx]
        data = data.frame(trait_values = subset_trait_values, env_groups = subset_env_groups, Covariate = subset_covariates)
      } else if (is.data.frame(covariate_values)) {
        subset_covariates = covariate_values[subset_idx, , drop = FALSE]
        data = data.frame(trait_values = subset_trait_values, env_groups = subset_env_groups, subset_covariates)
      }
      formula = as.formula("trait_values ~ env_groups + .")
      model = lm(formula, data = data)
    } else {
      data = data.frame(trait_values = subset_trait_values,
                        env_groups = subset_env_groups)
      model = lm(trait_values ~ env_groups, data = data)
    }

    # Calculate least square means (LSMs) for the current environment pair
    lsm = as.data.frame(emmeans::emmeans(model, ~ env_groups))

    # Calculate PPF for the current pair: 100 * (|LSM1 - LSM2| / LSM1)
    PPF_values[i] = 100 * (abs(lsm[1, "emmean"] - lsm[2, "emmean"]) / lsm[1, "emmean"])
  }

  close(pb)
  return(sum(PPF_values) / length(PPF_values))
}

################################

#' Phenotypic plasticity index (Pi)
#'
#' Computes the classic `(max - min) / max` phenotypic plasticity index, either
#' for a single numeric vector or for each of several trait columns in a data
#' frame or matrix.
#'
#' @param data A numeric vector of trait values, or a data frame / matrix
#'   containing one or more trait columns.
#' @param trait_cols Optional character or numeric vector identifying trait
#'   columns in `data`. When `NULL`, `data` is treated as a single numeric
#'   vector.
#' @return A single numeric value when `trait_cols` is `NULL`, or a named
#'   numeric vector (one value per trait) otherwise.
#' @examples
#' calculate_Phenotypic_Plasticity_Index(c(2, 4, 6, 8))
#' df <- data.frame(trait1 = c(2, 4, 6, 8), trait2 = c(1, 3, 5, 7))
#' calculate_Phenotypic_Plasticity_Index(df, trait_cols = c("trait1", "trait2"))
#' @export
calculate_Phenotypic_Plasticity_Index = function(data, trait_cols = NULL) {
  if (is.null(trait_cols)) {
    # If trait_col is not provided, assume data is a numeric vector
    max_value = max(data, na.rm = TRUE)
    min_value = min(data, na.rm = TRUE)

    # Calculate Pi for the single vector
    Pi = (max_value - min_value) / max_value
    return(Pi)
  } else {
    # If trait_cols are provided, assume data is a data frame or matrix
    Pi_values = numeric(length(trait_cols))  # Initialize vector to store Pi for each trait
    names(Pi_values) = trait_cols  # Name the vector by trait columns

    for (i in seq_along(trait_cols)) {
      trait_column = trait_cols[i]
      trait_data = if (is.numeric(trait_column)) data[[trait_column]] else data[[trait_column]]

      max_value = max(trait_data, na.rm = TRUE)
      min_value = min(trait_data, na.rm = TRUE)

      # Calculate Pi as (Max - Min) / Max
      Pi_values[i] = (max_value - min_value) / max_value
    }

    return(Pi_values)
  }
}

################################

#' Plasticity quotient (PQ)
#'
#' Computes the range of trait values standardized by the range of environment
#' values.
#'
#' @param trait_values Numeric vector of trait measurements across environments.
#' @param env_values Optional numeric vector of environment values, the same
#'   length as `trait_values`. Defaults to equidistant indices `1, 2, ..., n`.
#' @return A single numeric value: `range(trait_values) / range(env_values)`.
#'   Returns `NA` with a warning when the environment range is zero.
#' @examples
#' calculate_PQ(c(2, 4, 6, 8), env_values = c(1, 2, 3, 4))
#' @export
calculate_PQ = function(trait_values, env_values = NULL) {
  # Ensure trait_values is numeric
  if (!is.numeric(trait_values)) stop("trait_values must be a numeric vector.")

  n = length(trait_values)

  # If no explicit env_values are given, assume equidistant environments
  if (is.null(env_values)) {
    env_values = seq_len(n)
  }

  # Ensure env_values is numeric and of correct length
  if (!is.numeric(env_values)) stop("env_values must be a numeric vector.")
  if (length(env_values) != n) stop("trait_values and env_values must have the same length.")

  # Compute the range of trait values
  range_trait = abs(max(trait_values) - min(trait_values))

  # Compute the range of environment values (used for standardization)
  range_env = abs(max(env_values) - min(env_values))

  # Ensure env range is non-zero to avoid division errors
  if (range_env == 0) {
    warning("Environment range is zero. PQ calculation is not meaningful.")
    return(NA)
  }

  # Compute the Plasticity Quotient (PQ)
  PQ_value = range_trait / range_env

  return(PQ_value)
}

################################

#' Phenotypic Range (PR)
#'
#' Computes the range (max - min) of trait values, either across all
#' environments at once or separately within each unique environment.
#'
#' @param trait_values Numeric vector of trait measurements across environments.
#' @param env_values Optional vector of environment labels, the same length as
#'   `trait_values`. Defaults to equidistant indices `1, 2, ..., n`.
#' @param across Logical; if `TRUE` (default), compute a single range across all
#'   environments. If `FALSE`, compute the range within each unique environment.
#' @return A single numeric value when `across = TRUE`, or a numeric vector (one
#'   value per unique environment) when `across = FALSE`.
#' @examples
#' calculate_PR(c(2, 4, 6, 8), env_values = c(1, 2, 3, 4))
#' calculate_PR(c(2, 4, 6, 8), env_values = c(1, 1, 2, 2), across = FALSE)
#' @export
calculate_PR = function(trait_values, env_values = NULL, across = TRUE) {
  # Ensure trait_values is numeric
  if (!is.numeric(trait_values)) stop("trait_values must be a numeric vector.")

  n = length(trait_values)

  # If no explicit env_values are given, assume equidistant environments
  if (is.null(env_values)) {
    env_values = seq_len(n)
  }

  # Ensure env_values is numeric and matches trait_values in length
  if (!is.numeric(env_values)) stop("env_values must be a numeric vector.")
  if (length(env_values) != n) stop("trait_values and env_values must have the same length.")

  # If across = TRUE, compute PR across all environments
  if (across) {
    PR_value = max(trait_values, na.rm = TRUE) - min(trait_values, na.rm = TRUE)
    return(PR_value)
  }

  # Compute PR within each environment
  unique_envs = unique(env_values)
  PR_values = numeric(length(unique_envs))

  for (i in seq_along(unique_envs)) {
    env_mask = env_values == unique_envs[i]
    env_trait_values = trait_values[env_mask]
    PR_values[i] = max(env_trait_values, na.rm = TRUE) - min(env_trait_values, na.rm = TRUE)
  }

  return(PR_values)
}

################################

#' Norm-of-reaction width (NRW)
#'
#' Computes the range (max - min) of trait values, at the level of raw
#' observations, environment means, or per-group environment means.
#'
#' @param trait_values Numeric vector of trait measurements.
#' @param env_values Optional vector of environment labels, the same length as
#'   `trait_values`. Defaults to equidistant indices `1, 2, ..., n`.
#' @param group_values Optional vector of group labels (e.g. genotype), the
#'   same length as `trait_values`, to compute NRW separately per group.
#' @param across Logical; if `TRUE`, compute the range directly across raw
#'   `trait_values` (ignoring environment grouping). Defaults to `FALSE`.
#' @return A single numeric value (overall or per environment-mean range), or a
#'   named numeric vector of per-group values when `group_values` is supplied.
#' @examples
#' calculate_NRW(c(2, 4, 6, 8), env_values = c(1, 1, 2, 2))
#' @export
calculate_NRW = function(trait_values, env_values = NULL, group_values = NULL, across = FALSE) {
  if (!is.numeric(trait_values)) stop("trait_values must be numeric.")

  # Assume equidistant environments if not provided
  if (is.null(env_values)) env_values = seq_along(trait_values)
  if (length(env_values) != length(trait_values)) stop("trait_values and env_values must be the same length.")

  # Compute NRW across all environments
  if (across) {
    return(max(trait_values, na.rm = TRUE) - min(trait_values, na.rm = TRUE))
  }

  # Compute NRW per group (e.g., genotype)
  if (!is.null(group_values)) {
    unique_groups = unique(group_values)
    NRW_values = setNames(numeric(length(unique_groups)), unique_groups)

    for (g in unique_groups) {
      mask = group_values == g
      mean_per_env = tapply(trait_values[mask], env_values[mask], mean, na.rm = TRUE)
      NRW_values[g] = max(mean_per_env, na.rm = TRUE) - min(mean_per_env, na.rm = TRUE)
    }

    return(NRW_values)
  }

  # Default: Calculate NRW across environments
  mean_per_env = tapply(trait_values, env_values, mean, na.rm = TRUE)
  NRW_value = max(mean_per_env, na.rm = TRUE) - min(mean_per_env, na.rm = TRUE)

  return(NRW_value)
}

################################

#' Environmental Sensitivity Performance (ESP)
#'
#' Computes, for each environment, the relative deviation of its mean trait
#' value from the overall mean, and sums the absolute deviations across
#' environments.
#'
#' @param trait_values Numeric vector of trait measurements.
#' @param env_values Optional vector of environment labels, the same length as
#'   `trait_values`. Defaults to equidistant indices `1, 2, ..., n`.
#' @param env_subset Optional vector restricting the calculation to a subset of
#'   environment labels.
#' @return A single numeric value: the sum of absolute per-environment relative
#'   deviations from the overall mean.
#' @examples
#' calculate_ESP(c(2, 4, 6, 8), env_values = c(1, 1, 2, 2))
#' @export
calculate_ESP = function(trait_values, env_values = NULL, env_subset = NULL) {
  if (!is.numeric(trait_values)) stop("trait_values must be numeric.")

  # Assume equidistant environments if not provided
  if (is.null(env_values)) env_values = seq_along(trait_values)
  if (length(env_values) != length(trait_values)) stop("trait_values and env_values must be the same length.")

  # Use all unique environments unless specific ones are given
  unique_envs = unique(env_values)
  if (!is.null(env_subset)) unique_envs = unique_envs[unique_envs %in% env_subset]

  # Compute mean trait value across all environments
  mean_trait_all = mean(trait_values, na.rm = TRUE)

  # Compute ESP per environment
  ESP_values = setNames(numeric(length(unique_envs)), unique_envs)

  for (e in unique_envs) {
    mean_trait_env = mean(trait_values[env_values == e], na.rm = TRUE)
    ESP_values[as.character(e)] = (mean_trait_env - mean_trait_all) / mean_trait_all
  }

  y=abs(ESP_values)
  ESP_grand=sum(y)
  return(ESP_grand)
}

################################

#' Relative plasticity index (RPI)
#'
#' Computes the relative distance between trait values recorded in different
#' environments, either for a single specified environment pair or averaged
#' across all pairwise combinations of environments.
#'
#' @param trait_values Numeric vector of trait measurements.
#' @param env_values Optional vector of environment labels, the same length as
#'   `trait_values`. Defaults to equidistant indices `1, 2, ..., n`.
#' @param env1 Optional label of the first environment in a specific pair to
#'   compare.
#' @param env2 Optional label of the second environment in a specific pair to
#'   compare.
#' @return A single numeric value: the mean relative distance
#'   `abs(x1 - x2) / (x1 + x2)`, for the requested pair or averaged over all
#'   pairs. Returns `NA` with a warning when fewer than two unique environments
#'   are present.
#' @examples
#' calculate_RPI(c(2, 3, 6, 7), env_values = c(1, 1, 2, 2))
#' @export
calculate_RPI = function(trait_values, env_values = NULL, env1 = NULL, env2 = NULL) {
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
    warning("Only one unique environment present. RPI cannot be calculated.")
    return(NA)
  }

  # If specific environments are given, calculate RPI for that pair
  if (!is.null(env1) && !is.null(env2)) {
    if (!(env1 %in% unique_envs) || !(env2 %in% unique_envs)) stop("env1 and env2 must be in env_values.")

    idx1 = which(env_values == env1)
    idx2 = which(env_values == env2)

    # Ensure equal sample size
    min_len = min(length(idx1), length(idx2))
    idx1 = idx1[seq_len(min_len)]
    idx2 = idx2[seq_len(min_len)]

    RPI_values = abs(trait_values[idx1] - trait_values[idx2]) / (trait_values[idx1] + trait_values[idx2])
    return(mean(RPI_values, na.rm = TRUE))
  }

  # Otherwise, compute RPI for all pairs
  env_combinations = combn(unique_envs, 2, simplify = TRUE)
  RPI_results = numeric(ncol(env_combinations))

  for (i in seq_len(ncol(env_combinations))) {
    env1_current = env_combinations[1, i]
    env2_current = env_combinations[2, i]

    idx1 = which(env_values == env1_current)
    idx2 = which(env_values == env2_current)

    min_len = min(length(idx1), length(idx2))
    idx1 = idx1[seq_len(min_len)]
    idx2 = idx2[seq_len(min_len)]

    RPI_values = abs(trait_values[idx1] - trait_values[idx2]) / (trait_values[idx1] + trait_values[idx2])
    RPI_results[i] = mean(RPI_values, na.rm = TRUE)
  }
  RPI_grand=sum(RPI_results)/length(RPI_results)
  names(RPI_results) = apply(env_combinations, 2, paste, collapse = "-")
  return(RPI_grand)
}

################################

#' Phenotypic Flexibility Index (PFI)
#'
#' Identifies the sample with the largest absolute deviation from a baseline
#' and expresses that deviation relative to its own baseline.
#'
#' @param trait_values Numeric vector of trait measurements.
#' @param baseline_values Optional numeric vector of baseline values, the same
#'   length as `trait_values`. Defaults to the overall mean of `trait_values`
#'   repeated for every sample.
#' @return A single numeric value: `max(abs(trait - baseline)) / baseline` at
#'   the sample of maximum deviation. Returns `NA` with a warning when that
#'   baseline value is zero.
#' @examples
#' calculate_PFI(c(2, 4, 6, 8))
#' @export
calculate_PFI <- function(trait_values, baseline_values = NULL) {
  # Validate input for trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }

  n <- length(trait_values)

  # If baseline_values is not provided, use the overall mean for all samples.
  if (is.null(baseline_values)) {
    baseline_values <- rep(mean(trait_values, na.rm = TRUE), n)
  } else {
    # Validate baseline_values
    if (!is.numeric(baseline_values)) {
      stop("baseline_values must be a numeric vector.")
    }
    if (length(baseline_values) != n) {
      stop("baseline_values must have the same length as trait_values.")
    }
  }

  # Calculate the absolute deviations from the baseline for each sample
  deviations <- abs(trait_values - baseline_values)

  # Identify the sample with the maximum deviation
  max_deviation_index <- which.max(deviations)
  max_deviation <- deviations[max_deviation_index]

  # Retrieve the baseline corresponding to the maximum deviation sample
  baseline_for_max <- baseline_values[max_deviation_index]

  if (baseline_for_max == 0) {
    warning("The baseline value corresponding to the maximum deviation is 0; PFI is undefined. Returning NA.")
    return(NA)
  }

  # Calculate the PFI as the ratio of the maximum deviation to its corresponding baseline
  pfi_value <- max_deviation / baseline_for_max

  return(pfi_value)
}

################################

#' Absolute Plasticity Coefficient (APC)
#'
#' Orders observations by environment, computes per-environment means, and
#' averages the absolute differences between consecutive environment means.
#'
#' @param trait_values Numeric vector of trait measurements.
#' @param env_labels Optional vector of environment labels, the same length as
#'   `trait_values`. Defaults to equidistant indices `1, 2, ..., n`.
#' @param sequential_env Logical; if `TRUE` (default), order environments
#'   numerically/alphabetically (or by level for ordered factors). If `FALSE`,
#'   preserve the order of first appearance.
#' @return A single numeric value: the mean absolute difference between
#'   consecutive environment means. Returns `NA` with a warning when fewer than
#'   two unique environments are present.
#' @examples
#' calculate_APC(c(2, 4, 6, 8), env_labels = c(1, 1, 2, 2))
#' @export
calculate_APC <- function(trait_values, env_labels = NULL, sequential_env = TRUE) {
  # Input validation for trait_values
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }

  # If no environment labels are provided, assume equidistant environments.
  if (is.null(env_labels)) {
    env_labels <- seq_along(trait_values)
  } else {
    if (length(trait_values) != length(env_labels)) {
      stop("trait_values and env_labels must have the same length.")
    }
  }

  # Order the data based on the environment labels.
  if (sequential_env) {
    if (is.numeric(env_labels)) {
      ordering <- order(env_labels)
    } else if (is.factor(env_labels)) {
      if (is.ordered(env_labels)) {
        ordering <- order(as.numeric(env_labels))
      } else {
        ordering <- seq_along(env_labels)
      }
    } else {  # For character vectors
      ordering <- order(env_labels)
    }
  } else {
    # Use the order of first appearance
    uniq_env <- unique(env_labels)
    ordering <- order(match(env_labels, uniq_env))
  }

  trait_values <- trait_values[ordering]
  env_labels <- env_labels[ordering]

  # Calculate mean trait value for each unique environment
  env_means <- tapply(trait_values, env_labels, mean, na.rm = TRUE)
  if (length(env_means) < 2) {
    warning("Less than two unique environments found; APC is undefined. Returning NA.")
    return(NA)
  }

  # Compute absolute differences between consecutive environment means
  abs_diff <- abs(diff(env_means))

  # Calculate APC as the mean of these absolute differences
  apc_value <- mean(abs_diff, na.rm = TRUE)

  return(apc_value)
}

################################

#' Multivariate plasticity index (MVPi)
#'
#' Runs a PCA on a matrix of trait values, retains the first `n_axes` principal
#' components, computes per-environment centroids in that reduced space, and
#' summarizes plasticity as the mean pairwise Euclidean distance between
#' centroids.
#'
#' @param trait_data A numeric matrix (or object coercible to one) of trait
#'   values, with samples in rows and traits in columns.
#' @param env Optional vector of environment labels, the same length as
#'   `nrow(trait_data)`. Defaults to equidistant indices `1, 2, ..., n`.
#' @param n_axes Number of leading principal component axes to retain. Defaults
#'   to `1`.
#' @return A single numeric value: the mean pairwise Euclidean distance between
#'   environment centroids in PCA space. Returns `NA` with a warning when fewer
#'   than two environments are present.
#' @examples
#' trait_data <- matrix(c(1, 2, 3, 2, 3, 4, 5, 6, 5, 6, 7, 6), ncol = 3, byrow = TRUE)
#' env <- c("E1", "E1", "E2", "E2")
#' calculate_MVPi(trait_data, env = env, n_axes = 2)
#' @export
calculate_MVPi <- function(trait_data, env = NULL, n_axes = 1) {
  # Convert trait_data to a numeric matrix if not already one.
  if (!is.matrix(trait_data)) {
    trait_data <- as.matrix(trait_data)
  }
  if (!is.numeric(trait_data)) {
    stop("trait_data must be numeric.")
  }

  # Determine number of samples.
  n_samples <- nrow(trait_data)

  # If no environment vector is provided, assume equidistant environments.
  if (is.null(env)) {
    env <- seq_len(n_samples)
  } else {
    if (length(env) != n_samples) {
      stop("The length of env must equal the number of rows in trait_data.")
    }
  }

  # Convert env to a factor.
  env <- factor(env)

  # Perform PCA on the trait data.
  pca_result <- prcomp(trait_data, center = TRUE, scale. = FALSE)

  # Ensure that n_axes does not exceed the available dimensions.
  if (n_axes > ncol(pca_result$x)) {
    stop("n_axes exceeds the number of available principal components.")
  }

  # Retain the first n_axes of the PCA scores.
  scores <- pca_result$x[, 1:n_axes, drop = FALSE]

  # Compute centroids for each environment.
  env_levels <- levels(env)
  if (length(env_levels) < 2) {
    warning("Less than two environments found; MVPi is undefined. Returning NA.")
    return(NA)
  }
  centroids <- sapply(env_levels, function(x) {
    colMeans(scores[env == x, , drop = FALSE], na.rm = TRUE)
  })
  centroids <- t(centroids)  # rows = environment centroids

  # Compute pairwise Euclidean distances between centroids.
  distances <- as.vector(dist(centroids))

  # Calculate MVPi as the mean Euclidean distance.
  mvpi <- mean(distances, na.rm = TRUE)

  return(mvpi)
}

################################

#' Plasticity response index (PRI)
#'
#' Computes the difference between mean trait values under an "extreme"
#' condition and a "control" condition, relative to the overall mean.
#'
#' @param trait_values Numeric vector of trait measurements.
#' @param env_indicator Binary vector (`0` = control, `1` = extreme), the same
#'   length as `trait_values`.
#' @return A single numeric value:
#'   `(mean_extreme - mean_control) / overall_mean`. Returns `NA` with a
#'   warning when the overall mean is zero.
#' @examples
#' calculate_PRI(c(2, 4, 6, 8), env_indicator = c(0, 0, 1, 1))
#' @export
calculate_PRI <- function(trait_values, env_indicator) {
  # Validate input
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector.")
  }
  if (length(trait_values) != length(env_indicator)) {
    stop("trait_values and env_indicator must have the same length.")
  }
  if (!all(env_indicator %in% c(0, 1))) {
    stop("env_indicator must be a binary vector with values 0 (control) and 1 (extreme).")
  }

  # Calculate overall mean
  overall_mean <- mean(trait_values, na.rm = TRUE)
  if (overall_mean == 0) {
    warning("Overall mean of trait values is 0; PRI is undefined. Returning NA.")
    return(NA)
  }

  # Calculate means for extreme and control environments
  mean_extreme <- mean(trait_values[env_indicator == 1], na.rm = TRUE)
  mean_control <- mean(trait_values[env_indicator == 0], na.rm = TRUE)

  # Calculate PRI
  pri_value <- (mean_extreme - mean_control) / overall_mean

  return(pri_value)
}
