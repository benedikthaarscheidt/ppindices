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
      model = lm(subset_trait_values ~ subset_env_groups)
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
