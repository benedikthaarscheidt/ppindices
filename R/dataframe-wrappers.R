rdpi_mean_calculation = function(dataframe, trait_cols, sp = NULL, factors = NULL, factors_not_in_dataframe = NULL, stat_analysis = NULL) {

  # Convert column indices to names if necessary
  if (!is.null(sp)) {
    sp = if (is.numeric(sp)) names(dataframe)[sp] else sp
  }
  traits = if (is.numeric(trait_cols)) names(dataframe)[trait_cols] else trait_cols

  # Combine internal and external factors into a single dataframe
  if (!is.null(factors_not_in_dataframe)) {
    if (length(factors_not_in_dataframe[[1]]) != nrow(dataframe)) {
      stop("The length of external factors must match the number of rows in the dataframe.")
    }
    external_factors_df = as.data.frame(factors_not_in_dataframe)
    if (!is.null(factors)) {
      factors = if (is.numeric(factors)) names(dataframe)[factors] else factors
      combined_factors_df = cbind(dataframe[factors], external_factors_df)
    } else {
      combined_factors_df = external_factors_df
    }
    dataframe$Combined_Factors = interaction(combined_factors_df, drop = TRUE)
  } else if (!is.null(factors)) {
    factors = if (is.numeric(factors)) names(dataframe)[factors] else factors
    dataframe$Combined_Factors = interaction(dataframe[factors], drop = TRUE)
  } else {
    stop("You must provide either internal factors, external factors, or both.")
  }

  dataframe$Combined_Factors = as.factor(dataframe$Combined_Factors)

  all_rdpi_results = data.frame()  # Initialize a dataframe to store all RDPI values

  if (is.null(sp)) {
    unique_species = list("Single_Group" = dataframe)
  } else {
    unique_species = split(dataframe, dataframe[[sp]])
  }

  for (species_name in names(unique_species)) {
    species_data = unique_species[[species_name]]

    for (trait in traits) {
      RDPI_values = data.frame(sp = character(), env1 = character(), env2 = character(), rdpi = numeric())

      # Get unique environment combinations
      env_levels = levels(species_data$Combined_Factors)
      n_env_levels = length(env_levels)

      mean_values = aggregate(species_data[[trait]], list(species_data$Combined_Factors), mean)
      colnames(mean_values) = c("Env_Combination", "Mean_Trait")

      for (i in 1:(n_env_levels - 1)) {
        for (j in (i + 1):n_env_levels) {
          mean_i = mean_values$Mean_Trait[mean_values$Env_Combination == env_levels[i]]
          mean_j = mean_values$Mean_Trait[mean_values$Env_Combination == env_levels[j]]

          # Calculate RDPI for this trait between the two environments
          rdpi_value = abs(mean_i - mean_j) / min(mean_i, mean_j)

          # Append the RDPI value for this species, trait, and environment pair
          RDPI_values = rbind(RDPI_values, data.frame(sp = species_name, env1 = env_levels[i], env2 = env_levels[j], rdpi = rdpi_value))
        }
      }

      all_rdpi_results = rbind(all_rdpi_results, RDPI_values)
    }
  }

  # If statistical analysis is requested, perform ANOVA and Tukey's HSD
  if (!is.null(stat_analysis)) {
    all_trait_data = data.frame()

    for (species_name in names(unique_species)) {
      species_data = unique_species[[species_name]]

      for (trait in traits) {
        trait_data = data.frame(
          sp = species_name,
          trait = trait,
          Combined_Factors = species_data$Combined_Factors,
          Trait_Value = species_data[[trait]]
        )

        all_trait_data = rbind(all_trait_data, trait_data)
      }
    }

    # Perform ANOVA and Tukey's HSD test
    anova_results = list()
    tukey_results = list()

    for (trait in traits) {
      fit = aov(Trait_Value ~ Combined_Factors, data = subset(all_trait_data, trait == trait))
      anova_results[[trait]] = summary(fit)

      # Perform Tukey's HSD test
      Tukey = agricolae::HSD.test(fit, trt = "Combined_Factors")
      tukey_results[[trait]] = Tukey
    }

    # Create boxplots for the traits across environmental combinations
    boxplot_traits = ggplot2::ggplot(all_trait_data, ggplot2::aes(x = Combined_Factors, y = Trait_Value, fill = trait)) +
      ggplot2::geom_boxplot() +
      ggplot2::facet_wrap(~trait, scales = "free_y") +
      ggplot2::theme_bw() +
      ggplot2::xlab("Environmental Combinations") +
      ggplot2::ylab("Trait Values") +
      ggplot2::ggtitle("Boxplots of Trait Values Across Environmental Combinations") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    return(list(
      rdpi_results = all_rdpi_results,
      trait_boxplots = boxplot_traits,
      anova_results = anova_results,
      tukey_results = tukey_results
    ))
  }

  # If stat_analysis is NULL, return only the RDPI values
  return(all_rdpi_results)
}

################################

evwpi_calculation = function(dataframe, trait_cols, sp = NULL, env_col, factors, use_median = FALSE) {
  # NOTE: this was an early idea and needs further development

  # Convert traits and sp from column indices to names if necessary
  if (!is.null(sp) && is.numeric(sp)) {
    sp = names(dataframe)[sp]
  }
  traits = if (is.numeric(trait_cols)) names(dataframe)[trait_cols] else trait_cols

  # Handle environments: if it's a vector of indices, extract the columns; if it's a list of vectors, use directly
  if (is.numeric(env_col)) {
    environments = dataframe[[names(dataframe)[env_col]]]
  }

  # Handle factors: if factors is a vector of indices/names, extract the columns; if it's a list of vectors, use directly
  if (is.numeric(factors) || is.character(factors)) {
    factors = lapply(factors, function(f) dataframe[[f]])
  }

  # Ensure the factors are correctly handled and match the length of the dataframe
  if (any(sapply(factors, length) != nrow(dataframe))) {
    stop("Length of factors must match the number of rows in the dataframe.")
  }

  # Step 1: Calculate the variance induced by each factor
  variance_scores = sapply(factors, function(factor) {
    var(factor, na.rm = TRUE)
  })

  # Normalize the variance scores to get the weights
  total_variance = sum(variance_scores, na.rm = TRUE)
  weights = variance_scores / total_variance

  all_results = list()

  if (is.null(sp)) {
    unique_species = list("Single_Group" = dataframe)
  } else {
    unique_species = split(dataframe, dataframe[[sp]])
  }

  # Loop over each species or group
  for (species_name in names(unique_species)) {
    species_data = unique_species[[species_name]]

    EVWPI_results = list()

    # Initialize a data frame to store EVWPI values for each pair of environments
    EVWPI_values = data.frame(sp = character(), trait = character(), evwpi = numeric())
    unique_envs=unique(environments)

    # Loop over each trait
    for (trait in traits) {
      for (i in 1:(length(unique_envs) - 1)) {
        for (j in (i + 1):length(unique_envs)) {

          env_i = unique_envs[i]
          env_j = unique_envs[j]

          trait_i = species_data[[trait]][environments == env_i]
          trait_j = species_data[[trait]][environments == env_j]
          abs_diff = abs(outer(trait_i, trait_j, "-"))

          # Calculate mean or median absolute phenotypic distance
          mean_abs_diff = if (use_median) median(abs_diff, na.rm = TRUE) else mean(abs_diff, na.rm = TRUE)

          # Calculate the weighted environmental distance using the variance-based weights
          weighted_env_distance = sum(weights)


          # Calculate EVWPI for this pair of environments
          evwpi_value = if (weighted_env_distance != 0) mean_abs_diff / weighted_env_distance else NaN


          # Store the results
          EVWPI_values = rbind(EVWPI_values, data.frame(sp = species_name, trait = trait, evwpi = evwpi_value))
        }
      }
    }

    all_results[[species_name]] = EVWPI_values
  }

  return(all_results)
}

################################

calculate_SPI = function(data, env_col, trait_cols, env1, env2, reference_env) {

  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }

  # Initialize a list to store SPI values for each trait
  spi_results = numeric(length(trait_cols))
  names(spi_results) = trait_cols

  i=0
  # Loop through each trait column and calculate SPI
  for (trait_col in trait_cols) {
    i=i+1
    # Extract trait data
    trait_values = if (is.numeric(trait_col)) data[[trait_col]] else data[[trait_col]]

    # Calculate mean trait values for env1 and env2
    trait_env1_value = mean(trait_values[env_col == env1], na.rm = TRUE)
    trait_env2_value = mean(trait_values[env_col == env2], na.rm = TRUE)

    # Calculate the standard deviation of the reference environment's trait values
    sd_reference = sd(trait_values[env_col == reference_env], na.rm = TRUE)

    # Calculate the Standardized Plasticity Index (SPI) for this trait
    spi = (trait_env2_value - trait_env1_value) / sd_reference

    # Store the result
    spi_results[i] = spi
  }

  return(spi_results)
}

################################

calculate_SPM = function(data, env_col, trait_col, resident_env, nonresident_env) {

  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }

  # Extract trait data
  trait_values = data[[trait_col]]

  # Calculate mean trait values for resident and nonresident environments
  trait_resident = mean(trait_values[env_col == resident_env], na.rm = TRUE)
  trait_nonresident = mean(trait_values[env_col == nonresident_env], na.rm = TRUE)

  # Calculate the Standardized Plasticity Metric (SPM)
  spm = abs(trait_resident - trait_nonresident) / trait_resident

  return(spm)
}

################################

calculate_Plasticity_Ratio = function(data, env_col, trait_cols, pop_col) {

  # Handle env_col
  if (is.numeric(env_col) && length(env_col) == 1) {
    env_col = data[[env_col]]
  } else if (is.vector(env_col) && length(env_col) == nrow(data)) {
    env_col = env_col
  } else {
    env_col = data[[env_col]]
  }

  # Handle pop_col
  if (is.numeric(pop_col) && length(pop_col) == 1) {
    pop_col = data[[pop_col]]
  } else if (is.vector(pop_col) && length(pop_col) == nrow(data)) {
    pop_col = pop_col
  } else {
    pop_col = data[[pop_col]]
  }

  # Initialize a vector to store Plasticity Ratios for each trait
  plasticity_ratios = numeric(length(trait_cols))
  names(plasticity_ratios) = trait_cols

  # Loop over each trait column
  for (i in seq_along(trait_cols)) {
    # Extract the trait values for the current trait
    trait_values = if (is.numeric(trait_cols[i])) data[[trait_cols[i]]] else data[[trait_cols[i]]]

    # Create a temporary data frame for the ANOVA
    temp_data = data.frame(trait_values, pop_col, env_col)

    # Perform a one-way ANOVA to obtain SSpop and SStotal
    anova_result = aov(trait_values ~ pop_col + env_col, data = temp_data)

    # Extract the sum of squares between populations (SSpop) and total sum of squares (SStotal)
    ss_total = sum(anova(anova_result)[, "Sum Sq"])
    ss_pop = anova(anova_result)["pop_col", "Sum Sq"]

    # Calculate and store the Plasticity Ratio
    plasticity_ratios[i] = ss_pop / ss_total
  }

  return(plasticity_ratios)
}
