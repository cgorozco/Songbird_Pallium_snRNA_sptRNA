# Just the cell type proportion function

create_proportion_plots <- function(seuObj, metaUse, class_level, CellClass_numerator, CellClass_denominator, region_pairs, fig_dir = "./Figures", table_dir = "./Tables") {
  
  # Set up metadata
  MetaData2 <- seuObj[[c("orig.ident", "RegionAcronym", metaUse)]]
  
  # Calculate proportions
  CellCounts <- MetaData2 %>%
    group_by(orig.ident, !!sym(class_level), RegionAcronym) %>%
    summarise(Count = n(), .groups = 'drop')
  
  cell_classes_filter <- c(CellClass_numerator, CellClass_denominator) 
  
  proportion_col_name <- sprintf("%s_over__%s", paste(CellClass_numerator, collapse = "_n_"), paste(CellClass_denominator, collapse = "_n_"))
  
  summary_df <- CellCounts %>%
    filter(!!sym(class_level) %in% cell_classes_filter) %>%
    group_by(orig.ident, !!sym(class_level), RegionAcronym) %>%
    summarise(Total_Count = sum(Count), .groups = 'drop') %>%
    spread(key = !!sym(class_level), value = Total_Count, fill = 0) %>%
    mutate(!!sym(proportion_col_name) :=
             rowSums(select(., all_of(CellClass_numerator))) /
             rowSums(select(., all_of(CellClass_denominator))))
  
  # Create lists to store plots and results
  plots <- list()
  results <- list()
  p_values <- data.frame(Comparison = character(), Region_A = character(), RegA_prop_Median = numeric(), Region_B = character(), RegB_prop_Median = numeric(),P_Value = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each combination of Region Acronyms
  for (region_pair in region_pairs) {
    # Filter the data for the current combination
    df_filtered <- summary_df %>% filter(RegionAcronym %in% region_pair)
    df_filtered$RegionAcronym <- factor(df_filtered$RegionAcronym, levels = region_pair)
    
    # Create the formula for the wilcox.test
    formula <- as.formula(paste(proportion_col_name, "~ RegionAcronym"))
    
    # Perform Wilcoxon signed-rank test
    wilcox_test_result <- wilcox.test(formula, data = df_filtered)
    
    # Store the results
    results[[paste(region_pair, collapse = " and ")]] <- wilcox_test_result
    
    # Save the p-values a median
    medians <- df_filtered %>%
      group_by(RegionAcronym) %>%
      summarize(median_value = median(!!sym(proportion_col_name), na.rm = TRUE))
    
    p_values <- rbind(p_values, data.frame(Proportion_numerator = paste(CellClass_numerator, collapse = " & "), Proportion_denominator = paste(CellClass_denominator, collapse = " & "), Region_A = region_pair[1], RegA_prop_Median = medians$median_value[medians$RegionAcronym == region_pair[1]], Region_B = region_pair[2], RegB_prop_Median = medians$median_value[medians$RegionAcronym == region_pair[2]],P_Value = wilcox_test_result$p.value))
    
    # Round numeric values to the nearest fourth decimal
    numeric_cols <- sapply(p_values, is.numeric)
    p_values[numeric_cols] <- lapply(p_values[numeric_cols], function(x) round(x, 4))
    
    
    # Calculate the x-position for the p-value annotation
    x_pos <- mean(as.numeric(factor(region_pair, levels = region_pair)))
    
    # Create the ggplot for the current combination
    plot <- ggplot(df_filtered, aes(x = RegionAcronym, y = !!sym(proportion_col_name))) +
      # stat_summary(fun = median, geom = "bar", fill = c("deepskyblue4", "lightblue"), width = 0.5) +
geom_boxplot(fill = c("#0f0f0f", "#a6a6a6"), alpha = 0.7) +
      geom_jitter(width = 0.2, size = 2) +
      labs(title = paste("Boxplot\n", gsub("over__", "over_\n", proportion_col_name), "\n(", paste(region_pair, collapse = " and "), ")", sep = " "), 
           x = "Region Acronym", y = gsub("over__", "over_\n", proportion_col_name)) +
      theme_classic() +
      theme(text = element_text(size = 10)) +
      # Add the p-value to the plot
      annotate("text", x = x_pos, y = max(df_filtered[[proportion_col_name]], na.rm = TRUE) * 1.05, 
               label = paste("p =", signif(wilcox_test_result$p.value, 3)), vjust = 2, hjust = 0.5, size = 4, color = "black")
    
    # Save the plot
    ggsave(filename = file.path(fig_dir, sprintf("%s__wRegions_%s.tiff", proportion_col_name, paste(region_pair, collapse = "_"))), 
           plot = plot, bg = "white", width = 6, height = 6)
    
    # Add the plot to the list
    plots[[paste(region_pair, collapse = " and ")]] <- plot
  }
  
  # Determine the significance
  p_values$significance <- case_when(
    p_values$P_Value < 0.05 ~ "*",
    p_values$P_Value < 0.07 ~ "~",
    TRUE ~ "No"
  )
  fig_dir
  # Print the table
  print(p_values)
  
  # Save the table as a CSV file
  write_csv(p_values, file.path(table_dir, sprintf("%s__wRegions_%s_pValues.csv", proportion_col_name, paste(region_pair, collapse = "_"))))
  
  return(list(plots = plots, results = results, p_values = p_values))
}






#-- v2
create_proportion_plots <- function(seuObj, metaUse, class_level, CellClass_numerator, CellClass_denominator, region_pairs, fig_dir = "./Figures", table_dir = "./Tables", y_min = NULL, y_max = NULL) {
  
  # Set up metadata
  MetaData2 <- seuObj[[c("orig.ident", "RegionAcronym", metaUse)]]
  
  # Calculate proportions
  CellCounts <- MetaData2 %>%
    group_by(orig.ident, !!sym(class_level), RegionAcronym) %>%
    summarise(Count = n(), .groups = 'drop')
  
  cell_classes_filter <- c(CellClass_numerator, CellClass_denominator) 
  
  proportion_col_name <- sprintf("%s_over__%s", paste(CellClass_numerator, collapse = "_n_"), paste(CellClass_denominator, collapse = "_n_"))
  
  summary_df <- CellCounts %>%
    filter(!!sym(class_level) %in% cell_classes_filter) %>%
    group_by(orig.ident, !!sym(class_level), RegionAcronym) %>%
    summarise(Total_Count = sum(Count), .groups = 'drop') %>%
    spread(key = !!sym(class_level), value = Total_Count, fill = 0) %>%
    mutate(!!sym(proportion_col_name) :=
             rowSums(select(., all_of(CellClass_numerator))) /
             rowSums(select(., all_of(CellClass_denominator))))
  
  # Create lists to store plots and results
  plots <- list()
  results <- list()
  p_values <- data.frame(Comparison = character(), Region_A = character(), RegA_prop_Median = numeric(), Region_B = character(), RegB_prop_Median = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each combination of Region Acronyms
  for (region_pair in region_pairs) {
    # Filter the data for the current combination
    df_filtered <- summary_df %>% filter(RegionAcronym %in% region_pair)
    df_filtered$RegionAcronym <- factor(df_filtered$RegionAcronym, levels = region_pair)
    
    # Create the formula for the wilcox.test
    formula <- as.formula(paste(proportion_col_name, "~ RegionAcronym"))
    
    # Perform Wilcoxon signed-rank test
    wilcox_test_result <- wilcox.test(formula, data = df_filtered)
    
    # Store the results
    results[[paste(region_pair, collapse = " and ")]] <- wilcox_test_result
    
    # Save the p-values and medians
    medians <- df_filtered %>%
      group_by(RegionAcronym) %>%
      summarize(median_value = median(!!sym(proportion_col_name), na.rm = TRUE))
    
    p_values <- rbind(p_values, data.frame(Proportion_numerator = paste(CellClass_numerator, collapse = " & "), Proportion_denominator = paste(CellClass_denominator, collapse = " & "), Region_A = region_pair[1], RegA_prop_Median = medians$median_value[medians$RegionAcronym == region_pair[1]], Region_B = region_pair[2], RegB_prop_Median = medians$median_value[medians$RegionAcronym == region_pair[2]], P_Value = wilcox_test_result$p.value))
    
    # Round numeric values to the nearest fourth decimal
    numeric_cols <- sapply(p_values, is.numeric)
    p_values[numeric_cols] <- lapply(p_values[numeric_cols], function(x) round(x, 4))
    
    # Calculate the x-position for the p-value annotation
    x_pos <- mean(as.numeric(factor(region_pair, levels = region_pair)))
    
    # Create the ggplot for the current combination
    plot <- ggplot(df_filtered, aes(x = RegionAcronym, y = !!sym(proportion_col_name))) +
      geom_boxplot(fill = c("#0f0f0f", "#a6a6a6"), alpha = 0.7) +
      geom_jitter(width = 0.1, size = 2) +
      labs(title = paste("Boxplot\n", gsub("over__", "over_\n", proportion_col_name), "\n(", paste(region_pair, collapse = " and "), ")", sep = " "), 
           x = "Region Acronym", y = gsub("over__", "over_\n", proportion_col_name)) +
      theme_classic() +
      theme(text = element_text(size = 10),
            axis.text = element_text(size = 14, face = "bold", color = "black")) +
      # Add the p-value to the plot
      annotate("text", x = x_pos, y = max(df_filtered[[proportion_col_name]], na.rm = TRUE) * 1.05, 
               label = paste("p =", signif(wilcox_test_result$p.value, 3)), vjust = , hjust = 0.5, size = 4, color = "black")
    
    
    # Set y-axis limits if specified
    if (!is.null(y_min) & !is.null(y_max)) {
      plot <- plot + ylim(y_min, y_max)
    }
    
    # Save the plot
    ggsave(filename = file.path(fig_dir, sprintf("%s__wRegions_%s.tiff", proportion_col_name, paste(region_pair, collapse = "_"))), 
           plot = plot, bg = "white", width = 3.5, height = 4)
    
    # Add the plot to the list
    plots[[paste(region_pair, collapse = " and ")]] <- plot
  }
  
  # Determine the significance
  p_values$significance <- case_when(
    p_values$P_Value < 0.05 ~ "*",
    p_values$P_Value < 0.07 ~ "~",
    TRUE ~ "No"
  )
  
  # Print the table
  print(p_values)
  
  # Save the table as a CSV file
  write_csv(p_values, file.path(table_dir, sprintf("%s__wRegions_%s_pValues.csv", proportion_col_name, paste(region_pair, collapse = "_"))))
  
  return(list(plots = plots, results = results, p_values = p_values))
}



#--- v3 
create_proportion_plots <- function(seuObj, metaUse, class_level_numerator, class_level_denominator, CellClass_numerator, CellClass_denominator, region_pairs, fig_dir = "./Figures", table_dir = "./Tables", y_min = NULL, y_max = NULL) {
  
  # Set up metadata
  MetaData2 <- seuObj[[c("orig.ident", "RegionAcronym", metaUse)]]
  
  # Calculate cell counts for the numerator (using the first class level)
  CellCounts_numerator <- MetaData2 %>%
    group_by(orig.ident, !!sym(class_level_numerator), RegionAcronym) %>%
    summarise(Count = n(), .groups = 'drop')
  
  # Calculate cell counts for the denominator (using the second class level)
  CellCounts_denominator <- MetaData2 %>%
    group_by(orig.ident, !!sym(class_level_denominator), RegionAcronym) %>%
    summarise(Count = n(), .groups = 'drop')
  
  # Filter the relevant cell classes for both numerator and denominator
  cell_classes_filter_numerator <- c(CellClass_numerator)
  cell_classes_filter_denominator <- c(CellClass_denominator)
  
  proportion_col_name <- sprintf("%s_over__%s", paste(CellClass_numerator, collapse = "_n_"), paste(CellClass_denominator, collapse = "_n_"))
  
  # Summarize and calculate the proportion for each combination of orig.ident and RegionAcronym
  summary_df <- CellCounts_numerator %>%
    filter(!!sym(class_level_numerator) %in% cell_classes_filter_numerator) %>%
    group_by(orig.ident, RegionAcronym) %>%
    summarise(Total_Numerator_Count = sum(Count), .groups = 'drop') %>%
    left_join(
      CellCounts_denominator %>%
        filter(!!sym(class_level_denominator) %in% cell_classes_filter_denominator) %>%
        group_by(orig.ident, RegionAcronym) %>%
        summarise(Total_Denominator_Count = sum(Count), .groups = 'drop'),
      by = c("orig.ident", "RegionAcronym")
    ) %>%
    mutate(!!sym(proportion_col_name) := Total_Numerator_Count / Total_Denominator_Count)
  
  # Create lists to store plots and results
  plots <- list()
  results <- list()
  p_values <- data.frame(Comparison = character(), Region_A = character(), RegA_prop_Median = numeric(), Region_B = character(), RegB_prop_Median = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each combination of Region Acronyms
  for (region_pair in region_pairs) {
    # Filter the data for the current combination
    df_filtered <- summary_df %>% filter(RegionAcronym %in% region_pair)
    df_filtered$RegionAcronym <- factor(df_filtered$RegionAcronym, levels = region_pair)
    
    # Create the formula for the wilcox.test
    formula <- as.formula(paste(proportion_col_name, "~ RegionAcronym"))
    
    # Perform Wilcoxon signed-rank test
    wilcox_test_result <- wilcox.test(formula, data = df_filtered)
    
    # Store the results
    results[[paste(region_pair, collapse = " and ")]] <- wilcox_test_result
    
    # Save the p-values and medians
    medians <- df_filtered %>%
      group_by(RegionAcronym) %>%
      summarize(median_value = median(!!sym(proportion_col_name), na.rm = TRUE))
    
    p_values <- rbind(p_values, data.frame(Proportion_numerator = paste(CellClass_numerator, collapse = " & "), Proportion_denominator = paste(CellClass_denominator, collapse = " & "), Region_A = region_pair[1], RegA_prop_Median = medians$median_value[medians$RegionAcronym == region_pair[1]], Region_B = region_pair[2], RegB_prop_Median = medians$median_value[medians$RegionAcronym == region_pair[2]], P_Value = wilcox_test_result$p.value))
    
    # Round numeric values to the nearest fourth decimal
    numeric_cols <- sapply(p_values, is.numeric)
    p_values[numeric_cols] <- lapply(p_values[numeric_cols], function(x) round(x, 4))
    
    # Calculate the x-position for the p-value annotation
    x_pos <- mean(as.numeric(factor(region_pair, levels = region_pair)))
    
    # Create the ggplot for the current combination
    plot <- ggplot(df_filtered, aes(x = RegionAcronym, y = !!sym(proportion_col_name))) +
      geom_boxplot(fill = c("#0f0f0f", "#a6a6a6"), alpha = 0.7) +
      geom_jitter(width = 0.1, size = 2) +
      labs(title = paste("Boxplot\n", gsub("over__", "over_\n", proportion_col_name), "\n(", paste(region_pair, collapse = " and "), ")", sep = " "), 
           x = "Region Acronym", y = gsub("over__", "over_\n", proportion_col_name)) +
      theme_classic() +
      theme(text = element_text(size = 10),
            axis.text = element_text(size = 14, face = "bold", color = "black")) +
      # Add the p-value to the plot
      annotate("text", x = x_pos, y = max(df_filtered[[proportion_col_name]], na.rm = TRUE) * 1.05, 
               label = paste("p =", signif(wilcox_test_result$p.value, 3)), vjust = 0, hjust = 0.5, size = 4, color = "black")
    
    # Set y-axis limits if specified
    if (!is.null(y_min) & !is.null(y_max)) {
      plot <- plot + ylim(y_min, y_max)
    }
    
    # Save the plot
    ggsave(filename = file.path(fig_dir, sprintf("%s__wRegions_%s.tiff", proportion_col_name, paste(region_pair, collapse = "_"))), 
           plot = plot, bg = "white", width = 3.5, height = 4)
    ggsave(filename = file.path(fig_dir, sprintf("%s__wRegions_%s.svg", proportion_col_name, paste(region_pair, collapse = "_"))), 
           plot = plot, bg = "white", width = 3.5, height = 4)
    
    # Add the plot to the list
    plots[[paste(region_pair, collapse = " and ")]] <- plot
  }
  
  # Determine the significance
  p_values$significance <- case_when(
    p_values$P_Value < 0.05 ~ "*",
    p_values$P_Value < 0.07 ~ "~",
    TRUE ~ "No"
  )
  
  # Print the table
  print(p_values)
  
  # Save the table as a CSV file
  write_csv(p_values, file.path(table_dir, sprintf("%s__wRegions_%s_pValues.csv", proportion_col_name, paste(region_pair, collapse = "_"))))
  
  return(list(plots = plots, results = results, p_values = p_values))
}

create_proportion_plots2 <- function(seuObj, metaUse, class_level_numerator, class_level_denominator, CellClass_numerator, CellClass_denominator, region_pairs, fig_dir = "./Figures", table_dir = "./Tables", y_min = NULL, y_max = NULL) {
  
  MetaData2 <- seuObj[[c("orig.ident", "RegionAcronym", metaUse)]]
  
  # Calculate cell counts for the numerator (using the first class level)
  CellCounts_numerator <- MetaData2 %>%
    group_by(orig.ident, !!sym(class_level_numerator), RegionAcronym) %>%
    summarise(Count = n(), .groups = 'drop')
  
  # Calculate cell counts for the denominator (using the second class level)
  CellCounts_denominator <- MetaData2 %>%
    group_by(orig.ident, !!sym(class_level_denominator), RegionAcronym) %>%
    summarise(Count = n(), .groups = 'drop')
  
  # Filter the relevant cell classes for both numerator and denominator
  cell_classes_filter_numerator <- c(CellClass_numerator)
  cell_classes_filter_denominator <- c(CellClass_denominator)
  
  proportion_col_name <- sprintf("%s_over__%s", paste(CellClass_numerator, collapse = "_n_"), paste(CellClass_denominator, collapse = "_n_"))
  
  # Summarize and calculate the proportion for each combination of orig.ident and RegionAcronym
  summary_df <- CellCounts_numerator %>%
    filter(!!sym(class_level_numerator) %in% cell_classes_filter_numerator) %>%
    group_by(orig.ident, RegionAcronym) %>%
    summarise(Total_Numerator_Count = sum(Count), .groups = 'drop') %>%
    left_join(
      CellCounts_denominator %>%
        filter(!!sym(class_level_denominator) %in% cell_classes_filter_denominator) %>%
        group_by(orig.ident, RegionAcronym) %>%
        summarise(Total_Denominator_Count = sum(Count), .groups = 'drop'),
      by = c("orig.ident", "RegionAcronym")
    ) %>%
    mutate(!!sym(proportion_col_name) := Total_Numerator_Count / Total_Denominator_Count)
  
  # Create lists to store plots and results
  plots <- list()
  results <- list()
  p_values <- data.frame(Comparison = character(), Region_A = character(), RegA_prop_Median = numeric(), Region_B = character(), RegB_prop_Median = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each combination of Region Acronyms
  for (region_pair in region_pairs) {
    # Filter the data for the current combination
    df_filtered <- summary_df %>% filter(RegionAcronym %in% region_pair)
    df_filtered$RegionAcronym <- factor(df_filtered$RegionAcronym, levels = region_pair)
    
    # Create the formula for the wilcox.test
    formula <- as.formula(paste(proportion_col_name, "~ RegionAcronym"))
    
    # Perform Wilcoxon signed-rank test
    wilcox_test_result <- wilcox.test(formula, data = df_filtered)
    
    # Store the results
    results[[paste(region_pair, collapse = " and ")]] <- wilcox_test_result
    
    # Save the p-values and medians
    medians <- df_filtered %>%
      group_by(RegionAcronym) %>%
      summarize(median_value = median(!!sym(proportion_col_name), na.rm = TRUE))
    
    p_values <- rbind(p_values, data.frame(Proportion_numerator = paste(CellClass_numerator, collapse = " & "), Proportion_denominator = paste(CellClass_denominator, collapse = " & "), Region_A = region_pair[1], RegA_prop_Median = medians$median_value[medians$RegionAcronym == region_pair[1]], Region_B = region_pair[2], RegB_prop_Median = medians$median_value[medians$RegionAcronym == region_pair[2]], P_Value = wilcox_test_result$p.value))
    
    # Round numeric values to the nearest fourth decimal
    numeric_cols <- sapply(p_values, is.numeric)
    p_values[numeric_cols] <- lapply(p_values[numeric_cols], function(x) round(x, 4))
    
    # Calculate the x- and y- position for the p-value annotation
    x_pos <- mean(as.numeric(factor(region_pair, levels = region_pair)))
    y_line <- max(df_filtered[[proportion_col_name]], na.rm = TRUE) * 1.04
    
    # Determine annotation text based on p-value
    p_value_text <- ifelse(wilcox_test_result$p.value < 0.05, "*", paste("p =", signif(wilcox_test_result$p.value, 3)))
    
    plot <- ggplot(df_filtered, aes(x = RegionAcronym, y = !!sym(proportion_col_name))) +
      geom_boxplot(fill = c("#0f0f0f", "#a6a6a6"), alpha = 0.7, outlier.shape = NA) +  # Exclude outliers
      geom_jitter(width = 0.1, size = 2) +
      geom_segment(aes(x = 1, xend = 2, y = y_line, yend = y_line),
                   color = "#060606", size = 0.8) +
      labs(title = element_blank(), x = element_blank(), y = element_blank()) +
      theme_classic() +
      theme(text = element_text(size = 10),
            axis.text = element_text(size = 14, face = "bold", color = "#060606"),
            axis.line = element_line(color = "#060606")) +
      annotate("text", x = x_pos, y = max(df_filtered[[proportion_col_name]], na.rm = TRUE) * 1.05,
               label = p_value_text, vjust = 0, hjust = 0.5, size = 4, color = "#060606", fontface = "bold") +
      ylim(y_min, NA)

    # Save the plot
    ggsave(filename = file.path(fig_dir, sprintf("%s__wRegions_%s.tiff", proportion_col_name, paste(region_pair, collapse = "_"))), 
           plot = plot, bg = "white", width = 3.5, height = 4)
    ggsave(filename = file.path(fig_dir, sprintf("%s__wRegions_%s.svg", proportion_col_name, paste(region_pair, collapse = "_"))), 
           plot = plot, bg = "white", width = 3.5, height = 4)
    
    # Add the plot to the list
    plots[[paste(region_pair, collapse = " and ")]] <- plot
  }
  
  # Determine the significance
  p_values$significance <- case_when(
    p_values$P_Value < 0.05 ~ "*",
    p_values$P_Value < 0.07 ~ "~",
    TRUE ~ "No"
  )
  
  # Print the table
  print(p_values)
  
  # Save the table as a CSV file
  write_csv(p_values, file.path(table_dir, sprintf("%s__wRegions_%s_pValues.csv", proportion_col_name, paste(region_pair, collapse = "_"))))
  
  return(list(plots = plots, results = results, p_values = p_values))
}







####====== For Anova
create_proportion_plots_anova <- function(seuObj, metaUse, class_level_numerator, class_level_denominator, CellClass_numerator, CellClass_denominator, regions_to_compare, fig_dir = "./Figures", table_dir = "./Tables", y_min = NULL, y_max = NULL) {
  
  # Set up metadata
  MetaData2 <- seuObj[[c("orig.ident", "RegionAcronym", metaUse)]]
  
  # Calculate cell counts for the numerator (using the first class level)
  CellCounts_numerator <- MetaData2 %>%
    group_by(orig.ident, !!sym(class_level_numerator), RegionAcronym) %>%
    summarise(Count = n(), .groups = 'drop')
  
  # Calculate cell counts for the denominator (using the second class level)
  CellCounts_denominator <- MetaData2 %>%
    group_by(orig.ident, !!sym(class_level_denominator), RegionAcronym) %>%
    summarise(Count = n(), .groups = 'drop')
  
  # Filter the relevant cell classes for both numerator and denominator
  cell_classes_filter_numerator <- c(CellClass_numerator)
  cell_classes_filter_denominator <- c(CellClass_denominator)
  
  proportion_col_name <- sprintf("%s_over__%s", paste(CellClass_numerator, collapse = "_n_"), paste(CellClass_denominator, collapse = "_n_"))
  
  # Summarize and calculate the proportion for each combination of orig.ident and RegionAcronym
  summary_df <- CellCounts_numerator %>%
    filter(!!sym(class_level_numerator) %in% cell_classes_filter_numerator) %>%
    group_by(orig.ident, RegionAcronym) %>%
    summarise(Total_Numerator_Count = sum(Count), .groups = 'drop') %>%
    left_join(
      CellCounts_denominator %>%
        filter(!!sym(class_level_denominator) %in% cell_classes_filter_denominator) %>%
        group_by(orig.ident, RegionAcronym) %>%
        summarise(Total_Denominator_Count = sum(Count), .groups = 'drop'),
      by = c("orig.ident", "RegionAcronym")
    ) %>%
    mutate(!!sym(proportion_col_name) := Total_Numerator_Count / Total_Denominator_Count)
  
  # Filter for the selected regions to compare
  df_filtered <- summary_df %>% filter(RegionAcronym %in% regions_to_compare)
  df_filtered$RegionAcronym <- factor(df_filtered$RegionAcronym, levels = regions_to_compare)
  
  # Perform ANOVA test
  anova_test <- aov(as.formula(paste(proportion_col_name, "~ RegionAcronym")), data = df_filtered)
  anova_summary <- summary(anova_test)
  
  # Post-hoc Tukey HSD test to identify pairwise differences if significant
  tukey_test <- TukeyHSD(anova_test)
  
  # Create the plot
  plot <- ggplot(df_filtered, aes(x = factor(RegionAcronym, levels = c("Hyper", "Meso", "Av",  "An", "LMAN", "NC", "NCM", "HVC", "Arco", "RA")), 
                                  y = !!sym(proportion_col_name))) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.1, size = 2) +
    labs(title = paste("Boxplot\n", gsub("over__", "over_\n", proportion_col_name), "\n(", paste(regions_to_compare, collapse = ", "), ")", sep = " "), 
         x = "Region Acronym", y = gsub("over__", "over_\n", proportion_col_name)) +
    theme_classic() +
    theme(text = element_text(size = 10),
          axis.text = element_text(size = 14, face = "bold", color = "black"))
  
# Set y-axis limits if specified
  if (!is.null(y_min) & !is.null(y_max)) {
    plot <- plot + ylim(y_min, y_max)
  }
  
  # Save the plot
  ggsave(filename = file.path(fig_dir, sprintf("%s__wRegions_%s.tiff", proportion_col_name, paste(regions_to_compare, collapse = "_"))), 
         plot = plot, bg = NA, width = 7.5, height = 5)
  
  # Create the results table for the ANOVA and Tukey test
  anova_results <- data.frame(ANOVA_Summary = capture.output(anova_summary))
  tukey_results <- as.data.frame(tukey_test$RegionAcronym)
  colnames(tukey_results) <- c("Diff", "Lower_CI", "Upper_CI", "P_Value")
  
  # Print the summary tables
  print(anova_results)
  print(tukey_results)
  
  # Save the tables as CSV files
  write_csv(anova_results, file.path(table_dir, sprintf("%s__ANOVA_Results.csv", proportion_col_name)))
  write_csv(tukey_results, file.path(table_dir, sprintf("%s__TukeyHSD_Results.csv", proportion_col_name)))
  
  return(list(plot = plot, anova_summary = anova_summary, tukey_test = tukey_test, anova_results = anova_results, tukey_results = tukey_results))
}

