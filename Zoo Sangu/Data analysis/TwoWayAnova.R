TwoWayAnova <- function(data) {
  # Load necessary packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("emmeans", quietly = TRUE)) install.packages("emmeans")
  if (!requireNamespace("multcomp", quietly = TRUE)) install.packages("multcomp")
  
  library(ggplot2)
  library(emmeans)
  library(multcomp)
  
  # Check if data has the correct structure
  if (ncol(data) != 3) {
    stop("Data should have exactly three columns: two factors and one numeric variable.")
  }
  
  # Define column names for easy access
  factor_var1 <- names(data)[1]
  factor_var2 <- names(data)[2]
  numeric_var <- names(data)[3]
  
  # Ensure the factor columns are factors
  data[[factor_var1]] <- as.factor(data[[factor_var1]])
  data[[factor_var2]] <- as.factor(data[[factor_var2]])
  
  # Check that the numeric variable is numeric
  if (!is.numeric(data[[numeric_var]])) {
    stop(paste(numeric_var, "must be a numeric variable."))
  }
  
  # Run two-way ANOVA
  aov_model <- aov(as.formula(paste(numeric_var, "~", factor_var1, "*", factor_var2)), data = data)
  aov_summary <- summary(aov_model) # Summary of the ANOVA model
  
  # Check for significant interaction and perform Tukey's HSD if applicable
  TukeySHD <- NULL
  if (aov_summary[[1]][["Pr(>F)"]][3] < 0.05) { # Check interaction significance
    TukeySHD <- TukeyHSD(aov_model, which = interaction(data[[factor_var1]], data[[factor_var2]]))
  }
  
  # Plotting results
  plot <- ggplot(data, aes_string(x = factor_var1, y = numeric_var, fill = factor_var2)) +
    geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.7) +
    geom_jitter(size = 2, position = position_jitterdodge(jitter.width = 0.2), color = "blue", alpha = 0.6) +
    theme_bw() +
    labs(x = factor_var1, y = numeric_var, title = "Two-Way ANOVA Results") +
    theme(text = element_text(size = 12), plot.title = element_text(hjust = 0.5))
  
  # Return all relevant information as a list
  return(list(
    ANOVA_Summary = aov_summary,
    TukeyHSD = TukeySHD,
    Plot = plot
  ))
}


# Example dataset
data_example <- iris[, c("Species", "Sepal.Width", "Sepal.Length")]
data_example <- data_example %>% mutate(Sepal.Length = Sepal.Length + rnorm(nrow(data_example), 0, 0.5))

# Running the TwoWayAnova function
result <- TwoWayAnova(data_example)

# Checking the results
print(result$ANOVA_Summary)
if (!is.null(result$TukeyHSD)) {
  print(result$TukeyHSD)
}
print(result$Plot)










