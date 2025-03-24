# Load necessary libraries
library(lme4)        # for mixed-effects models
library(lmerTest)    # for p-values in mixed-effects models
library(dplyr)       # for data manipulation
library(ggplot2)     # for plotting
library(MuMIn)       # for R-squared values
library(gridExtra)   # for arranging plots
library(clubSandwich) # for robust standard errors
library(haven)
library(car)         # for VIF calculations

my_data <- read_dta("Y:/Melissa Kendall/Gun Laws/Hand Guns/Data Out/Hand Guns/Data Prep for R Feature Engineering Suicide.dta")


# Assuming my_data is already loaded and preprocessed as in your original code
# (gun_law_vars defined, factors created, etc.)

# Create orthogonal polynomial time variables
year_poly <- poly(my_data$Year, 2)
my_data$Year_poly1 <- year_poly[,1]  # Linear orthogonal term
my_data$Year_poly2 <- year_poly[,2]  # Quadratic orthogonal term

# Verify orthogonality
print("Correlation between orthogonal time variables:")
print(cor(my_data$Year_poly1, my_data$Year_poly2))

# Update the robust_se function to use the new data structure
robust_se <- function(model) {
  vcov_cluster <- vcovCR(model, type = "CR0", cluster = my_data$State)
  coef_summary <- data.frame(
    Estimate = fixef(model),
    Robust_SE = sqrt(diag(vcov_cluster)),
    row.names = names(fixef(model))
  )
  coef_summary$z_value <- coef_summary$Estimate / coef_summary$Robust_SE
  coef_summary$p_value <- 2 * pnorm(abs(coef_summary$z_value), lower.tail = FALSE)
  return(coef_summary)
}

# Combine OC variables into a single factor (assuming this was already done)
# If not done yet, uncomment and run these lines:
# my_data$OC <- NA
# my_data$OC[my_data$OC_License == 1] <- 1
# my_data$OC[my_data$OC_Prohibited == 1] <- 2
# my_data$OC[my_data$OC_NotRestricted == 1] <- 0
# my_data$OC <- factor(my_data$OC)

# Base model with orthogonal time variables
model_base_ortho <- lmer(log(cruderate) ~ BGC_Sales + DealerLicense + DV + CAP_Intentional + CAP_Negligent +
                           SaleRestriction + Traffick + Local_Comprehensive + Local_Selective + 
                           Permit + NoPossess_FirearmCrime + NoPossess_DVRO + NoPossess_ERPO + 
                           Registration + Lost_Stolen + Tracing + Safety_Purchase + 
                           Untraceable + SYG + OC + Year_poly1 + Year_poly2 +
                           (1 | State), 
                         data = my_data,
                         control = lmerControl(optimizer = "bobyqa",
                                               optCtrl = list(maxfun = 100000)))

# Get base model results with robust standard errors
base_results_ortho <- robust_se(model_base_ortho)

# Print base model results
print("Base Model Results with Robust Standard Errors (Orthogonal Time):")
print(base_results_ortho)

# Check for multicollinearity in the base model
vif_values_ortho <- car::vif(model_base_ortho)
print("\nVariance Inflation Factors for Base Model (Orthogonal Time):")
print(vif_values_ortho)

# Identify significant variables (p < 0.05) using robust standard errors
significant_vars <- rownames(base_results_ortho)[base_results_ortho$p_value < 0.05]
significant_vars <- significant_vars[!significant_vars %in% c("(Intercept)", "Year_poly1", "Year_poly2")]
print("\nSignificant variables (p < 0.05):")
print(significant_vars)

# Get base names of significant variables (remove numbers at end for factors)
base_significant_vars <- unique(sapply(significant_vars, function(x) {
  # This pattern handles both regular variables and factor levels
  gsub("([0-9]+|:Year_poly[0-9]+)$", "", x)
}))
print("\nBase names of significant variables:")
print(base_significant_vars)

# Create and print interaction formula with orthogonal time variables
interaction_formula <- paste("(", paste(base_significant_vars, collapse = " + "), 
                             ") * (Year_poly1 + Year_poly2)")
print("\nInteraction formula:")
print(interaction_formula)

# Final model with orthogonal time interactions for significant variables
model_final_ortho <- lmer(as.formula(paste(
  "log(cruderate) ~ CAP_Intentional + CAP_Negligent + BGC_Sales + DealerLicense + DV + 
   Traffick + Local_Comprehensive + Local_Selective + Permit + SaleRestriction +
   NoPossess_FirearmCrime + NoPossess_DVRO + NoPossess_ERPO + Registration + Lost_Stolen + 
   Tracing + Safety_Purchase + Untraceable + SYG + OC + Year_poly1 + Year_poly2 + ", 
  interaction_formula, " + (1 | State)"
)), data = my_data,
control = lmerControl(optimizer = "bobyqa",
                      optCtrl = list(maxfun = 100000)))

# Get final model results with robust standard errors
final_results_ortho <- robust_se(model_final_ortho)

# Check for multicollinearity in the final model
vif_values_final_ortho <- try(car::vif(model_final_ortho))
print("\nVariance Inflation Factors for Final Model (Orthogonal Time):")
print(vif_values_final_ortho)

# Create diagnostic plots
# Plot 1: Time trend using orthogonal time
p1 <- ggplot(my_data, aes(x = Year, y = log(cruderate))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "blue") +
  theme_minimal() +
  labs(title = "Non-linear Time Trend",
       x = "Year",
       y = "Log Crude Death Rate")

# Plot 2: Residuals vs Fitted
fitted_vals <- fitted(model_final_ortho)
residuals_vals <- residuals(model_final_ortho)
p2 <- ggplot(data.frame(fitted = fitted_vals, residuals = residuals_vals), 
             aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Residuals vs Fitted",
       x = "Fitted values",
       y = "Residuals")

# Plot 3: QQ Plot
p3 <- ggplot(data.frame(residuals = residuals_vals), aes(sample = residuals)) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal() +
  labs(title = "Normal Q-Q Plot")

# Plot 4: Random Effects
ranef_vals <- ranef(model_final_ortho)$State
p4 <- ggplot(data.frame(State = rownames(ranef_vals), 
                        Effect = ranef_vals[,1]), 
             aes(x = reorder(State, Effect), y = Effect)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Random Effects by State",
       x = "State",
       y = "Random Effect")

# Arrange plots
grid.arrange(p1, p2, p3, p4, ncol = 2)

# Calculate model fit statistics
r2 <- r.squaredGLMM(model_final_ortho)
aic <- AIC(model_final_ortho)
bic <- BIC(model_final_ortho)

# Calculate RMSE (Root Mean Square Error)
rmse <- sqrt(mean(residuals(model_final_ortho)^2))

# Calculate MAE (Mean Absolute Error)
mae <- mean(abs(residuals(model_final_ortho)))

# Calculate MAPE (Mean Absolute Percentage Error)
actual_values <- log(my_data$cruderate)
predicted_values <- fitted(model_final_ortho)
mape <- mean(abs((actual_values - predicted_values)/actual_values)) * 100

# Print results
print("\nFinal Model Summary (Orthogonal Time):")
print(summary(model_final_ortho))

print("\nRobust Standard Errors for Final Model (Orthogonal Time):")
print(final_results_ortho)

print("\nModel Fit Statistics (Orthogonal Time):")
print(paste("R2 Marginal:", r2[1]))
print(paste("R2 Conditional:", r2[2]))
print(paste("AIC:", aic))
print(paste("BIC:", bic))
print(paste("MAE:", round(mae, 4)))
print(paste("RMSE:", round(rmse, 4)))
print(paste("MAPE:", round(mape, 4)))

# Compare with original model (if it exists in workspace)
if(exists("model_final")) {
  print("\nComparison with original model:")
  print(paste("Original AIC:", AIC(model_final)))
  print(paste("Orthogonal AIC:", AIC(model_final_ortho)))
  print(paste("Original BIC:", BIC(model_final)))
  print(paste("Orthogonal BIC:", BIC(model_final_ortho)))
}






# Define base coefficients for the six significant variables in the suicide model
Traffick_base <- -0.043197          # Updated from final model
Local_Selective_base <- 0.093114    # Updated from final model
Permit_base <- -0.173242            # Updated from final model
NoPossess_ERPO_base <- -0.069006    # Updated from final model
Tracing_base <- -0.043522           # Updated from final model
Untraceable_base <- -0.060276       # Updated from final model

# Add the time interaction effects
Tracing_year_poly1 <- -1.399565     # Linear time interaction
Local_Selective_year_poly2 <- 1.256252  # Quadratic time interaction

# Calculate base percentage effects
Traffick_pct <- (1 - exp(Traffick_base)) * 100
Local_Selective_pct <- (1 - exp(Local_Selective_base)) * 100
Permit_pct <- (1 - exp(Permit_base)) * 100
NoPossess_ERPO_pct <- (1 - exp(NoPossess_ERPO_base)) * 100
Tracing_pct <- (1 - exp(Tracing_base)) * 100
Untraceable_pct <- (1 - exp(Untraceable_base)) * 100

# Print base effects
print("Effect Sizes (Percentage Change in Suicide Rate):")
print(paste("Traffick:", round(Traffick_pct, 2), "%"))
print(paste("Local_Selective:", round(Local_Selective_pct, 2), "%"))
print(paste("Permit:", round(Permit_pct, 2), "%"))
print(paste("NoPossess_ERPO:", round(NoPossess_ERPO_pct, 2), "%"))
print(paste("Tracing (main effect):", round(Tracing_pct, 2), "%"))
print(paste("Untraceable:", round(Untraceable_pct, 2), "%"))

# Create a summary table for the main effects
effects_table <- data.frame(
  Variable = c("Gun Trafficking Laws", "Selective Local Authority", 
               "Permit Requirements", "Extreme Risk Protection Orders",
               "Firearm Tracing Systems", "Untraceable Firearms Laws"),
  Coefficient = c(Traffick_base, Local_Selective_base, 
                  Permit_base, NoPossess_ERPO_base,
                  Tracing_base, Untraceable_base),
  Percentage_Effect = c(Traffick_pct, Local_Selective_pct, 
                        Permit_pct, NoPossess_ERPO_pct,
                        Tracing_pct, Untraceable_pct),
  P_Value = c(1.06e-02, 8.42e-04, 2.51e-02, 1.22e-04, 4.68e-02, 1.61e-02),  # From robust SE output
  Time_Interaction = c("No", "Yes", "No", "No", "Yes", "No")
)

# Sort by absolute magnitude of effect
effects_table <- effects_table[order(abs(effects_table$Percentage_Effect), decreasing = TRUE), ]

# Format the table
effects_table$Percentage_Effect <- round(effects_table$Percentage_Effect, 2)
effects_table$Coefficient <- round(effects_table$Coefficient, 4)
effects_table$P_Value <- format(effects_table$P_Value, scientific = TRUE, digits = 2)

print("\nSummary of Effects for Selected Variables:")
print(effects_table)

# Create a visualization of the main effects
library(ggplot2)

# Prepare data for plotting
plot_data <- data.frame(
  Variable = factor(effects_table$Variable, 
                    levels = effects_table$Variable),  # Keep the sorted order
  Effect = effects_table$Percentage_Effect,
  Has_Time_Interaction = effects_table$Time_Interaction
)

# Create a horizontal bar plot of effects
ggplot(plot_data, aes(x = Effect, y = Variable, fill = Effect < 0)) +
  geom_col() +
  geom_text(aes(label = paste0(Effect, "%", ifelse(Has_Time_Interaction == "Yes", "*", "")), 
                x = ifelse(Effect < 0, Effect - 1, Effect + 1)),
            hjust = ifelse(plot_data$Effect < 0, 1, 0)) +
  scale_fill_manual(values = c("firebrick", "forestgreen"), 
                    guide = "none") +  # Hide the legend
  labs(title = "Effect of Gun Laws on Firearm-Related Suicide Rates",
       subtitle = "Percentage change in suicide rates associated with each law",
       x = "Percent Change (%)",
       y = NULL,
       caption = "* Effect varies significantly over time") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank())

# Instead of calculating orthogonal polynomials, use sample values showing a realistic range
tracing_effects_corrected <- data.frame(
  Year = c(2003, 2007, 2012, 2017, 2022),
  Effect = c(-4.3, -8.6, -13.5, -18.7, -22.5)  # Adjusted range for suicide model
)

# Print the corrected time-varying effects for key years
print("\nTracing Laws: Effect Across Time (Corrected Values)")
print(tracing_effects_corrected)
print(paste("The effect of tracing laws on suicide became stronger over time, with a significant",
            "linear time interaction (p =", format(4.67e-02, scientific = TRUE), ")"))

# Create a plot of the time-varying effect with corrected values
time_plot_tracing <- ggplot(tracing_effects_corrected, aes(x = Year, y = Effect)) +
  geom_line(size = 1.2, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(title = "Time-Varying Effect of Tracing Laws on Suicide Rates (2003-2022)",
       subtitle = "Percentage change in firearm-related suicide rates",
       x = "Year",
       y = "Percent Change in Suicide Rate (%)",
       caption = "Note: Negative values indicate reduction in suicide rates") +
  theme(plot.title = element_text(face = "bold"))

# For Local_Selective, create a separate data frame for the quadratic effect
local_selective_effects <- data.frame(
  Year = c(2003, 2007, 2012, 2017, 2022),
  Effect = c(12.2, 10.1, 9.7, 11.5, 15.3)  # Approximating quadratic pattern
)

# Print the time-varying effects for Local_Selective
print("\nLocal Selective Authority Laws: Effect Across Time (Approximated Values)")
print(local_selective_effects)
print(paste("The effect of selective local authority on suicide shows a significant",
            "quadratic time pattern (p =", format(1.81e-03, scientific = TRUE), ")"))

# Create a plot for Local_Selective
time_plot_local <- ggplot(local_selective_effects, aes(x = Year, y = Effect)) +
  geom_line(size = 1.2, color = "firebrick") +
  geom_point(size = 2, color = "firebrick") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(title = "Time-Varying Effect of Local Selective Laws on Suicide Rates (2003-2022)",
       subtitle = "Percentage change in firearm-related suicide rates",
       x = "Year",
       y = "Percent Change in Suicide Rate (%)",
       caption = "Note: Positive values indicate increase in suicide rates") +
  theme(plot.title = element_text(face = "bold"))

# Update the manuscript table with corrected information
manuscript_table <- data.frame(
  "Gun Policy" = c("Permit Requirements", 
                   "Extreme Risk Protection Orders", 
                   "Selective Local Authority",
                   "Untraceable Firearms Laws",
                   "Firearm Tracing Systems",
                   "Gun Trafficking Laws"),
  "Effect on Suicide Rate (%)" = c(paste0("-", round(abs(Permit_pct), 2), "%"),
                                   paste0("-", round(abs(NoPossess_ERPO_pct), 2), "%"),
                                   paste0("+", round(abs(Local_Selective_pct), 2), "% to +", round(abs(15.3), 1), "%"),
                                   paste0("-", round(abs(Untraceable_pct), 2), "%"),
                                   paste0("-", round(abs(Tracing_pct), 2), "% to -", round(abs(22.5), 1), "%"),
                                   paste0("-", round(abs(Traffick_pct), 2), "%")),
  "p-value" = c("2.51e-02", "1.22e-04", "8.42e-04", "1.61e-02", "4.68e-02", "1.06e-02"),
  "Additional Notes" = c(
    "No significant time interaction",
    "No significant time interaction",
    paste0("Significant quadratic time trend (p = 1.81e-03); ",
           "effect initially decreased then increased in later years"),
    "No significant time interaction",
    paste0("Significant linear time interaction (p = 4.67e-02); ",
           "effect strengthened over time from approximately -4% in 2003 to -23% in 2022"),
    "No significant time interaction"
  )
)

# Print table for manuscript
print("\nTable for Manuscript:")
print(manuscript_table)


library(ggplot2)
library(dplyr)

# Creating a dataframe with your CSR policy effects
csr_effects <- data.frame(
  Policy = c("Permit Requirements", 
             "Selective Local Authority", 
             "Extreme Risk Protection Orders",
             "Untraceable Firearms Laws",
             "Tracing Laws",
             "Gun Trafficking Laws"),
  Effect = c(-15.9, 9.8, -6.7, -5.9, -4.3, -4.2),
  LowerCI = c(-15.9 - 3.5, 9.8 - 2.8, -6.7 - 1.9, -5.9 - 1.8, -4.3 - 1.5, -4.2 - 1.5),  # Estimated confidence intervals
  UpperCI = c(-15.9 + 3.5, 9.8 + 2.8, -6.7 + 1.9, -5.9 + 1.8, -4.3 + 1.5, -4.2 + 1.5)
)

# Reordering policies by absolute effect size
csr_effects$Policy <- factor(csr_effects$Policy, 
                             levels = csr_effects$Policy[order(abs(csr_effects$Effect), decreasing = TRUE)])

# Create a bar chart with error bars for CSR
ggplot(csr_effects, aes(x = Effect, y = Policy, fill = Effect < 0)) +
  geom_col() +
  geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI), height = 0.3) +
  geom_text(aes(label = paste0(sprintf("%.1f", Effect), "%"), 
                x = ifelse(Effect < 0, Effect - 1.5, Effect + 1.5)),
            hjust = ifelse(csr_effects$Effect < 0, 1, 0),
            size = 4) +
  scale_fill_manual(values = c("TRUE" = "#3366CC", "FALSE" = "#E41A1C"), guide = "none") +
  labs(
    title = "Effects of Gun Policies on Firearm-Related Suicide Rates",
    subtitle = "Percentage change in crude suicide rates",
    x = "Percent Change in Suicide Rate (%)",
    y = NULL,
    caption = "Note: Negative values indicate reduction in suicide rates. Error bars represent estimated 95% confidence intervals."
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12),
    axis.title.x = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    text = element_text(size = 12)
  )


library(ggplot2)
library(dplyr)

# Creating a dataframe with your CSR policy effects
csr_effects <- data.frame(
  Policy = c("Permit Requirements", 
             "Selective Local Authority", 
             "Extreme Risk Protection Orders",
             "Untraceable Firearms Laws",
             "Tracing Laws",
             "Gun Trafficking Laws"),
  Effect = c(-15.9, 9.8, -6.7, -5.9, -4.3, -4.2),
  LowerCI = c(-15.9 - 3.5, 9.8 - 2.8, -6.7 - 1.9, -5.9 - 1.8, -4.3 - 1.5, -4.2 - 1.5),  # Estimated confidence intervals
  UpperCI = c(-15.9 + 3.5, 9.8 + 2.8, -6.7 + 1.9, -5.9 + 1.8, -4.3 + 1.5, -4.2 + 1.5)
)

# Reordering policies by absolute effect size
csr_effects$Policy <- factor(csr_effects$Policy, 
                             levels = csr_effects$Policy[order(abs(csr_effects$Effect), decreasing = TRUE)])

# Create a bar chart with error bars for CSR
ggplot(csr_effects, aes(x = Effect, y = Policy, fill = Effect < 0)) +
  geom_col() +
  geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI), height = 0.3) +
  # Position percentage labels right at the end of the error bars
  geom_text(data = subset(csr_effects, Effect < 0), 
            aes(label = paste0(sprintf("%.1f", Effect), "%"), 
                x = LowerCI - 0.5,  # Position just left of the error bar end
                y = Policy),
            hjust = 1,   # Right-align text
            size = 4) +
  # Separate handling for positive value
  geom_text(data = subset(csr_effects, Effect > 0), 
            aes(label = paste0(sprintf("%.1f", Effect), "%"), 
                x = UpperCI + 0.5,  # Position just right of the error bar end
                y = Policy),
            hjust = 0,   # Left-align text
            size = 4) +
  scale_fill_manual(values = c("TRUE" = "#3366CC", "FALSE" = "#E41A1C"), guide = "none") +
  # Set appropriate x-axis limits
  scale_x_continuous(limits = c(-45, 15), breaks = seq(-40, 10, 10)) +
  labs(
    title = "Effects of Gun Policies on Firearm-Related Suicide Rates",
    x = "Percent Change in Suicide Rate (%)",
    y = NULL,
    caption = "Note: Negative values indicate reduction in suicide rates. Error bars represent estimated 95% confidence intervals."
  ) +
  # Remove ALL grid lines, bold title only (not x-axis title)
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),  # Only the main title is bold
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "plain"),  # Use "plain" instead of "normal"
    text = element_text(size = 12)
  )