library(tidyverse)

#From Lee et al (parsed by deepseek)
df <- data.frame(
  risk_factor = c("Less education", "Hearing loss", "TBI", "Hypertension",
                  "Excessive alcohol", "Obesity", "Smoking", "Depression",
                  "Social isolation", "Physical inactivity", "Diabetes",
                  "Air pollution"),
  RR = c(1.59, 1.94, 1.84, 1.61, 1.18, 1.60, 1.59, 1.90, 1.57, 1.39, 1.54, 1.09),
  lower_CI = c(1.26, 1.38, 1.54, 1.16, 1.06, 1.34, 1.15, 1.55, 1.32, 1.16, 1.33, 1.07),
  upper_CI = c(2.01, 2.73, 2.20, 2.24, 1.31, 1.92, 2.20, 2.33, 1.85, 1.67, 1.79, 1.11),
  communality = c(72.0, 83.5, 85.4, 53.4, 59.8, 65.9, 63.5, 56.6, 43.3, 30.0, 54.5, 63.8)
)

# Create the data frame
df_prevalence <- data.frame(
  risk_factor = c("Less education", "Hearing loss", "TBI", "Hypertension",
                  "Excessive alcohol", "Obesity", "Smoking", "Depression",
                  "Social isolation", "Physical inactivity", "Diabetes",
                  "Air pollution"),
  total = c(10.7, 10.8, 17.1, 42.2, 3.6, 44.0, 8.5, 7.4, 11.9,
                       62.8, 28.6, 22.8),
  hispanic = c(27.1, 13.1, 10.3, 38.5, 2.0, 48.3, 6.9, 10.7, 24.0, 68.6, 41.0,
               44.4),
  asian = c(6.4, 6.9, 6.0, 38.5, 0.7, 14.6, 4.9, 4.3, 8.0, 56.6,
                         44.1, 55.2),
  black = c(10.6, 6.5, 9.2, 61.0, 2.7, 54.3, 11.7, 6.6, 12.1,
                         73.2, 37.2, 41.3),
  white = c(5.5, 10.6, 20.1, 39.8, 4.2, 43.5, 8.4, 7.2, 10.8,
                         61.3, 25.4, 17.2)
)
dementiarisk <- df_prevalence |>
  left_join(df) |>
  mutate(logrr = log(RR)) |>
  mutate(sdlog = (log(upper_CI) - log(lower_CI)) / (2*qnorm(0.975))) |>
  select(risk_factor, RR, lower_CI, upper_CI, logrr, sdlog, everything()) |>
  select(-communality)

usethis::use_data(dementiarisk, overwrite = TRUE)

# Create a vector of risk factor names
risk_factors <- c("Less education", "Hearing loss", "TBI", "Hypertension",
                  "Excessive alcohol", "Obesity", "Smoking", "Depression",
                  "Social isolation", "Physical inactivity", "Diabetes",
                  "Air pollution")

# Create the correlation matrix
cor_matrix <- matrix(c(
  1.00, 0.11, -0.15, 0.25, -0.06, 0.18, 0.10, 0.31, 0.24, 0.16, 0.22, 0.52,
  0.11, 1.00, 0.09, -0.02, 0.09, 0.01, 0.03, 0.00, -0.09, 0.05, 0.01, -0.08,
  -0.15, 0.09, 1.00, -0.16, 0.09, -0.03, -0.01, 0.04, -0.11, -0.08, -0.02, -0.13,
  0.25, -0.02, -0.16, 1.00, 0.05, 0.35, -0.10, 0.17, 0.11, 0.14, 0.29, 0.21,
  -0.06, 0.09, 0.09, 0.05, 1.00, -0.12, 0.21, 0.03, -0.03, -0.02, 0.00, -0.10,
  0.18, 0.01, -0.03, 0.35, -0.12, 1.00, -0.15, 0.20, 0.13, 0.21, 0.43, 0.12,
  0.10, 0.03, -0.01, -0.10, 0.21, -0.15, 1.00, 0.15, 0.10, 0.09, -0.05, -0.01,
  0.31, 0.00, 0.04, 0.17, 0.03, 0.20, 0.15, 1.00, 0.21, 0.21, 0.20, 0.29,
  0.24, -0.09, -0.11, 0.11, -0.03, 0.13, 0.10, 0.21, 1.00, 0.07, 0.08, 0.23,
  0.16, 0.05, -0.08, 0.14, -0.02, 0.21, 0.09, 0.21, 0.07, 1.00, 0.15, 0.09,
  0.22, 0.01, -0.02, 0.29, 0.00, 0.43, -0.05, 0.20, 0.08, 0.15, 1.00, 0.20,
  0.52, -0.08, -0.13, 0.21, -0.10, 0.12, -0.01, 0.29, 0.23, 0.09, 0.20, 1.00
), nrow = 12, byrow = TRUE)

# Set row and column names
rownames(cor_matrix) <- risk_factors
colnames(cor_matrix) <- risk_factors

cov_df <- as.data.frame(cor_matrix) |>
  mutate(risk_factor = colnames(cor_matrix)) |>
  left_join(
    dementiarisk |>
      select(risk_factor, sdlog)
  ) |>
  rename(risk_factor_1 = risk_factor) |>
  rename(sd1 = sdlog) |>
  pivot_longer(cols =  `Less education`:`Air pollution`,
               names_to = "risk_factor_2", values_to = "correlation") |>
  left_join(
    dementiarisk |>
      select(risk_factor, sdlog),
    by = c("risk_factor_2" = "risk_factor")
  ) |>
  rename(sd2 = sdlog) |>
  mutate(covariance = correlation*sd1*sd2) |>
  select(risk_factor_1, risk_factor_2, covariance) |>
  pivot_wider(id_cols = risk_factor_1, names_from = risk_factor_2, values_from = covariance) |>
  select(-risk_factor_1) |>
  as.matrix()
rownames(cov_df) <- colnames(cov_df)
dementiacov <- cov_df
usethis::use_data(dementiacov, overwrite = TRUE)

