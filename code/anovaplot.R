# ============================================================
# Linear regression and visualization of MemToolbox mu
# ============================================================

library(readxl)
library(dplyr)
library(ggplot2)


# 1. Read data ------------------------------------------------------------

mem_results <- read_excel(
  path = "memresults.xlsx",
  sheet = "Sheet1"
)

# 2. Prepare variables ----------------------------------------------------

# td is kept as a numeric variable for linear regression.
analysis_data <- mem_results %>%
  select(SubjectID, td, mu) %>%
  mutate(
    SubjectID = factor(SubjectID),
    td = as.numeric(td)
  )

# 3. Fit an overall least-squares regression line ------------------------

linear_model <- lm(
  mu ~ td,
  data = analysis_data
)

# Display the regression results
summary(linear_model)

# Display the intercept and slope
coef(linear_model)

# 4. Descriptive statistics for each temporal location -------------------

summary_data <- analysis_data %>%
  group_by(td) %>%
  summarise(
    n = n(),
    mean_mu = mean(mu),
    sd_mu = sd(mu),
    se_mu = sd_mu / sqrt(n),
    .groups = "drop"
  )

print(summary_data)

# 5. Plot -----------------------------------------------------------------

intercept <- coef(linear_model)[1]
slope <- coef(linear_model)[2]

my_colors <- c( "0.2" = "#6DB290",  "0.4" = "#44948F",  "0.6" = "#24768B",  "0.8" = "#215584")

mu_plot <- ggplot(analysis_data, aes(x = td, y = mu, fill = factor(td))) +
  stat_summary(fun = mean, geom = "bar", width = 0.12, color = "white", linewidth = 0.5, alpha = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.04, color = "black", linewidth = 1.2) +
  geom_jitter(width = 0.03, height = 0, alpha = 0.5, color = "grey30", size = 2) +
  annotate("segment", x = 0.2, xend = 0.8, y = intercept + slope * 0.2, yend = intercept + slope * 0.8, color = "black", linewidth = 0.9, lineend = "butt") +
  scale_fill_manual(values = my_colors) +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8), labels = c("20%", "40%", "60%", "80%"), limits = c(0.1, 0.9)) +
  labs(x = "Location", y = "Temporal Memory Bias (μ)") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(size = 13, color = "black", face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none"
  )

print(mu_plot)

ggsave("mu_linear_fit.png", plot = mu_plot, width = 6, height = 5, dpi = 600)
ggsave("mu_linear_fit.pdf", plot = mu_plot, width = 6, height = 5)