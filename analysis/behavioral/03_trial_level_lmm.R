# Trial-level linear mixed-effects analysis of temporal memory error
#
# Required packages:
#   lme4, lmerTest, performance, effectsize, ggplot2
#
# Fits the original trial-level model with signed boundary distance.
# Default input/output paths are relative to this script's directory.
# Optional command-line arguments: input CSV, output directory.

required_packages <- c(
  "lme4",
  "lmerTest",
  "performance",
  "effectsize",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    paste0(
      "Missing required R package(s): ",
      paste(missing_packages, collapse = ", "),
      ". Install them before running this script."
    ),
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_argument)) {
  dirname(normalizePath(sub("^--file=", "", script_argument[1]), mustWork = TRUE))
} else {
  getwd()
}
arguments <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(arguments) >= 1L) arguments[1] else file.path(
  script_dir, "fmri_trials_with_signed_boundary_distance.csv"
)
output_dir <- if (length(arguments) >= 2L) arguments[2] else file.path(
  script_dir, "signed_lmm_results"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file not found: ", input_file,
      ". Supply the input CSV as the first command-line argument."
    ),
    call. = FALSE
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Read and validate the cleaned fMRI trial-level data
# -----------------------------------------------------------------------------

trial_data <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM",
  na.strings = c("", "NA", "NaN")
)

required_columns <- c(
  "subject",
  "trial",
  "signed_error",
  "temporal_location(%)",
  "interval_width",
  "retention_delay",
  "image_similarity",
  "distance_to_boundary",
  "signed_distance_to_boundary",
  "target_frame",
  "nearest_boundary"
)

missing_columns <- setdiff(required_columns, names(trial_data))
if (length(missing_columns) > 0L) {
  stop(
    paste0(
      "The fMRI trial-level file is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

numeric_columns <- setdiff(required_columns, "subject")
for (column_name in numeric_columns) {
  original_values <- trial_data[[column_name]]
  converted_values <- suppressWarnings(as.numeric(original_values))

  conversion_failures <- !is.na(original_values) & is.na(converted_values)
  if (any(conversion_failures)) {
    stop(
      paste0(
        "Column '", column_name,
        "' contains values that cannot be converted to numeric."
      ),
      call. = FALSE
    )
  }

  trial_data[[column_name]] <- converted_values
}

missing_by_column <- vapply(
  trial_data[required_columns],
  function(x) sum(is.na(x)),
  integer(1)
)

if (any(missing_by_column > 0L)) {
  missing_message <- paste0(
    names(missing_by_column)[missing_by_column > 0L],
    "=",
    missing_by_column[missing_by_column > 0L],
    collapse = ", "
  )
  stop(
    paste0(
      "Missing values were found in model variables: ", missing_message,
      ". The script does not silently drop model rows."
    ),
    call. = FALSE
  )
}

expected_locations <- c(20, 40, 60, 80)
observed_locations <- sort(unique(trial_data[["temporal_location(%)"]]))
if (!isTRUE(all.equal(observed_locations, expected_locations))) {
  stop(
    paste0(
      "Unexpected temporal locations. Expected ",
      paste(expected_locations, collapse = ", "),
      "; found ", paste(observed_locations, collapse = ", "), "."
    ),
    call. = FALSE
  )
}

expected_widths <- c(12, 60, 300)
observed_widths <- sort(unique(trial_data$interval_width))
if (!isTRUE(all.equal(observed_widths, expected_widths))) {
  stop(
    paste0(
      "Unexpected interval widths. Expected ",
      paste(expected_widths, collapse = ", "),
      "; found ", paste(observed_widths, collapse = ", "), "."
    ),
    call. = FALSE
  )
}

if (anyDuplicated(trial_data[c("subject", "trial")]) > 0L) {
  stop(
    "Duplicated subject-by-trial rows were found in the fMRI trial data.",
    call. = FALSE
  )
}

if (any(trial_data$signed_distance_to_boundary !=
        trial_data$target_frame - trial_data$nearest_boundary) ||
    any(abs(trial_data$signed_distance_to_boundary) !=
        trial_data$distance_to_boundary)) {
  stop("Signed or absolute boundary-distance validation failed.", call. = FALSE)
}
if (nrow(trial_data) != 4200L || length(unique(trial_data$subject)) != 18L) {
  stop("Expected the original 4200 trials from 18 subjects.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# Prepare the variables exactly as described in the manuscript
# -----------------------------------------------------------------------------

z_score <- function(x, variable_name) {
  variable_sd <- stats::sd(x)
  if (!is.finite(variable_sd) || variable_sd == 0) {
    stop(
      paste0("Cannot standardize '", variable_name, "': SD is zero."),
      call. = FALSE
    )
  }
  as.numeric(scale(x))
}

model_data <- data.frame(
  subject = factor(trial_data$subject),
  error_scaled = trial_data$signed_error / 100,
  temporal_location_ordinal = match(
    trial_data[["temporal_location(%)"]],
    expected_locations
  ),
  interval_width_ordinal = match(
    trial_data$interval_width,
    expected_widths
  ),
  retention_delay_z = z_score(
    trial_data$retention_delay,
    "retention_delay"
  ),
  image_similarity_z = z_score(
    trial_data$image_similarity,
    "image_similarity"
  ),
  signed_distance_to_boundary_z = z_score(
    trial_data$signed_distance_to_boundary,
    "signed_distance_to_boundary"
  )
)

if (anyNA(model_data$temporal_location_ordinal) ||
    anyNA(model_data$interval_width_ordinal)) {
  stop("Ordinal coding produced missing values.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# Fit the trial-level random-intercept model
# -----------------------------------------------------------------------------

lmm_formula <- error_scaled ~
  signed_distance_to_boundary_z +
  retention_delay_z +
  image_similarity_z +
  temporal_location_ordinal +
  interval_width_ordinal +
  (1 | subject)

lmm_model <- lmerTest::lmer(
  formula = lmm_formula,
  data = model_data,
  REML = TRUE,
  na.action = na.fail
)

# Reproduce the original absolute-distance model on the same trials.
baseline_data <- model_data
baseline_data$distance_to_boundary_z <- z_score(
  trial_data$distance_to_boundary, "distance_to_boundary"
)
baseline_formula <- error_scaled ~ distance_to_boundary_z +
  retention_delay_z + image_similarity_z + temporal_location_ordinal +
  interval_width_ordinal + (1 | subject)
baseline_model <- lmerTest::lmer(
  baseline_formula, data = baseline_data, REML = TRUE, na.action = na.fail
)
reference_file <- "G:/dataforgithub_formal/results/behavioral/trial_level_lmm_fixed_effects.csv"
baseline_coefficients <- as.data.frame(summary(baseline_model)$coefficients)
reference_max_difference <- NA_real_
if (file.exists(reference_file)) {
  reference <- read.csv(reference_file, stringsAsFactors = FALSE)
  reproduced <- baseline_coefficients[reference$term, "Estimate"]
  reference_max_difference <- max(abs(reproduced - reference$estimate))
  if (!is.finite(reference_max_difference) || reference_max_difference > 1e-7) {
    stop("The original model's coefficients were not reproduced.", call. = FALSE)
  }
}
convergence_messages <- lmm_model@optinfo$conv$lme4$messages
if (length(convergence_messages) || any(lmm_model@optinfo$conv$opt != 0)) {
  stop("Signed-distance model failed convergence checks.", call. = FALSE)
}
write.csv(model_data, file.path(output_dir, "trial_level_lmm_model_data.csv"), row.names = FALSE)
scaling_variables <- c("signed_distance_to_boundary", "retention_delay", "image_similarity")
scaling_table <- data.frame(
  variable = scaling_variables,
  mean = vapply(trial_data[scaling_variables], mean, numeric(1)),
  sd = vapply(trial_data[scaling_variables], sd, numeric(1))
)
write.csv(scaling_table, file.path(output_dir, "trial_level_lmm_scaling.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Fixed effects, confidence intervals, model fit, and effect sizes
# -----------------------------------------------------------------------------

coefficient_matrix <- as.data.frame(summary(lmm_model)$coefficients)
coefficient_matrix$term <- rownames(coefficient_matrix)
rownames(coefficient_matrix) <- NULL

fixed_effects <- data.frame(
  term = coefficient_matrix$term,
  estimate = coefficient_matrix[["Estimate"]],
  std_error = coefficient_matrix[["Std. Error"]],
  df = coefficient_matrix[["df"]],
  statistic = coefficient_matrix[["t value"]],
  p_value = coefficient_matrix[["Pr(>|t|)"]],
  stringsAsFactors = FALSE
)

critical_t <- stats::qt(0.975, df = fixed_effects$df)
fixed_effects$conf_low <-
  fixed_effects$estimate - critical_t * fixed_effects$std_error
fixed_effects$conf_high <-
  fixed_effects$estimate + critical_t * fixed_effects$std_error

term_labels <- c(
  "(Intercept)" = "Intercept",
  "signed_distance_to_boundary_z" = "Signed Distance to Boundary",
  "distance_to_boundary_z" = "Absolute Distance to Boundary",
  "retention_delay_z" = "Retention Delay",
  "image_similarity_z" = "Image Similarity",
  "temporal_location_ordinal" = "Temporal Location",
  "interval_width_ordinal" = "Interval Width"
)

fixed_effects$predictor <- unname(term_labels[fixed_effects$term])
fixed_effects$predictor[is.na(fixed_effects$predictor)] <-
  fixed_effects$term[is.na(fixed_effects$predictor)]

r2_result <- performance::r2(lmm_model)
r2_table <- data.frame(
  marginal_r2 = unname(r2_result$R2_marginal),
  conditional_r2 = unname(r2_result$R2_conditional),
  n_observations = stats::nobs(lmm_model),
  n_subjects = nlevels(model_data$subject),
  singular_fit = lme4::isSingular(lmm_model, tol = 1e-4)
)

eta_result <- effectsize::eta_squared(
  lmm_model,
  partial = TRUE,
  ci = 0.95
)
eta_table <- as.data.frame(eta_result)
anova_table <- as.data.frame(anova(lmm_model))
anova_table$term <- rownames(anova_table)
rownames(anova_table) <- NULL
write.csv(anova_table, file.path(output_dir, "trial_level_lmm_anova.csv"), row.names = FALSE)

baseline_terms <- fixed_effects$term
baseline_terms[baseline_terms == "signed_distance_to_boundary_z"] <- "distance_to_boundary_z"
baseline_matched <- baseline_coefficients[baseline_terms, , drop = FALSE]
comparison <- data.frame(
  term = fixed_effects$term,
  predictor = fixed_effects$predictor,
  absolute_model_estimate = baseline_matched[["Estimate"]],
  signed_model_estimate = fixed_effects$estimate,
  absolute_model_p = baseline_matched[["Pr(>|t|)"]],
  signed_model_p = fixed_effects$p_value,
  stringsAsFactors = FALSE
)
write.csv(comparison, file.path(output_dir, "absolute_vs_signed_fixed_effects.csv"), row.names = FALSE)

if ("Parameter" %in% names(eta_table)) {
  eta_table$predictor <- unname(term_labels[eta_table$Parameter])
  missing_labels <- is.na(eta_table$predictor)
  eta_table$predictor[missing_labels] <-
    eta_table$Parameter[missing_labels]
}

# -----------------------------------------------------------------------------
# Save numerical results and the fitted model
# -----------------------------------------------------------------------------

write.csv(
  fixed_effects,
  file.path(output_dir, "trial_level_lmm_fixed_effects.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  eta_table,
  file.path(output_dir, "trial_level_lmm_partial_eta_squared.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  r2_table,
  file.path(output_dir, "trial_level_lmm_model_fit.csv"),
  row.names = FALSE,
  na = ""
)

saveRDS(
  lmm_model,
  file.path(output_dir, "trial_level_lmm_model.rds")
)

report_file <- file.path(output_dir, "trial_level_lmm_report.txt")
report_lines <- capture.output({
  cat("Trial-level linear mixed-effects model: signed boundary distance\n")
  cat("======================================\n\n")
  cat("Input:", input_file, "\n")
  cat("Participants:", nlevels(model_data$subject), "\n")
  cat("Observations:", nrow(model_data), "\n")
  cat("Temporal-location coding: 20, 40, 60, 80 -> 1, 2, 3, 4\n")
  cat("Interval-width coding: 12, 60, 300 s -> 1, 2, 3\n")
  cat("Signed error was divided by 100.\n")
  cat("Signed boundary distance = target_frame - nearest_boundary (frames).\n")
  cat("Negative raw distance: target before boundary; positive: target after boundary.\n")
  cat("Boundary distances validated against both frame columns and the stored absolute distance.\n")
  cat("Maximum baseline coefficient difference from archived original results:", reference_max_difference, "\n")
  cat("Model optimizer convergence code:", lmm_model@optinfo$conv$opt, "\n")
  cat(
    "Continuous covariates were z-standardized across all included trials.\n\n"
  )

  cat("Model formula\n")
  cat("-------------\n")
  print(lmm_formula)

  cat("\nModel summary\n")
  cat("-------------\n")
  print(summary(lmm_model))

  cat("\nFixed-effect tests\n")
  cat("------------------\n")
  print(anova(lmm_model))

  cat("\nModel R-squared\n")
  cat("---------------\n")
  print(r2_result)

  cat("\nPartial eta squared\n")
  cat("-------------------\n")
  print(eta_result)

  cat("\nScaling parameters\n")
  print(scaling_table, row.names = FALSE)
  cat("\nOriginal absolute-distance versus signed-distance coefficients\n")
  print(comparison, row.names = FALSE)
  cat("\nThe two models use different boundary definitions. Their coefficient comparison is descriptive.\n")
  cat("Partial eta squared is an effect size, not an additive share of total variance.\n")

  cat("\nSession information\n")
  cat("-------------------\n")
  print(sessionInfo())
})

writeLines(report_lines, con = report_file, useBytes = TRUE)

# -----------------------------------------------------------------------------
# Reproduce the coefficient figure from the original analysis script
# -----------------------------------------------------------------------------

plot_fixed_effects <- fixed_effects[fixed_effects$term != "(Intercept)", ]
plot_fixed_effects$sig_label <- ifelse(
  plot_fixed_effects$p_value < 0.001,
  "***",
  ifelse(
    plot_fixed_effects$p_value < 0.01,
    "**",
    ifelse(plot_fixed_effects$p_value < 0.05, "*", "")
  )
)

coefficient_plot <- ggplot2::ggplot(
  plot_fixed_effects,
  ggplot2::aes(
    x = reorder(predictor, estimate),
    y = estimate,
    fill = estimate
  )
) +
  ggplot2::geom_bar(
    stat = "identity",
    width = 0.7,
    color = "white"
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = estimate - std_error,
      ymax = estimate + std_error
    ),
    width = 0.2,
    color = "black"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sig_label,
      y = estimate + sign(estimate) * (std_error + 0.02)
    ),
    size = 8,
    vjust = 0.7,
    color = "black"
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#377eb8",
    mid = "#f7f7f7",
    high = "#e41a1c",
    midpoint = 0
  ) +
  ggplot2::coord_flip() +
  ggplot2::theme_minimal(base_size = 18) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 22, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 16),
    axis.title.x = ggplot2::element_text(size = 18),
    axis.title.y = ggplot2::element_text(size = 18),
    axis.text.x = ggplot2::element_text(size = 16),
    axis.text.y = ggplot2::element_text(size = 16),
    legend.position = "none"
  ) +
  ggplot2::labs(
    title = "LMM Coefficients: Signed Boundary Distance",
    subtitle = "Continuous covariates standardized; *p<.05, **p<.01, ***p<.001",
    x = "Predictors",
    y = "Estimate"
  )

ggplot2::ggsave(
  filename = file.path(output_dir, "LMM_coefficients.png"),
  plot = coefficient_plot,
  width = 11,
  height = 5,
  dpi = 600,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(output_dir, "LMM_coefficients.pdf"),
  plot = coefficient_plot,
  width = 11,
  height = 5,
  units = "in",
  bg = "white"
)

# -----------------------------------------------------------------------------
# Reproduce the partial eta-squared figure from the original analysis script
# -----------------------------------------------------------------------------

eta_column <- intersect(c("Eta2_partial", "Eta2"), names(eta_table))
if (length(eta_column) != 1L || !"predictor" %in% names(eta_table)) {
  stop(
    "Could not identify the partial eta-squared column for plotting.",
    call. = FALSE
  )
}

eta_plot_data <- eta_table
eta_plot_data$effect_size <- eta_plot_data[[eta_column]]

eta_plot <- ggplot2::ggplot(
  eta_plot_data,
  ggplot2::aes(
    x = reorder(predictor, effect_size),
    y = effect_size,
    fill = effect_size
  )
) +
  ggplot2::geom_bar(stat = "identity", width = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f%%", effect_size * 100)),
    hjust = -0.1,
    size = 6
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_gradient(low = "#91bfdb", high = "#4575b4") +
  ggplot2::theme_minimal(base_size = 18) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 22, face = "bold"),
    axis.title.x = ggplot2::element_text(size = 18),
    axis.title.y = ggplot2::element_text(size = 18),
    axis.text.x = ggplot2::element_text(size = 16),
    axis.text.y = ggplot2::element_text(size = 16),
    legend.position = "none"
  ) +
  ggplot2::labs(
    title = "LMM Effect Sizes (Partial Eta Squared)",
    x = "Predictors",
    y = expression("Partial Eta Squared (" * eta[p]^2 * ")")
  ) +
  ggplot2::ylim(0, max(eta_plot_data$effect_size, na.rm = TRUE) * 1.2)

ggplot2::ggsave(
  filename = file.path(output_dir, "LMM_eta_squared.png"),
  plot = eta_plot,
  width = 11,
  height = 5,
  dpi = 600,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(output_dir, "LMM_eta_squared.pdf"),
  plot = eta_plot,
  width = 11,
  height = 5,
  units = "in",
  bg = "white"
)

cat("Saved trial-level LMM outputs to:", output_dir, "\n")
