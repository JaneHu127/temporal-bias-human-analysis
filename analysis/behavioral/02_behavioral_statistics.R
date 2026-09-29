# Behavioral statistics replacing the original GraphPad Prism workflow
#
# Required packages:
#   readr, dplyr, tidyr, rstatix, afex
#
# Run this script from the repository root. It uses the published
# MemToolbox parameter table and does not modify any source data.
#
# Reaction-time analyses are intentionally omitted from this script.

required_packages <- c("readr", "dplyr", "tidyr", "rstatix", "afex")

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
# Paths and constants
# -----------------------------------------------------------------------------

input_file <- file.path(
  "data", "behavioral", "memtoolbox_parameters.csv"
)
output_dir <- file.path("results", "behavioral")

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file not found: ", input_file,
      ". Run this script from the repository root."
    ),
    call. = FALSE
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

location_levels <- c(20, 40, 60, 80)
group_levels <- c("fMRI", "replication", "control")

# -----------------------------------------------------------------------------
# Read the actual published columns, including the Greek mu and sigma headers
# -----------------------------------------------------------------------------

parameter_data <- readr::read_csv(
  input_file,
  locale = readr::locale(encoding = "UTF-8"),
  name_repair = "minimal",
  show_col_types = FALSE
)
parameter_data <- as.data.frame(parameter_data, stringsAsFactors = FALSE)

expected_columns <- c(
  "subject",
  "group",
  "temporal_location(%)",
  "\u03bc",
  "\u03c3",
  "g"
)

missing_columns <- setdiff(expected_columns, names(parameter_data))
if (length(missing_columns) > 0L) {
  stop(
    paste0(
      "The parameter file is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

parameter_data <- parameter_data[expected_columns]
names(parameter_data) <- c(
  "subject", "group", "temporal_location", "mu", "sigma", "g"
)

numeric_columns <- c("temporal_location", "mu", "sigma", "g")
for (column_name in numeric_columns) {
  original_values <- parameter_data[[column_name]]
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

  parameter_data[[column_name]] <- converted_values
}

# -----------------------------------------------------------------------------
# Data-integrity checks
# -----------------------------------------------------------------------------

missing_by_column <- vapply(
  parameter_data,
  function(x) sum(is.na(x)),
  integer(1)
)

duplicate_rows <- duplicated(
  parameter_data[c("subject", "group", "temporal_location")]
) | duplicated(
  parameter_data[c("subject", "group", "temporal_location")],
  fromLast = TRUE
)

observed_groups <- unique(parameter_data$group)
missing_groups <- setdiff(group_levels, observed_groups)
unexpected_groups <- setdiff(observed_groups, group_levels)

if (length(missing_groups) > 0L || length(unexpected_groups) > 0L) {
  stop(
    paste0(
      "Unexpected group labels. Missing: ",
      ifelse(length(missing_groups) == 0L, "none", paste(missing_groups, collapse = ", ")),
      "; unexpected: ",
      ifelse(
        length(unexpected_groups) == 0L,
        "none",
        paste(unexpected_groups, collapse = ", ")
      ),
      "."
    ),
    call. = FALSE
  )
}

observed_locations <- sort(unique(parameter_data$temporal_location))
if (!isTRUE(all.equal(observed_locations, location_levels))) {
  stop(
    paste0(
      "Unexpected temporal locations. Expected ",
      paste(location_levels, collapse = ", "),
      "; found ", paste(observed_locations, collapse = ", "), "."
    ),
    call. = FALSE
  )
}

if (any(missing_by_column > 0L)) {
  stop(
    paste0(
      "Missing values were found: ",
      paste0(
        names(missing_by_column), "=", missing_by_column,
        collapse = ", "
      ),
      "."
    ),
    call. = FALSE
  )
}

if (any(duplicate_rows)) {
  duplicated_keys <- unique(
    parameter_data[duplicate_rows, c("subject", "group", "temporal_location")]
  )
  stop(
    paste0(
      "Duplicated subject-by-group-by-location rows were found: ",
      paste(
        apply(duplicated_keys, 1, paste, collapse = "/"),
        collapse = ", "
      ),
      "."
    ),
    call. = FALSE
  )
}

observations_per_participant <- stats::aggregate(
  temporal_location ~ group + subject,
  data = parameter_data,
  FUN = length
)
names(observations_per_participant)[3] <- "n_observations"

incomplete_participants <- observations_per_participant[
  observations_per_participant$n_observations != length(location_levels),
]

if (nrow(incomplete_participants) > 0L) {
  stop(
    paste0(
      "Each participant must have exactly four parameter rows. Problematic IDs: ",
      paste(
        paste0(
          incomplete_participants$group,
          "/",
          incomplete_participants$subject,
          " (n=",
          incomplete_participants$n_observations,
          ")"
        ),
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

location_completeness <- stats::aggregate(
  temporal_location ~ group + subject,
  data = parameter_data,
  FUN = function(x) identical(sort(x), location_levels)
)

if (!all(location_completeness$temporal_location)) {
  stop(
    "At least one participant does not have the complete 20/40/60/80 design.",
    call. = FALSE
  )
}

participant_counts <- stats::aggregate(
  subject ~ group,
  data = unique(parameter_data[c("group", "subject")]),
  FUN = length
)
names(participant_counts)[2] <- "n_participants"
participant_counts$group <- factor(participant_counts$group, levels = group_levels)
participant_counts <- participant_counts[order(participant_counts$group), ]
participant_counts$group <- as.character(participant_counts$group)

parameter_data$temporal_location <- factor(
  parameter_data$temporal_location,
  levels = location_levels,
  ordered = TRUE
)

# -----------------------------------------------------------------------------
# Participant-level early/late and overall bias measures
# -----------------------------------------------------------------------------

participant_summary <- parameter_data |>
  tidyr::pivot_wider(
    id_cols = c(subject, group),
    names_from = temporal_location,
    values_from = c(mu, sigma, g),
    names_glue = "{.value}_{temporal_location}"
  ) |>
  dplyr::mutate(
    early_abs_mu = (abs(mu_20) + abs(mu_40)) / 2,
    late_abs_mu = (abs(mu_60) + abs(mu_80)) / 2,
    early_late_abs_mu_difference = early_abs_mu - late_abs_mu,
    early_sigma = (sigma_20 + sigma_40) / 2,
    late_sigma = (sigma_60 + sigma_80) / 2,
    early_g = (g_20 + g_40) / 2,
    late_g = (g_60 + g_80) / 2,
    overall_bias = (mu_20 + mu_40 + mu_60 + mu_80) / 4
  ) |>
  dplyr::mutate(group_order = match(group, group_levels)) |>
  dplyr::arrange(group_order, subject) |>
  dplyr::select(-group_order)

summary_output <- participant_summary |>
  dplyr::select(
    subject,
    group,
    early_abs_mu,
    late_abs_mu,
    early_sigma,
    late_sigma,
    early_g,
    late_g,
    overall_bias,
    early_late_abs_mu_difference
  )

readr::write_csv(
  summary_output,
  file.path(output_dir, "early_late_bias_summary.csv"),
  na = ""
)

# -----------------------------------------------------------------------------
# Output helpers
# -----------------------------------------------------------------------------

results_list <- list()
report_lines <- character()

add_report <- function(...) {
  report_lines <<- c(report_lines, paste0(...))
}

format_number <- function(x, digits = 4L) {
  if (length(x) == 0L || is.na(x)) {
    return("NA")
  }
  formatC(x, format = "f", digits = digits)
}

format_p <- function(p) {
  if (length(p) == 0L || is.na(p)) {
    return("NA")
  }
  if (p < 0.0001) {
    return(formatC(p, format = "e", digits = 6))
  }
  formatC(p, format = "f", digits = 6)
}

add_result <- function(
    analysis,
    group,
    dependent_variable,
    comparison,
    test,
    statistic = NA_real_,
    df1 = NA_real_,
    df2 = NA_real_,
    p_value = NA_real_,
    effect_size = NA_real_,
    effect_size_type = NA_character_,
    mean_1 = NA_real_,
    mean_2 = NA_real_,
    sd_1 = NA_real_,
    sd_2 = NA_real_,
    n = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    sem = NA_real_,
    p_adjustment = "none",
    correction = "none",
    notes = "") {
  results_list[[length(results_list) + 1L]] <<- data.frame(
    analysis = analysis,
    group = group,
    dependent_variable = dependent_variable,
    comparison = comparison,
    test = test,
    statistic = statistic,
    df1 = df1,
    df2 = df2,
    p_value = p_value,
    effect_size = effect_size,
    effect_size_type = effect_size_type,
    mean_1 = mean_1,
    mean_2 = mean_2,
    sd_1 = sd_1,
    sd_2 = sd_2,
    n = n,
    ci_low = ci_low,
    ci_high = ci_high,
    sem = sem,
    p_adjustment = p_adjustment,
    correction = correction,
    notes = notes,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Repeated-measures ANOVA and pairwise comparisons for mu
# -----------------------------------------------------------------------------

run_mu_rm_anova <- function(group_name) {
  group_data <- parameter_data[parameter_data$group == group_name, ]
  group_data$subject <- factor(group_data$subject)
  group_data$temporal_location <- factor(
    group_data$temporal_location,
    levels = location_levels,
    ordered = FALSE
  )

  anova_result <- rstatix::anova_test(
    data = group_data,
    dv = mu,
    wid = subject,
    within = temporal_location,
    effect.size = "pes"
  )

  mauchly <- anova_result[["Mauchly's Test for Sphericity"]]
  sphericity <- anova_result[["Sphericity Corrections"]]
  violated <- nrow(mauchly) > 0L && mauchly$p[1] < 0.05
  selected_anova <- rstatix::get_anova_table(
    anova_result,
    correction = "auto"
  )

  f_value <- selected_anova$F[1]
  df1 <- selected_anova$DFn[1]
  df2 <- selected_anova$DFd[1]
  p_value <- stats::pf(f_value, df1, df2, lower.tail = FALSE)
  partial_eta_squared <- selected_anova$pes[1]
  correction_label <- if (violated) "Greenhouse-Geisser" else "none"
  gg_epsilon <- if (nrow(sphericity) > 0L) sphericity$GGe[1] else NA_real_

  add_result(
    analysis = paste0(group_name, " temporal-location effect on mu"),
    group = group_name,
    dependent_variable = "mu",
    comparison = "20 vs 40 vs 60 vs 80",
    test = "one-factor repeated-measures ANOVA",
    statistic = f_value,
    df1 = df1,
    df2 = df2,
    p_value = p_value,
    effect_size = partial_eta_squared,
    effect_size_type = "partial eta squared",
    n = length(unique(group_data$subject)),
    correction = correction_label,
    notes = ifelse(
      is.na(gg_epsilon),
      "",
      paste0("Greenhouse-Geisser epsilon=", format_number(gg_epsilon, 6L))
    )
  )

  if (nrow(mauchly) > 0L) {
    add_result(
      analysis = paste0(group_name, " temporal-location effect on mu"),
      group = group_name,
      dependent_variable = "mu",
      comparison = "sphericity",
      test = "Mauchly's test",
      statistic = mauchly$W[1],
      p_value = mauchly$p[1],
      n = length(unique(group_data$subject)),
      notes = ifelse(violated, "sphericity violated", "sphericity not violated")
    )
  }

  add_report("")
  add_report(group_name, ": mu across temporal locations")
  add_report(strrep("-", nchar(group_name) + 30L))
  add_report(
    "Mauchly W=", format_number(mauchly$W[1]),
    ", p=", format_p(mauchly$p[1]), "."
  )
  add_report(
    if (violated) {
      paste0(
        "Sphericity was violated; Greenhouse-Geisser correction was used ",
        "(epsilon=", format_number(gg_epsilon, 6L), ")."
      )
    } else {
      "Sphericity was not violated; unadjusted degrees of freedom were used."
    }
  )
  add_report(
    "F(", format_number(df1, 3L), ", ", format_number(df2, 3L), ")=",
    format_number(f_value), ", p=", format_p(p_value),
    ", partial eta squared=", format_number(partial_eta_squared), "."
  )

  pair_indices <- utils::combn(location_levels, 2L)
  pairwise_rows <- vector("list", ncol(pair_indices))

  for (pair_index in seq_len(ncol(pair_indices))) {
    location_1 <- pair_indices[1, pair_index]
    location_2 <- pair_indices[2, pair_index]

    values_1 <- group_data[
      as.numeric(as.character(group_data$temporal_location)) == location_1,
      c("subject", "mu")
    ]
    values_2 <- group_data[
      as.numeric(as.character(group_data$temporal_location)) == location_2,
      c("subject", "mu")
    ]
    names(values_1)[2] <- "value_1"
    names(values_2)[2] <- "value_2"
    paired_values <- merge(values_1, values_2, by = "subject", sort = TRUE)

    test_result <- stats::t.test(
      paired_values$value_1,
      paired_values$value_2,
      paired = TRUE,
      alternative = "two.sided"
    )
    difference_values <- paired_values$value_1 - paired_values$value_2
    dz <- mean(difference_values) / stats::sd(difference_values)

    pairwise_rows[[pair_index]] <- data.frame(
      location_1 = location_1,
      location_2 = location_2,
      statistic = unname(test_result$statistic),
      df = unname(test_result$parameter),
      p_raw = test_result$p.value,
      effect_size = dz,
      mean_1 = mean(paired_values$value_1),
      mean_2 = mean(paired_values$value_2),
      sd_1 = stats::sd(paired_values$value_1),
      sd_2 = stats::sd(paired_values$value_2),
      n = nrow(paired_values)
    )
  }

  pairwise_table <- do.call(rbind, pairwise_rows)
  pairwise_table$p_adjusted <- stats::p.adjust(
    pairwise_table$p_raw,
    method = "holm"
  )

  add_report("Pairwise paired t-tests (Holm-adjusted p values):")

  for (row_index in seq_len(nrow(pairwise_table))) {
    pair_row <- pairwise_table[row_index, ]

    add_result(
      analysis = paste0(group_name, " mu pairwise comparison"),
      group = group_name,
      dependent_variable = "mu",
      comparison = paste0(
        pair_row$location_1, "% vs ", pair_row$location_2, "%"
      ),
      test = "paired-samples t-test",
      statistic = pair_row$statistic,
      df1 = pair_row$df,
      p_value = pair_row$p_adjusted,
      effect_size = pair_row$effect_size,
      effect_size_type = "Cohen's dz",
      mean_1 = pair_row$mean_1,
      mean_2 = pair_row$mean_2,
      sd_1 = pair_row$sd_1,
      sd_2 = pair_row$sd_2,
      n = pair_row$n,
      p_adjustment = "Holm",
      notes = paste0("unadjusted p=", format_p(pair_row$p_raw))
    )

    add_report(
      "  ", pair_row$location_1, "% vs ", pair_row$location_2,
      "%: t(", format_number(pair_row$df, 0L), ")=",
      format_number(pair_row$statistic),
      ", raw p=", format_p(pair_row$p_raw),
      ", Holm p=", format_p(pair_row$p_adjusted),
      ", dz=", format_number(pair_row$effect_size), "."
    )
  }

  list(
    statistic = f_value,
    df1 = df1,
    df2 = df2,
    p_value = p_value,
    effect_size = partial_eta_squared
  )
}

# -----------------------------------------------------------------------------
# Participant-level paired t-tests
# -----------------------------------------------------------------------------

run_paired_test <- function(
    group_name,
    early_variable,
    late_variable,
    dependent_variable,
    analysis_label) {
  group_summary <- participant_summary[
    participant_summary$group == group_name,
  ]

  early_values <- group_summary[[early_variable]]
  late_values <- group_summary[[late_variable]]
  difference_values <- early_values - late_values

  test_result <- stats::t.test(
    early_values,
    late_values,
    paired = TRUE,
    alternative = "two.sided"
  )
  dz <- mean(difference_values) / stats::sd(difference_values)

  add_result(
    analysis = analysis_label,
    group = group_name,
    dependent_variable = dependent_variable,
    comparison = "early vs late",
    test = "paired-samples t-test",
    statistic = unname(test_result$statistic),
    df1 = unname(test_result$parameter),
    p_value = test_result$p.value,
    effect_size = dz,
    effect_size_type = "Cohen's dz",
    mean_1 = mean(early_values),
    mean_2 = mean(late_values),
    sd_1 = stats::sd(early_values),
    sd_2 = stats::sd(late_values),
    n = length(early_values),
    ci_low = test_result$conf.int[1],
    ci_high = test_result$conf.int[2]
  )

  add_report("")
  add_report(analysis_label)
  add_report(strrep("-", nchar(analysis_label)))
  add_report(
    "Early: M=", format_number(mean(early_values)),
    ", SD=", format_number(stats::sd(early_values)),
    "; late: M=", format_number(mean(late_values)),
    ", SD=", format_number(stats::sd(late_values)), "."
  )
  add_report(
    "t(", format_number(unname(test_result$parameter), 0L), ")=",
    format_number(unname(test_result$statistic)),
    ", p=", format_p(test_result$p.value),
    ", dz=", format_number(dz),
    ", 95% CI of early-late difference [",
    format_number(test_result$conf.int[1]), ", ",
    format_number(test_result$conf.int[2]), "]."
  )

  unname(test_result$statistic)
}

# -----------------------------------------------------------------------------
# Control one-sample test of the manuscript-compatible overall signed bias
# -----------------------------------------------------------------------------

run_control_one_sample_test <- function() {
  control_summary <- participant_summary[
    participant_summary$group == "control",
  ]
  bias_values <- control_summary$overall_bias

  test_result <- stats::t.test(
    bias_values,
    mu = 0,
    alternative = "two.sided"
  )
  cohen_d <- mean(bias_values) / stats::sd(bias_values)
  sem_value <- stats::sd(bias_values) / sqrt(length(bias_values))

  add_result(
    analysis = "control overall bias vs zero",
    group = "control",
    dependent_variable = "overall signed mu",
    comparison = "mean(mu20, mu40, mu60, mu80) vs zero",
    test = "one-sample t-test",
    statistic = unname(test_result$statistic),
    df1 = unname(test_result$parameter),
    p_value = test_result$p.value,
    effect_size = cohen_d,
    effect_size_type = "Cohen's d",
    mean_1 = mean(bias_values),
    sd_1 = stats::sd(bias_values),
    n = length(bias_values),
    ci_low = test_result$conf.int[1],
    ci_high = test_result$conf.int[2],
    sem = sem_value,
    notes = "overall bias is the participant mean of signed mu at all four locations"
  )

  add_report("")
  add_report("Control group: overall signed bias vs zero")
  add_report("------------------------------------------")
  add_report(
    "Overall bias was defined as mean(mu20, mu40, mu60, mu80)."
  )
  add_report(
    "M=", format_number(mean(bias_values)),
    ", SD=", format_number(stats::sd(bias_values)),
    ", SEM=", format_number(sem_value),
    ", t(", format_number(unname(test_result$parameter), 0L), ")=",
    format_number(unname(test_result$statistic)),
    ", p=", format_p(test_result$p.value),
    ", 95% CI [", format_number(test_result$conf.int[1]),
    ", ", format_number(test_result$conf.int[2]), "]."
  )
  add_report(
    "Early absolute mu: M=",
    format_number(mean(control_summary$early_abs_mu)),
    ", SD=", format_number(stats::sd(control_summary$early_abs_mu)),
    "; late absolute mu: M=",
    format_number(mean(control_summary$late_abs_mu)),
    ", SD=", format_number(stats::sd(control_summary$late_abs_mu)), "."
  )

  unname(test_result$statistic)
}

# -----------------------------------------------------------------------------
# Two Group x Location mixed ANOVAs on early/late absolute mu
# -----------------------------------------------------------------------------

run_mixed_anova <- function(memory_group) {
  comparison_groups <- c("control", memory_group)
  comparison_summary <- participant_summary[
    participant_summary$group %in% comparison_groups,
  ]

  mixed_data <- rbind(
    data.frame(
      participant_id = paste(
        comparison_summary$group,
        comparison_summary$subject,
        sep = "::"
      ),
      group = comparison_summary$group,
      location = "early",
      absolute_mu = comparison_summary$early_abs_mu
    ),
    data.frame(
      participant_id = paste(
        comparison_summary$group,
        comparison_summary$subject,
        sep = "::"
      ),
      group = comparison_summary$group,
      location = "late",
      absolute_mu = comparison_summary$late_abs_mu
    )
  )

  mixed_data$group <- factor(
    mixed_data$group,
    levels = c("control", memory_group)
  )
  mixed_data$location <- factor(
    mixed_data$location,
    levels = c("early", "late")
  )

  mixed_model <- suppressMessages(
    afex::aov_ez(
      id = "participant_id",
      dv = "absolute_mu",
      data = mixed_data,
      between = "group",
      within = "location",
      type = 3,
      anova_table = list(correction = "none", es = "pes")
    )
  )

  anova_table <- as.data.frame(mixed_model$anova_table)
  anova_table$effect <- rownames(anova_table)
  rownames(anova_table) <- NULL

  effect_labels <- c(
    "group" = "Group main effect",
    "location" = "Location main effect",
    "group:location" = "Group x Location interaction"
  )

  memory_values <- comparison_summary[
    comparison_summary$group == memory_group,
  ]
  control_values <- comparison_summary[
    comparison_summary$group == "control",
  ]

  cell_note <- paste0(
    memory_group, " early=", format_number(mean(memory_values$early_abs_mu)),
    ", ", memory_group, " late=", format_number(mean(memory_values$late_abs_mu)),
    "; control early=", format_number(mean(control_values$early_abs_mu)),
    ", control late=", format_number(mean(control_values$late_abs_mu))
  )

  add_report("")
  add_report(memory_group, " vs control: Group x Location mixed ANOVA")
  add_report(strrep("-", nchar(memory_group) + 42L))
  add_report(cell_note, ".")

  for (row_index in seq_len(nrow(anova_table))) {
    effect_name <- anova_table$effect[row_index]
    p_value <- anova_table[["Pr(>F)"]][row_index]

    add_result(
      analysis = paste0(memory_group, " vs control mixed ANOVA"),
      group = paste0(memory_group, " vs control"),
      dependent_variable = "early/late absolute mu",
      comparison = unname(effect_labels[effect_name]),
      test = "mixed-design ANOVA",
      statistic = anova_table$F[row_index],
      df1 = anova_table[["num Df"]][row_index],
      df2 = anova_table[["den Df"]][row_index],
      p_value = p_value,
      effect_size = anova_table$pes[row_index],
      effect_size_type = "partial eta squared",
      n = nrow(comparison_summary),
      notes = cell_note
    )

    add_report(
      unname(effect_labels[effect_name]), ": F(",
      format_number(anova_table[["num Df"]][row_index], 0L), ", ",
      format_number(anova_table[["den Df"]][row_index], 0L), ")=",
      format_number(anova_table$F[row_index]),
      ", p=", format_p(p_value),
      ", partial eta squared=",
      format_number(anova_table$pes[row_index]), "."
    )
  }

  interaction_row <- anova_table[anova_table$effect == "group:location", ]
  interaction_row$F[1]
}

# -----------------------------------------------------------------------------
# Between-group comparisons of the overall signed bias measure
# -----------------------------------------------------------------------------

run_between_group_test <- function(memory_group) {
  memory_values <- participant_summary$overall_bias[
    participant_summary$group == memory_group
  ]
  control_values <- participant_summary$overall_bias[
    participant_summary$group == "control"
  ]

  test_result <- stats::t.test(
    memory_values,
    control_values,
    paired = FALSE,
    alternative = "two.sided",
    var.equal = TRUE
  )

  n_memory <- length(memory_values)
  n_control <- length(control_values)
  pooled_sd <- sqrt(
    ((n_memory - 1) * stats::var(memory_values) +
      (n_control - 1) * stats::var(control_values)) /
      (n_memory + n_control - 2)
  )
  cohen_d <- (mean(memory_values) - mean(control_values)) / pooled_sd

  add_result(
    analysis = paste0(memory_group, " vs control overall bias"),
    group = paste0(memory_group, " vs control"),
    dependent_variable = "overall signed mu",
    comparison = paste0(memory_group, " vs control"),
    test = "independent-samples Student t-test",
    statistic = unname(test_result$statistic),
    df1 = unname(test_result$parameter),
    p_value = test_result$p.value,
    effect_size = cohen_d,
    effect_size_type = "Cohen's d (pooled SD)",
    mean_1 = mean(memory_values),
    mean_2 = mean(control_values),
    sd_1 = stats::sd(memory_values),
    sd_2 = stats::sd(control_values),
    n = n_memory + n_control,
    ci_low = test_result$conf.int[1],
    ci_high = test_result$conf.int[2],
    notes = paste0(
      "n_", memory_group, "=", n_memory,
      "; n_control=", n_control,
      "; overall bias=mean signed mu across four locations"
    )
  )

  add_report("")
  add_report(memory_group, " vs control: overall signed bias")
  add_report(strrep("-", nchar(memory_group) + 35L))
  add_report(
    memory_group, ": M=", format_number(mean(memory_values)),
    ", SD=", format_number(stats::sd(memory_values)),
    "; control: M=", format_number(mean(control_values)),
    ", SD=", format_number(stats::sd(control_values)), "."
  )
  add_report(
    "Student t(", format_number(unname(test_result$parameter), 0L), ")=",
    format_number(unname(test_result$statistic)),
    ", p=", format_p(test_result$p.value),
    ", Cohen's d=", format_number(cohen_d),
    ", 95% CI of mean difference [",
    format_number(test_result$conf.int[1]), ", ",
    format_number(test_result$conf.int[2]), "]."
  )

  unname(test_result$statistic)
}

# -----------------------------------------------------------------------------
# Run all requested analyses except reaction time
# -----------------------------------------------------------------------------

add_report("Behavioral statistical analyses")
add_report("===============================")
add_report("")
add_report("DATA CHECKS")
add_report("-----------")
add_report("Input: ", input_file)
add_report("Participants by group:")
for (row_index in seq_len(nrow(participant_counts))) {
  add_report(
    "  ", participant_counts$group[row_index], ": ",
    participant_counts$n_participants[row_index]
  )
}
add_report(
  "Observations per participant: min=",
  min(observations_per_participant$n_observations),
  ", max=", max(observations_per_participant$n_observations), "."
)
add_report(
  "Available temporal locations: ",
  paste(location_levels, collapse = ", "), "."
)
add_report(
  "Missing values: ",
  paste0(names(missing_by_column), "=", missing_by_column, collapse = ", "),
  "."
)
add_report("Duplicated subject x location rows: ", sum(duplicate_rows), ".")
add_report("Reaction-time analysis: omitted by instruction.")

fmri_anova <- run_mu_rm_anova("fMRI")
fmri_abs_mu_t <- run_paired_test(
  "fMRI",
  "early_abs_mu",
  "late_abs_mu",
  "absolute mu",
  "fMRI early vs late absolute mu"
)
fmri_sigma_t <- run_paired_test(
  "fMRI",
  "early_sigma",
  "late_sigma",
  "sigma",
  "fMRI early vs late sigma"
)
fmri_g_t <- run_paired_test(
  "fMRI",
  "early_g",
  "late_g",
  "g",
  "fMRI early vs late guess rate"
)

replication_anova <- run_mu_rm_anova("replication")
replication_abs_mu_t <- run_paired_test(
  "replication",
  "early_abs_mu",
  "late_abs_mu",
  "absolute mu",
  "replication early vs late absolute mu"
)

control_t <- run_control_one_sample_test()
fmri_control_interaction <- run_mixed_anova("fMRI")
replication_control_interaction <- run_mixed_anova("replication")
fmri_control_t <- run_between_group_test("fMRI")
replication_control_t <- run_between_group_test("replication")

# -----------------------------------------------------------------------------
# Validation against approximate manuscript values (never used in calculations)
# -----------------------------------------------------------------------------

validation_table <- data.frame(
  analysis = c(
    "fMRI mu ANOVA F",
    "fMRI early-late absolute mu t",
    "fMRI sigma t",
    "fMRI g t",
    "replication mu ANOVA F",
    "replication early-late absolute mu t",
    "control bias vs zero t",
    "fMRI-control interaction F",
    "replication-control interaction F",
    "fMRI-control overall bias t",
    "replication-control overall bias t"
  ),
  manuscript_target = c(
    48.06,
    5.06,
    3.25,
    -1.35,
    182.15,
    8.33,
    0.92,
    4.643,
    7.178,
    2.78,
    3.21
  ),
  observed = c(
    fmri_anova$statistic,
    fmri_abs_mu_t,
    fmri_sigma_t,
    fmri_g_t,
    replication_anova$statistic,
    replication_abs_mu_t,
    control_t,
    fmri_control_interaction,
    replication_control_interaction,
    fmri_control_t,
    replication_control_t
  ),
  stringsAsFactors = FALSE
)

validation_tolerance <- pmax(
  0.10,
  abs(validation_table$manuscript_target) * 0.10
)
validation_table$status <- ifelse(
  abs(validation_table$observed - validation_table$manuscript_target) <=
    validation_tolerance,
  "approximately reproduced",
  "differs from manuscript target"
)

add_report("")
add_report("VALIDATION AGAINST MANUSCRIPT")
add_report("-----------------------------")
add_report(
  "Approximate agreement is defined here only as an absolute difference ",
  "within max(0.10, 10% of the manuscript statistic)."
)

for (row_index in seq_len(nrow(validation_table))) {
  add_report(
    validation_table$analysis[row_index],
    ": observed=", format_number(validation_table$observed[row_index]),
    ", target=", format_number(validation_table$manuscript_target[row_index]),
    " -> ", validation_table$status[row_index], "."
  )
}

if (any(validation_table$status == "differs from manuscript target")) {
  add_report("")
  add_report("Analyses requiring attention:")
  for (analysis_name in validation_table$analysis[
    validation_table$status == "differs from manuscript target"
  ]) {
    add_report("  - ", analysis_name)
  }
  add_report(
    "Possible reasons include a different parameter-data version, participant ",
    "exclusions, Prism/R implementation differences, or rounding. No values ",
    "were altered to force agreement."
  )
} else {
  add_report("All validation statistics were approximately reproduced.")
}

add_report("")
add_report("DECISIONS")
add_report("---------")
add_report(
  "Early/late paired comparisons and mixed ANOVAs use absolute mu: ",
  "early=mean(abs(mu20), abs(mu40)); ",
  "late=mean(abs(mu60), abs(mu80))."
)
add_report(
  "Control-vs-zero and between-group tests use overall signed bias: ",
  "mean(mu20, mu40, mu60, mu80)."
)
add_report("All t-tests are two-tailed.")
add_report("Independent-samples tests use the pooled-variance Student t-test.")
add_report("Pairwise location tests use Holm correction across six comparisons.")
add_report("Reaction-time testing is not included in this script.")

# -----------------------------------------------------------------------------
# Save the statistical table and readable report
# -----------------------------------------------------------------------------

statistics_table <- do.call(rbind, results_list)

readr::write_csv(
  statistics_table,
  file.path(output_dir, "behavioral_statistics.csv"),
  na = ""
)

report_lines <- c(
  report_lines,
  "",
  "SESSION INFORMATION",
  "-------------------",
  capture.output(sessionInfo())
)

report_file <- file.path(output_dir, "behavioral_statistics.txt")
writeLines(report_lines, con = report_file, useBytes = TRUE)

cat(paste(report_lines, collapse = "\n"), "\n")
cat("\nSaved behavioral statistics outputs to:", output_dir, "\n")
