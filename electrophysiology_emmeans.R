# =============================================================================
# Electrophysiology Estimated Marginal Means Analysis
# Model-based treatment comparisons at each voltage using emmeans
# =============================================================================

# Load required packages
required_packages <- c("mgcv", "tidyverse", "emmeans", "ggplot2", "patchwork")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# Resolve namespace conflicts
select <- dplyr::select
filter <- dplyr::filter
summarise <- dplyr::summarise
mutate <- dplyr::mutate
arrange <- dplyr::arrange

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

datos <- read.csv("data_nav_veneno.csv")

datos$Canal <- as.factor(datos$Canal)
datos$Tratamiento <- as.factor(datos$Tratamiento)
datos$Concentracion <- as.factor(datos$Concentracion)
datos$Celula <- as.factor(datos$Celula)
datos$Trat <- datos$Tratamiento

cat("=== DATA LOADED ===\n")
cat("Full data:", nrow(datos), "observations\n\n")

# =============================================================================
# CREATE OUTPUT DIRECTORIES
# =============================================================================

# Main results directory
results_dir <- "results_emmeans"
dir.create(results_dir, showWarnings = FALSE)

# Subdirectories for organization
dirs <- list(
  tables = file.path(results_dir, "tables"),
  plots_emmeans = file.path(results_dir, "plots_emmeans"),
  plots_difference = file.path(results_dir, "plots_difference"),
  plots_heatmaps = file.path(results_dir, "plots_heatmaps"),
  diagnostics = file.path(results_dir, "diagnostics")
)

for (d in dirs) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

cat("Output directories created under:", results_dir, "\n\n")

# =============================================================================
# 2. FIT GAM AND EXTRACT EMMEANS FOR A SINGLE CHANNEL/CONCENTRATION/VARIABLE
# =============================================================================

fit_and_emmeans <- function(data, channel, concentration, response_var,
                            voltage_grid = NULL) {

  # Subset data
  subset_data <- data %>%
    filter(Canal == channel, Concentracion == concentration,
           !is.na(!!sym(response_var))) %>%
    droplevels()

  # Check sufficient data
  if (nrow(subset_data) < 20 ||
      length(unique(subset_data$Tratamiento)) < 2) {
    return(list(
      model = NULL, emm = NULL, pairs = NULL,
      channel = channel, concentration = concentration,
      variable = response_var,
      error = "Insufficient data"
    ))
  }

  # Create cell-level random effect ID
  subset_data$LocalCell <- droplevels(
    interaction(subset_data$Celula, subset_data$Tratamiento)
  )

  n_voltages <- length(unique(subset_data$Voltaje))
  k_main <- min(15, n_voltages - 2)

  # Default voltage grid: every voltage in the data
  if (is.null(voltage_grid)) {
    voltage_grid <- sort(unique(subset_data$Voltaje))
  }
  # Only keep grid voltages within the data range
  v_range <- range(subset_data$Voltaje)
  voltage_grid <- voltage_grid[voltage_grid >= v_range[1] &
                                voltage_grid <= v_range[2]]

  k_var <- min(8, floor(n_voltages / 3))

  # --- Fit the model ---
  # For Corriente: use gaulss() to model heteroscedasticity
  # For Var1/Var2: use gaussian()
  gaulss_model <- NULL
  gaussian_model <- NULL
  family_used <- "gaussian"

  if (response_var == "Corriente") {
    # Try gaulss() first
    gaulss_model <- tryCatch({
      gam(
        list(
          as.formula(paste0(
            response_var,
            " ~ Trat + s(Voltaje, by = Trat, k = ", k_main,
            ") + s(LocalCell, bs = 're')"
          )),
          ~ Trat + s(Voltaje, k = k_var)
        ),
        family = gaulss(),
        data = subset_data,
        method = "REML",
        control = gam.control(maxit = 200, trace = FALSE)
      )
    }, error = function(e) NULL)

    if (!is.null(gaulss_model)) family_used <- "gaulss"
  }

  # Always fit gaussian model (needed for emmeans, also used for Var1/Var2)
  gaussian_model <- tryCatch({
    gam(
      as.formula(paste0(
        response_var,
        " ~ Trat + s(Voltaje, by = Trat, k = ", k_main,
        ") + s(LocalCell, bs = 're')"
      )),
      family = gaussian(),
      data = subset_data,
      method = "REML"
    )
  }, error = function(e) NULL)

  if (is.null(gaussian_model) && is.null(gaulss_model)) {
    return(list(
      model = NULL, gaulss_model = NULL,
      emm = NULL, pairs = NULL,
      channel = channel, concentration = concentration,
      variable = response_var,
      family_used = NA,
      error = "Both gaulss and gaussian models failed to fit"
    ))
  }

  # Use gaulss model for diagnostics when available, gaussian for emmeans
  diag_model <- if (!is.null(gaulss_model)) gaulss_model else gaussian_model
  emm_model <- gaussian_model  # emmeans requires standard family

  # --- Extract emmeans ---
  tryCatch({
    emm <- emmeans(
      emm_model,
      ~ Trat | Voltaje,
      at = list(Voltaje = voltage_grid)
    )

    pair_results <- pairs(emm)

    list(
      model = diag_model,          # gaulss when available (for diagnostics)
      emm_model = emm_model,       # gaussian (for emmeans)
      gaulss_model = gaulss_model, # may be NULL for Var1/Var2
      emm = emm,
      pairs = pair_results,
      channel = channel,
      concentration = concentration,
      variable = response_var,
      family_used = family_used,
      data = subset_data,
      voltage_grid = voltage_grid,
      error = NULL
    )
  }, error = function(e) {
    list(
      model = diag_model,
      emm_model = emm_model,
      gaulss_model = gaulss_model,
      emm = NULL, pairs = NULL,
      channel = channel, concentration = concentration,
      variable = response_var,
      family_used = family_used,
      error = paste("emmeans failed:", as.character(e))
    )
  })
}

# =============================================================================
# 3. EXTRACT RESULTS INTO TIDY DATA FRAMES
# =============================================================================

# Convert emmeans summary to data.frame
tidy_emm <- function(result) {
  if (is.null(result$emm)) return(NULL)

  emm_df <- as.data.frame(summary(result$emm))
  emm_df$Channel <- result$channel
  emm_df$Concentration <- result$concentration
  emm_df$Variable <- result$variable
  emm_df
}

# Convert pairwise comparisons to data.frame with significance
tidy_pairs <- function(result) {
  if (is.null(result$pairs)) return(NULL)

  pair_df <- as.data.frame(summary(result$pairs))
  pair_df$Channel <- result$channel
  pair_df$Concentration <- result$concentration
  pair_df$Variable <- result$variable

  # Add significance stars
  pair_df$sig <- ""
  pair_df$sig[pair_df$p.value < 0.05] <- "*"
  pair_df$sig[pair_df$p.value < 0.01] <- "**"
  pair_df$sig[pair_df$p.value < 0.001] <- "***"

  pair_df
}

# =============================================================================
# 4. PLOTTING FUNCTIONS
# =============================================================================

# Plot estimated marginal means for one channel/concentration/variable
plot_emmeans <- function(result) {
  if (is.null(result$emm)) return(NULL)

  emm_df <- as.data.frame(summary(result$emm))
  pair_df <- as.data.frame(summary(result$pairs))
  pair_df$sig <- ""
  pair_df$sig[pair_df$p.value < 0.05] <- "*"

  # Y-position for significance markers
  y_max <- max(emm_df$emmean + emm_df$SE, na.rm = TRUE)
  y_min <- min(emm_df$emmean - emm_df$SE, na.rm = TRUE)
  y_range <- y_max - y_min

  sig_data <- pair_df %>%
    filter(sig != "") %>%
    mutate(y_pos = y_max + 0.05 * y_range)

  p <- ggplot(emm_df, aes(x = Voltaje, y = emmean, color = Trat, fill = Trat)) +
    geom_ribbon(aes(ymin = emmean - SE, ymax = emmean + SE),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("Control" = "#2166AC", "Veneno" = "#B2182B"),
                       name = "Treatment") +
    scale_fill_manual(values = c("Control" = "#2166AC", "Veneno" = "#B2182B"),
                      name = "Treatment") +
    labs(
      title = paste0(result$channel, " - ", result$concentration, " mM"),
      subtitle = paste0(result$variable, " | Sig voltages: ",
                        nrow(sig_data), "/", nrow(pair_df)),
      x = "Voltage (mV)",
      y = paste0("Estimated ", result$variable)
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8, color = "gray30"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4)
    )

  if (nrow(sig_data) > 0) {
    p <- p + geom_point(data = sig_data, aes(x = Voltaje, y = y_pos),
                        shape = 8, size = 2, color = "black", inherit.aes = FALSE)
  }

  p
}

# Plot pairwise differences (Veneno - Control) with CIs
plot_pairwise_diff <- function(result) {
  if (is.null(result$pairs)) return(NULL)

  # Use confint to get confidence intervals
  pair_df <- tryCatch(
    as.data.frame(confint(result$pairs)),
    error = function(e) as.data.frame(summary(result$pairs))
  )

  # Ensure p.value column exists (confint may not include it)
  if (!"p.value" %in% names(pair_df)) {
    summ_df <- as.data.frame(summary(result$pairs))
    pair_df$p.value <- summ_df$p.value
  }

  pair_df$sig <- ifelse(pair_df$p.value < 0.05, "p < 0.05", "n.s.")

  # Note: emmeans pairs gives Control - Veneno by default
  # Flip sign so positive = Veneno effect
  pair_df$estimate_flipped <- -pair_df$estimate

  # Compute CI from SE if lower.CL/upper.CL are not available
  if ("lower.CL" %in% names(pair_df)) {
    pair_df$lower_flipped <- -pair_df$upper.CL
    pair_df$upper_flipped <- -pair_df$lower.CL
  } else {
    pair_df$lower_flipped <- pair_df$estimate_flipped - 1.96 * pair_df$SE
    pair_df$upper_flipped <- pair_df$estimate_flipped + 1.96 * pair_df$SE
  }

  ggplot(pair_df, aes(x = Voltaje, y = estimate_flipped)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = lower_flipped, ymax = upper_flipped),
                alpha = 0.15, fill = "#7570B3") +
    geom_line(color = "#7570B3", linewidth = 1) +
    geom_point(aes(color = sig), size = 2.5) +
    scale_color_manual(values = c("p < 0.05" = "#E7298A", "n.s." = "#7570B3"),
                       name = "") +
    labs(
      title = paste0(result$channel, " - ", result$concentration, " mM"),
      subtitle = paste0(result$variable, " | Difference (Veneno - Control)"),
      x = "Voltage (mV)",
      y = "Estimated difference"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8, color = "gray30"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4)
    )
}

# =============================================================================
# 5. RUN ANALYSIS FOR A DATASET
# =============================================================================

run_emmeans_analysis <- function(data, data_name, output_suffix) {

  cat(sprintf("\n=============================================================\n"))
  cat(sprintf("EMMEANS ANALYSIS: %s\n", data_name))
  cat(sprintf("=============================================================\n\n"))

  channels <- levels(data$Canal)
  concentrations <- levels(data$Concentracion)
  responses <- c("Corriente", "Var1", "Var2")

  # Voltage grid: every 5 mV within data range
  v_range <- range(data$Voltaje)
  voltage_grid <- seq(v_range[1], v_range[2], by = 5)

  # --- Fit all models and collect results ---
  all_results <- list()
  all_emm_df <- data.frame()
  all_pairs_df <- data.frame()

  for (ch in channels) {
    for (conc in concentrations) {
      for (resp in responses) {
        cat(sprintf("  Fitting %s | %s | %s...\n", ch, conc, resp))

        res <- fit_and_emmeans(data, ch, conc, resp, voltage_grid)
        key <- paste(ch, conc, resp, sep = "_")
        all_results[[key]] <- res

        if (is.null(res$model) && is.null(res$emm_model)) {
          cat(sprintf("    SKIPPED: %s\n", res$error))
          next
        }

        if (!is.null(res$family_used)) {
          cat(sprintf("    Family: %s\n", res$family_used))
        }
        if (!is.null(res$error)) {
          cat(sprintf("    Note: %s\n", res$error))
        }

        emm_tidy <- tidy_emm(res)
        if (!is.null(emm_tidy)) all_emm_df <- rbind(all_emm_df, emm_tidy)

        pair_tidy <- tidy_pairs(res)
        if (!is.null(pair_tidy)) all_pairs_df <- rbind(all_pairs_df, pair_tidy)
      }
    }
  }

  # --- Save tables ---
  write.csv(all_emm_df,
            file.path(dirs$tables, paste0("emmeans_estimates_", output_suffix, ".csv")),
            row.names = FALSE)
  write.csv(all_pairs_df,
            file.path(dirs$tables, paste0("emmeans_pairwise_", output_suffix, ".csv")),
            row.names = FALSE)

  # --- GAM DIAGNOSTICS with gam.check ---
  cat("\n  Creating GAM diagnostics (gam.check)...\n")

  diag_pdf <- file.path(dirs$diagnostics, paste0("gam_diagnostics_", output_suffix, ".pdf"))
  diag_txt <- file.path(dirs$diagnostics, paste0("gam_diagnostics_", output_suffix, ".txt"))

  pdf(diag_pdf, width = 10, height = 10)
  sink(diag_txt)  # Capture text output to file

  for (key in names(all_results)) {
    res <- all_results[[key]]
    if (is.null(res$model)) next

    base_title <- paste0(res$channel, " | Conc ", res$concentration,
                         " mM | ", res$variable)

    # --- Diagnostics for gaulss model (Corriente only, when available) ---
    if (!is.null(res$gaulss_model)) {
      model_title <- paste0(base_title, " [gaulss]")
      cat(sprintf("\n========================================\n"))
      cat(sprintf("MODEL: %s\n", model_title))
      cat(sprintf("========================================\n"))

      par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))
      gam.check(res$gaulss_model)
      mtext(model_title, outer = TRUE, cex = 1.1, font = 2, line = 0.5)
      cat("\n")
    }

    # --- Diagnostics for gaussian model (used for emmeans) ---
    emm_mod <- if (!is.null(res$emm_model)) res$emm_model else res$model
    family_label <- if (!is.null(res$gaulss_model)) "gaussian (emmeans)" else "gaussian"
    model_title <- paste0(base_title, " [", family_label, "]")

    cat(sprintf("\n========================================\n"))
    cat(sprintf("MODEL: %s\n", model_title))
    cat(sprintf("========================================\n"))

    par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))
    gam.check(emm_mod)
    mtext(model_title, outer = TRUE, cex = 1.1, font = 2, line = 0.5)
    cat("\n")
  }

  sink()  # Stop capturing text
  dev.off()
  par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))  # Reset par

  cat(sprintf("  Diagnostic plots saved to: %s\n", diag_pdf))
  cat(sprintf("  Diagnostic text saved to:  %s\n", diag_txt))

  # --- Print pairwise comparison summary ---
  cat("\n=== PAIRWISE COMPARISONS (Control vs Veneno at each voltage) ===\n")

  if (nrow(all_pairs_df) > 0) {
    sig_pairs <- all_pairs_df %>% filter(p.value < 0.05)

    cat(sprintf("Total comparisons: %d\n", nrow(all_pairs_df)))
    cat(sprintf("Significant (p < 0.05): %d\n\n", nrow(sig_pairs)))

    if (nrow(sig_pairs) > 0) {
      print(
        sig_pairs %>%
          select(Channel, Concentration, Variable, Voltaje,
                 estimate, SE, p.value, sig) %>%
          arrange(Channel, Concentration, Variable, Voltaje)
      )
    }
  }

  # --- Create emmeans plots per variable ---
  for (resp in responses) {
    cat(sprintf("\n  Creating emmeans plots for %s...\n", resp))

    emm_plots <- list()
    diff_plots <- list()

    for (ch in channels) {
      for (conc in concentrations) {
        key <- paste(ch, conc, resp, sep = "_")
        res <- all_results[[key]]

        if (!is.null(res$model)) {
          emm_plots[[paste(ch, conc, sep = "_")]] <- plot_emmeans(res)
          diff_plots[[paste(ch, conc, sep = "_")]] <- plot_pairwise_diff(res)
        }
      }
    }

    # Estimated means plots
    emm_plot_list <- emm_plots[!sapply(emm_plots, is.null)]
    if (length(emm_plot_list) > 0) {
      combined <- wrap_plots(emm_plot_list, ncol = 2) +
        plot_annotation(
          title = paste0("Estimated Marginal Means: ", resp, " (", data_name, ")"),
          subtitle = "Mean +/- SE from GAM | * = p < 0.05 (Veneno vs Control) | Blue = Control, Red = Veneno",
          theme = theme(plot.title = element_text(face = "bold", size = 14))
        )

      fname_base <- file.path(dirs$plots_emmeans, paste0("emmeans_", tolower(resp), "_", output_suffix))
      ggsave(paste0(fname_base, ".pdf"), combined, width = 14, height = 20)
      ggsave(paste0(fname_base, ".png"), combined, width = 14, height = 20, dpi = 150)
    }

    # Difference plots
    diff_plot_list <- diff_plots[!sapply(diff_plots, is.null)]
    if (length(diff_plot_list) > 0) {
      combined_diff <- wrap_plots(diff_plot_list, ncol = 2) +
        plot_annotation(
          title = paste0("Pairwise Difference (Veneno - Control): ", resp,
                         " (", data_name, ")"),
          subtitle = "Estimate +/- 95% CI | Pink = p < 0.05",
          theme = theme(plot.title = element_text(face = "bold", size = 14))
        )

      fname_base <- file.path(dirs$plots_difference, paste0("difference_", tolower(resp), "_", output_suffix))
      ggsave(paste0(fname_base, ".pdf"), combined_diff, width = 14, height = 20)
      ggsave(paste0(fname_base, ".png"), combined_diff, width = 14, height = 20, dpi = 150)
    }
  }

  # --- Significance summary heatmap ---
  cat("\n  Creating significance heatmap...\n")

  if (nrow(all_pairs_df) > 0) {
    heatmap_data <- all_pairs_df %>%
      mutate(
        log_p = -log10(p.value),
        log_p = ifelse(is.infinite(log_p), 10, log_p),
        log_p = ifelse(is.na(log_p), 0, log_p),
        Channel_Conc = paste0(Channel, " (", Concentration, "mM)")
      )

    for (resp in responses) {
      resp_data <- heatmap_data %>% filter(Variable == resp)
      if (nrow(resp_data) == 0) next

      p_heatmap <- ggplot(resp_data,
                          aes(x = Voltaje, y = Channel_Conc, fill = log_p)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = sig), size = 3, color = "black") +
        scale_fill_gradient2(
          low = "white", mid = "lightyellow", high = "red",
          midpoint = -log10(0.05), name = "-log10(p)"
        ) +
        scale_x_continuous(breaks = seq(-80, 60, by = 10)) +
        labs(
          title = paste0("Significance Heatmap (emmeans): ", resp,
                         " (", data_name, ")"),
          subtitle = "* p<0.05, ** p<0.01, *** p<0.001",
          x = "Voltage (mV)", y = "Channel (Concentration)"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 10),
          panel.grid = element_blank()
        )

      fname_base <- file.path(dirs$plots_heatmaps, paste0("heatmap_", tolower(resp), "_", output_suffix))
      ggsave(paste0(fname_base, ".pdf"), p_heatmap, width = 14, height = 8)
      ggsave(paste0(fname_base, ".png"), p_heatmap, width = 14, height = 8, dpi = 150)
    }
  }

  return(list(emm = all_emm_df, pairs = all_pairs_df, models = all_results))
}

# =============================================================================
# 6. RUN ANALYSIS
# =============================================================================

results <- run_emmeans_analysis(
  datos, "Full Range (-90 to +65 mV)", "full"
)

# =============================================================================
# 7. SUMMARY
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=============================================================================\n\n")

cat("Models fitted:\n")
cat("  Corriente: gaulss() for diagnostics + gaussian() for emmeans\n")
cat("  Var1/Var2: gaussian()\n")
cat("  Structure: response ~ Trat + s(Voltaje, by=Trat) + s(Cell, bs='re')\n")
cat("  Method: REML\n\n")

cat("emmeans extracts estimated marginal means at each voltage, accounting\n")
cat("for the cell random effect. pairs() tests Control vs Veneno at each\n")
cat("voltage, with proper standard errors from the fitted model.\n\n")

cat("Output directory:", results_dir, "\n\n")
cat("Directory structure:\n")
cat("  tables/\n")
cat("    - emmeans_estimates_*.csv          : Estimated means at each voltage\n")
cat("    - emmeans_pairwise_*.csv           : Pairwise comparisons\n")
cat("  plots_emmeans/\n")
cat("    - emmeans_corriente_*              : EMM plot for current\n")
cat("    - emmeans_var1_*                   : EMM plot for activation\n")
cat("    - emmeans_var2_*                   : EMM plot for inactivation\n")
cat("  plots_difference/\n")
cat("    - difference_corriente_*           : Difference plot for current\n")
cat("    - difference_var1_*                : Difference plot for activation\n")
cat("    - difference_var2_*                : Difference plot for inactivation\n")
cat("  plots_heatmaps/\n")
cat("    - heatmap_*                        : Significance heatmaps\n")
cat("  diagnostics/\n")
cat("    - gam_diagnostics_*.pdf            : gam.check plots (QQ, residuals, etc)\n")
cat("    - gam_diagnostics_*.txt            : gam.check text output (k-index, etc)\n\n")

cat("DIAGNOSTICS NOTE:\n")
cat("  For Corriente, diagnostics are generated for BOTH the gaulss() model\n")
cat("  (which models heteroscedasticity) and the gaussian() model (used for\n")
cat("  emmeans). 
