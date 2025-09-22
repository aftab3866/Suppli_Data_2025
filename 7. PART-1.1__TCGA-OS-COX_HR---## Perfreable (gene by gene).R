# ==========================
# Publication-ready Cox Model - ENHANCED with Clinical Forest Plots
# ==========================

# 0. Load libraries with pacman for efficient loading
if (!require(pacman)) install.packages("pacman")
pacman::p_load(survival, survminer, dplyr, forestmodel, rms, ggplot2, stringr, gridExtra, broom, purrr, ggpubr)

# Set seed for reproducibility
set.seed(42)

# Create results directory
if (!dir.exists("results")) {
  dir.create("results")
  message("Created 'results' directory")
}

#====================================
# 0. Data Preparation for Analysis
#====================================

# Load data with error handling
tryCatch({
  tcga_data <- readRDS("tcga_data.RDS")
  limma_res <- readRDS("limma_res.RDS")
}, error = function(e) {
  stop("Error loading RDS files: ", e$message)
})

# Define gene list
my_genes <- c("RRM2","MYBL2","TPX2","DES","KIF18B","SYNM","SLC7A11",
              "TNS4","FOLR1","ABCG2","SORBS1","MASP1","CAV1","MYL9",
              "FHL1","RXRG","SKA3","CRYAB","HJURP","OTX1","IQGAP3","RHPN1",
              "TCEAL2","SIM2","PYGM","NRG1")

# Clean clinical data
clinical_data <- as.data.frame(colData(tcga_data))

clinical_clean <- clinical_data %>%
  mutate(
    os_time = coalesce(days_to_death, days_to_last_follow_up),
    os_status = as.numeric(vital_status == "Dead")
  ) %>%
  filter(!is.na(os_time) & !is.na(os_status) & os_time > 0) %>%
  dplyr::select(
    patient = bcr_patient_barcode,
    os_time,
    os_status,
    age = age_at_diagnosis,
    gender,
    tumor_stage = ajcc_pathologic_stage
  ) %>%
  mutate(
    tumor_stage_simple = case_when(
      str_detect(tumor_stage, "^Stage I$|Stage IA|Stage IB") ~ "I",
      str_detect(tumor_stage, "^Stage II$|Stage IIA|Stage IIB") ~ "II",
      str_detect(tumor_stage, "^Stage III$|Stage IIIA|Stage IIIB|Stage IIIC") ~ "III",
      str_detect(tumor_stage, "^Stage IV$") ~ "IV",
      TRUE ~ NA_character_
    ),
    age = as.numeric(age),
    gender = factor(gender),
    tumor_stage_simple = factor(tumor_stage_simple, levels = c("I", "II", "III", "IV"))
  ) %>%
  filter(complete.cases(age, gender, tumor_stage_simple))

# Extract expression data
expr_matrix <- log2(assay(tcga_data, "tpm_unstrand") + 1)
message("Using TPM unstranded assay with log2 transformation")
message(paste("Expression matrix dimensions:", dim(expr_matrix)[1], "genes x", dim(expr_matrix)[2], "samples"))

# Map gene symbols to ENSEMBL IDs
gene_annotations <- as.data.frame(rowData(tcga_data))
ensembl_ids <- gene_annotations %>%
  filter(gene_name %in% my_genes) %>%
  dplyr::select(gene_id, gene_name)

missing_genes <- setdiff(my_genes, ensembl_ids$gene_name)
if (length(missing_genes) > 0) {
  warning("Missing genes: ", paste(missing_genes, collapse = ", "))
}

# Extract and transpose expression data
gene_expression_data <- expr_matrix[ensembl_ids$gene_id, , drop = FALSE]
gene_expression_data_t <- t(gene_expression_data)
colnames(gene_expression_data_t) <- ensembl_ids$gene_name[match(colnames(gene_expression_data_t), ensembl_ids$gene_id)]

# Integrate data
common_patients <- intersect(rownames(gene_expression_data_t), rownames(clinical_clean))
final_df <- data.frame(
  clinical_clean[common_patients, c("os_time", "os_status", "age", "gender", "tumor_stage_simple")],
  gene_expression_data_t[common_patients, , drop = FALSE]
)

colnames(final_df)[5] <- "tumor_stage"

# Data validation and cleaning
final_df <- final_df %>%
  mutate(
    os_status = as.numeric(os_status),
    tumor_stage = factor(tumor_stage)
  ) %>%
  filter(os_time > 0 & os_time <= quantile(os_time, 0.99, na.rm = TRUE))

# ==========================
# 1. Clinical Variables Only Analysis - FIXED to include tumor_stage
# ==========================
message("\n=== Analyzing Clinical Variables Only ===")

# Create a separate model JUST for forest plot that includes tumor_stage as covariate
# This is for visualization purposes only
clinical_covariates <- c("tumor_stage", "age", "gender")
clinical_formula <- as.formula(paste("Surv(os_time, os_status) ~", paste(clinical_covariates, collapse = " + ")))
clinical_model_for_plot <- coxph(clinical_formula, data = final_df, x = TRUE, y = TRUE)

# Create enhanced forest plot for clinical variables with p-values
create_enhanced_forest_plot <- function(cox_model, title = "Forest Plot") {
  # Extract model results
  model_summary <- broom::tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)
  
  # Clean variable names for better display
  model_summary <- model_summary %>%
    mutate(
      term_clean = case_when(
        term == "age" ~ "Age",
        term == "genderMale" ~ "Gender (Male)",
        term == "genderfemale" ~ "Gender (Female)",
        term == "genderFemale" ~ "Gender (Female)",
        grepl("tumor_stage", term) ~ gsub("tumor_stage", "Stage ", term),
        TRUE ~ term
      )
    )
  
  # Prepare data for plotting
  plot_data <- model_summary %>%
    mutate(
      variable = term_clean,
      hr = estimate,
      lower_ci = conf.low,
      upper_ci = conf.high,
      p_value = p.value,
      significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "NS"
      )
    )
  
  # Create the forest plot
  forest_plot <- ggplot(plot_data, aes(x = hr, y = reorder(variable, hr))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_point(size = 3, shape = 18, color = "blue") +
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2, color = "blue") +
    geom_text(aes(label = paste0("HR: ", round(hr, 2), 
                                 " (", round(lower_ci, 2), "-", round(upper_ci, 2), ")\n",
                                 "p = ", format.pval(p_value, digits = 2))),
              x = max(plot_data$upper_ci, na.rm = TRUE) * 1.1, hjust = 0, size = 3.5) +
    scale_x_log10() +
    labs(title = title,
         x = "Hazard Ratio (log scale)",
         y = "Variables") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          panel.grid.major = element_line(linetype = "dotted", color = "grey80"))
  
  return(forest_plot)
}

# Generate enhanced forest plot for clinical variables (including tumor_stage)
clinical_forest_plot <- create_enhanced_forest_plot(clinical_model_for_plot, "Clinical Variables Forest Plot (Including Tumor Stage)")
ggsave("results/clinical_forest_plot_enhanced.pdf", plot = clinical_forest_plot, width = 12, height = 8, dpi = 1200)

# Also generate ggforest plot for clinical variables
tryCatch({
  clinical_ggforest <- ggforest(clinical_model_for_plot, data = final_df, fontsize = 1.0)
  ggsave("results/clinical_ggforest_plot.pdf", plot = clinical_ggforest, width = 10, height = 6, dpi = 1200)
}, error = function(e) {
  message("ggforest plot failed: ", e$message)
})

# ==========================
# FUNCTION: Run analysis for a single gene - FIXED to include tumor_stage in forest plots
# ==========================
analyze_single_gene <- function(gene_name, data) {
  
  message(paste("\n=== Analyzing gene:", gene_name, "==="))
  
  # Check if gene exists in data
  if (!gene_name %in% colnames(data)) {
    warning(paste("Gene", gene_name, "not found in data. Skipping."))
    return(NULL)
  }
  
  # Create gene-specific directory
  gene_dir <- file.path("results", gene_name)
  if (!dir.exists(gene_dir)) {
    dir.create(gene_dir)
  }
  
  # Define covariates (clinical + single gene)
  covariates <- c("tumor_stage", "age", "gender", gene_name)
  covariates <- covariates[covariates %in% names(data)]
  
  # Fit Cox model for analysis (may stratify if needed)
  formula <- as.formula(paste("Surv(os_time, os_status) ~", paste(covariates, collapse = " + ")))
  cox_model <- coxph(formula, data = data, x = TRUE, y = TRUE)
  
  # Check proportional hazards
  ph_test <- cox.zph(cox_model)
  ph_violators <- rownames(ph_test$table)[ph_test$table[, "p"] < 0.05]
  
  # Create a separate model for visualization that always includes tumor_stage as covariate
  # This ensures tumor_stage appears in forest plots even if we stratify for the actual analysis
  formula_for_plot <- as.formula(paste("Surv(os_time, os_status) ~", 
                                       paste(c("tumor_stage", "age", "gender", gene_name), collapse = " + ")))
  cox_model_for_plot <- coxph(formula_for_plot, data = data, x = TRUE, y = TRUE)
  
  # Stratify if tumor_stage violates PH (for analysis model only)
  if ("tumor_stage" %in% ph_violators) {
    message("Stratifying by tumor_stage due to PH violation")
    covariates_no_stage <- setdiff(covariates, "tumor_stage")
    formula_strat <- as.formula(paste("Surv(os_time, os_status) ~",
                                      paste(covariates_no_stage, collapse = " + "),
                                      "+ strata(tumor_stage)"))
    cox_model <- coxph(formula_strat, data = data, x = TRUE, y = TRUE)
  }
  
  # Generate risk scores using the analysis model
  data$risk_score <- predict(cox_model, type = "risk")
  risk_quantiles <- quantile(data$risk_score, probs = c(0, 0.5, 1), na.rm = TRUE)
  data$risk_group <- cut(data$risk_score,
                         breaks = risk_quantiles,
                         labels = c("Low Risk", "High Risk"),
                         include.lowest = TRUE)
  
  # Kaplan-Meier plot
  km_fit <- survfit(Surv(os_time, os_status) ~ risk_group, data = data)
  km_plot <- ggsurvplot(km_fit,
                        data = data,
                        pval = TRUE,
                        conf.int = FALSE,
                        risk.table = TRUE,
                        legend.labs = levels(data$risk_group),
                        legend.title = "Risk Group",
                        xlab = "Time (Days)",
                        ylab = "Overall Survival Probability",
                        title = paste("Prognostic Model -", gene_name),
                        subtitle = "Clinical + Single Gene",
                        palette = c("blue", "red"),
                        ggtheme = theme_bw())
  
  # Save KM plot
  ggsave(file.path(gene_dir, paste0("km_plot_", gene_name, ".pdf")), 
         plot = km_plot$plot, width = 6, height = 6, dpi = 1200)
  
  # Enhanced forest plot with p-values (using the plot model that includes tumor_stage)
  tryCatch({
    forest_plot_enhanced <- create_enhanced_forest_plot(cox_model_for_plot, paste("Forest Plot -", gene_name))
    ggsave(file.path(gene_dir, paste0("forest_plot_enhanced_", gene_name, ".pdf")), 
           plot = forest_plot_enhanced, width = 12, height = 8, dpi = 1200)
  }, error = function(e) {
    message(paste("Enhanced forest plot failed for", gene_name, ":", e$message))
  })
  
  # ggforest plot (using the plot model that includes tumor_stage)
  tryCatch({
    forest_plot_ggforest <- ggforest(cox_model_for_plot, data = data, fontsize = 0.8)
    ggsave(file.path(gene_dir, paste0("forest_plot_ggforest_", gene_name, ".pdf")), 
           plot = forest_plot_ggforest, width = 10, height = 8, dpi = 1200)
  }, error = function(e) {
    message(paste("ggforest plot failed for", gene_name, ":", e$message))
  })
  
  # Standard forestmodel plot (using the plot model that includes tumor_stage)
  tryCatch({
    forest_plot_standard <- forest_model(cox_model_for_plot)
    ggsave(file.path(gene_dir, paste0("forest_plot_standard_", gene_name, ".pdf")), 
           plot = forest_plot_standard, width = 8, height = 6, dpi = 1200)
  }, error = function(e) {
    message(paste("Standard forest plot failed for", gene_name, ":", e$message))
  })
  
  # Model summary (from the analysis model)
  results_summary <- broom::tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)
  write.csv(results_summary, file.path(gene_dir, paste0("cox_results_", gene_name, ".csv")), row.names = FALSE)
  
  # Extract clinical variable information from the plot model (which always includes all clinical variables)
  clinical_results <- broom::tidy(cox_model_for_plot, exponentiate = TRUE, conf.int = TRUE)
  
  # Extract tumor stage information
  stage_II_hr <- if (any(clinical_results$term == "tumor_stageII")) clinical_results$estimate[clinical_results$term == "tumor_stageII"] else NA
  stage_II_p <- if (any(clinical_results$term == "tumor_stageII")) clinical_results$p.value[clinical_results$term == "tumor_stageII"] else NA
  stage_III_hr <- if (any(clinical_results$term == "tumor_stageIII")) clinical_results$estimate[clinical_results$term == "tumor_stageIII"] else NA
  stage_III_p <- if (any(clinical_results$term == "tumor_stageIII")) clinical_results$p.value[clinical_results$term == "tumor_stageIII"] else NA
  stage_IV_hr <- if (any(clinical_results$term == "tumor_stageIV")) clinical_results$estimate[clinical_results$term == "tumor_stageIV"] else NA
  stage_IV_p <- if (any(clinical_results$term == "tumor_stageIV")) clinical_results$p.value[clinical_results$term == "tumor_stageIV"] else NA
  
  # Extract gender information
  gender_hr <- if (any(grepl("gender", clinical_results$term))) clinical_results$estimate[grepl("gender", clinical_results$term)] else NA
  gender_p <- if (any(grepl("gender", clinical_results$term))) clinical_results$p.value[grepl("gender", clinical_results$term)] else NA
  
  # Extract age information
  age_hr <- if (any(clinical_results$term == "age")) clinical_results$estimate[clinical_results$term == "age"] else NA
  age_p <- if (any(clinical_results$term == "age")) clinical_results$p.value[clinical_results$term == "age"] else NA
  
  # Internal validation (using the analysis model)
  c_index <- NA
  tryCatch({
    # For stratified models, use appropriate formula
    if ("strata(tumor_stage)" %in% as.character(formula(cox_model))) {
      cph_formula <- as.formula(paste("Surv(os_time, os_status) ~", 
                                      paste(setdiff(covariates, "tumor_stage"), collapse = " + "),
                                      "+ strat(tumor_stage)"))
    } else {
      cph_formula <- formula
    }
    
    cph_model <- cph(cph_formula, data = data, x = TRUE, y = TRUE, surv = TRUE)
    val_result <- validate(cph_model, method = "boot", B = 100, dxy = TRUE)
    c_index <- (val_result["Dxy", "index.corrected"] + 1) / 2
    
    # Add validation metrics to results
    results_summary$c_index <- c_index
    results_summary$gene <- gene_name
    write.csv(results_summary, file.path(gene_dir, paste0("cox_results_", gene_name, ".csv")), row.names = FALSE)
  }, error = function(e) {
    message(paste("Validation failed for", gene_name, ":", e$message))
  })
  
  # Return summary statistics including clinical variable information
  gene_result <- list(
    gene = gene_name,
    hr = if (any(results_summary$term == gene_name)) results_summary$estimate[results_summary$term == gene_name] else NA,
    hr_ci_low = if (any(results_summary$term == gene_name)) results_summary$conf.low[results_summary$term == gene_name] else NA,
    hr_ci_high = if (any(results_summary$term == gene_name)) results_summary$conf.high[results_summary$term == gene_name] else NA,
    p_value = if (any(results_summary$term == gene_name)) results_summary$p.value[results_summary$term == gene_name] else NA,
    c_index = c_index,
    # Add clinical variable information
    stage_II_hr = stage_II_hr,
    stage_II_p = stage_II_p,
    stage_III_hr = stage_III_hr,
    stage_III_p = stage_III_p,
    stage_IV_hr = stage_IV_hr,
    stage_IV_p = stage_IV_p,
    gender_hr = gender_hr,
    gender_p = gender_p,
    age_hr = age_hr,
    age_p = age_p
  )
  
  return(gene_result)
}

# ==========================
# Run analysis for each gene
# ==========================
all_results <- list()

for (gene in my_genes) {
  result <- analyze_single_gene(gene, final_df)
  if (!is.null(result)) {
    all_results[[gene]] <- result
  }
}

# ==========================
# Create summary report with clinical variable information
# ==========================
if (length(all_results) > 0) {
  summary_df <- do.call(rbind, lapply(all_results, as.data.frame))
  
  # Reorder columns to have clinical variable information at the end
  summary_df <- summary_df %>%
    dplyr::select(gene, hr, hr_ci_low, hr_ci_high, p_value, c_index, 
                  stage_II_hr, stage_II_p, stage_III_hr, stage_III_p, stage_IV_hr, stage_IV_p,
                  gender_hr, gender_p, age_hr, age_p)
  
  write.csv(summary_df, "results/gene_analysis_summary.csv", row.names = FALSE)
  
  # Create summary plot of HRs
  hr_plot <- ggplot(summary_df, aes(x = reorder(gene, hr), y = hr)) +
    geom_point(size = 3, color = "blue") +
    geom_errorbar(aes(ymin = hr_ci_low, ymax = hr_ci_high), width = 0.2, color = "blue") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    coord_flip() +
    scale_y_log10() +
    labs(title = "Hazard Ratios for Individual Genes (Adjusted for Clinical Variables)",
         x = "Gene", y = "Hazard Ratio (95% CI, log scale)") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 10, face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave("results/gene_hr_summary.pdf", plot = hr_plot, width = 12, height = 10, dpi = 1200)
  
  # Create significance plot
  sig_plot <- ggplot(summary_df, aes(x = reorder(gene, -log10(p_value)), y = -log10(p_value))) +
    geom_bar(stat = "identity", aes(fill = p_value < 0.05), width = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "blue")) +
    labs(title = "Statistical Significance of Genes",
         x = "Gene", y = "-log10(p-value)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  ggsave("results/gene_significance_plot.pdf", plot = sig_plot, width = 12, height = 8, dpi = 1200)
}

# ==========================
# Final report
# ==========================
cat("=== ANALYSIS COMPLETED ===\n")
cat("Total genes analyzed:", length(all_results), "\n")
cat("Sample size:", nrow(final_df), "\n")
cat("Events:", sum(final_df$os_status), "\n")
cat("Clinical model results saved in 'results/clinical_forest_plot_enhanced.pdf'\n")
cat("Results saved in 'results' directory\n")

if (exists("summary_df")) {
  cat("\nTop significant genes (p < 0.05):\n")
  sig_genes <- summary_df[summary_df$p_value < 0.05 & !is.na(summary_df$p_value), ]
  if (nrow(sig_genes) > 0) {
    print(sig_genes[order(sig_genes$p_value), c("gene", "hr", "p_value")])
  } else {
    cat("No significant genes found at p < 0.05\n")
  }
}
message("All analyses completed successfully!")