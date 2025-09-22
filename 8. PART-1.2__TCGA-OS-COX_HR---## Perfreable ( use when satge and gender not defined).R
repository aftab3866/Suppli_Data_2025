# ==========================
# Publication-ready Cox Model - OPTIMIZED
# ==========================

# 0. Load libraries with pacman for efficient loading
if (!require(pacman)) install.packages("pacman")
pacman::p_load(survival, survminer, dplyr, forestmodel, rms, ggplot2, stringr, gridExtra)

# Set seed for reproducibility
set.seed(42)

# ==========================
# Create date-wise directory
# ==========================
analysis_date <- format(Sys.Date(), "%Y%m%d")
output_dir <- paste0("Analysis_", analysis_date)

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  message("Created directory: ", output_dir)
} else {
  message("Directory already exists: ", output_dir)
}

# Function to create full output paths
output_path <- function(filename) {
  file.path(output_dir, filename)
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
my_genes <- c("RRM2", "MYBL2", "TPX2", "DES", "KIF18B", "SYNM", "SLC7A11", "TNS4", 
              "FOLR1", "ABCG2", "SORBS1", "MASP1", "CAV1", "MYL9", "FHL1", "RXRG", 
              "SKA3", "CRYAB", "HJURP", "OTX1", "IQGAP3", "RHPN1", "TCEAL2", "SIM2", 
              "PYGM", "NRG1")

# Clean clinical data - REMOVED GENDER AND TUMOR_STAGE
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
    age = age_at_diagnosis
    # REMOVED: gender, tumor_stage
  ) %>%
  mutate(
    age = as.numeric(age)
  ) %>%
  filter(complete.cases(age)) # REMOVED: gender, tumor_stage_simple from complete.cases

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
  clinical_clean[common_patients, c("os_time", "os_status", "age")], # REMOVED: gender, tumor_stage_simple
  gene_expression_data_t[common_patients, , drop = FALSE]
)

# ==========================
# 1. Data validation
# ==========================
required_cols <- c("os_time", "os_status", "age") # REMOVED: tumor_stage
if (!all(required_cols %in% names(final_df))) {
  stop("Missing required columns: ", paste(setdiff(required_cols, names(final_df)), collapse = ", "))
}

final_df <- final_df %>%
  mutate(
    os_status = as.numeric(os_status)
    # REMOVED: tumor_stage factor conversion
  )

# Remove patients with extreme survival times
final_df <- final_df %>% filter(os_time > 0 & os_time <= quantile(os_time, 0.99, na.rm = TRUE))

# ==========================
# 2. Define covariates - ONLY AGE + GENES
# ==========================
gene_signature <- intersect(my_genes, names(final_df))
covariates <- c("age", gene_signature) # REMOVED: tumor_stage, gender
covariates <- covariates[covariates %in% names(final_df)]

# ==========================
# 3. Fit Cox model with PH check
# ==========================
formula <- as.formula(paste("Surv(os_time, os_status) ~", paste(covariates, collapse = " + ")))
cox_model <- coxph(formula, data = final_df, x = TRUE, y = TRUE)

# Check proportional hazards
ph_test <- cox.zph(cox_model)
ph_violators <- rownames(ph_test$table)[ph_test$table[, "p"] < 0.05]

# REMOVED: Stratification logic since tumor_stage is removed

# ==========================
# 4. Generate risk scores
# ==========================
final_df$risk_score <- predict(cox_model, type = "risk")
final_df$risk_group <- cut(final_df$risk_score,
                           breaks = quantile(final_df$risk_score, probs = c(0, 0.5, 1), na.rm = TRUE),
                           labels = c("Low Risk", "High Risk"),
                           include.lowest = TRUE)

#===========================
# Create consistent color palettes
binary_palette <- c("#E41A1C", "#377EB8")  # Red, Blue (colorblind-friendly)

# ==========================
# 5a. Kaplan-Meier curve for Risk Groups (Combined Model)
# ==========================
km_fit <- survfit(Surv(os_time, os_status) ~ risk_group, data = final_df)
km_plot <- ggsurvplot(km_fit,
                      data = final_df,
                      pval = TRUE,
                      conf.int = F,
                      risk.table = TRUE,
                      legend.labs = c("Low Risk","High Risk"),
                      legend.title = "Risk Group",
                      xlab = "Time (Days)",
                      ylab = "Overall Survival Probability",
                      title = "Prognostic Model",
                      subtitle = "Age + Gene Signature", # UPDATED subtitle
                      palette = c("blue", "red"),
                      ggtheme = theme_bw())
print(km_plot)
ggsave(output_path("kaplan_meier_plot.pdf"), plot = km_plot$plot, width = 5, height = 5, dpi = 1200)

# ==========================
# 5b. KM plots for individual variables - ONLY AGE REMAINS
# ==========================

# Create age groups (median split)
final_df$age_group <- cut(final_df$age,
                          breaks = quantile(final_df$age, probs = c(0, 0.5, 1), na.rm = TRUE),
                          labels = c("Older","Younger"),
                          include.lowest = TRUE)

# KM plot for Age
km_age <- survfit(Surv(os_time, os_status) ~ age_group, data = final_df)
km_plot_age <- ggsurvplot(km_age,
                          data = final_df,
                          pval = TRUE,
                          conf.int = F,
                          risk.table = TRUE,
                          legend.labs = levels(final_df$age_group),
                          legend.title = "Age Group",
                          xlab = "Time (Days)",
                          ylab = "Overall Survival Probability",
                          title = "Survival by Age",
                          palette = binary_palette,
                          ggtheme = theme_bw())
print(km_plot_age)
ggsave(output_path("km_plot_age.pdf"), plot = km_plot_age$plot, width = 6, height = 6, dpi = 1200)

# REMOVED: KM plots for Gender and Tumor Stage

# Combined plot with remaining variables
combined_plots <- grid.arrange(
  km_plot$plot + labs(title = "Risk Group"),
  km_plot_age$plot + labs(title = "Age"),
  ncol = 2
)
ggsave(output_path("km_plots_combined.pdf"), plot = combined_plots, width = 12, height = 6, dpi = 1200)

# ==========================
# 6. Forest plot
# ==========================
forest_plot <- forest_model(cox_model)
print(forest_plot)
ggsave(output_path("cox_forest_plot.pdf"), plot = forest_plot, width = 8, height = 10, dpi = 1200)

# Alternative forest plot with ggforest
forest_plot2 <- ggforest(cox_model, data = final_df, fontsize = 1)
print(forest_plot2)
ggsave(output_path("cox_forest_plot_ggforest.pdf"), plot = forest_plot2, width = 10, height = 11)

# ==========================
# 7. Internal validation
# ==========================
dd <- datadist(final_df)
options(datadist = "dd")

cph_model <- cph(formula, data = final_df, x = TRUE, y = TRUE, surv = TRUE)
val_result <- validate(cph_model, method = "boot", B = 200, dxy = TRUE)

# ==========================
# 8. Save results efficiently
# ==========================
results_summary <- broom::tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)
write.csv(results_summary, output_path("final_cox_results.csv"), row.names = FALSE)

# Save minimal necessary objects
saveRDS(list(cox_model = cox_model, final_df = final_df), output_path("analysis_results.rds"))

# ==========================
# 9. Generate comprehensive report
# ==========================
cat("=== ANALYSIS COMPLETED ===\n")
cat("Sample size:", nrow(final_df), "\n")
cat("Events:", sum(final_df$os_status), "\n")
cat("Median follow-up:", median(final_df$os_time), "days\n")
cat("Clinical variables included: age only\n") # UPDATED
cat("Gene signature size:", length(gene_signature), "genes\n")
cat("PH violations:", ifelse(length(ph_violators) > 0, paste(ph_violators, collapse = ", "), "None"), "\n")
cat("Validation Dxy:", round(val_result["Dxy", "index.corrected"], 3), "\n")

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), output_path("session_info.txt"))

message("All analyses completed successfully! Results saved in: ", output_dir)