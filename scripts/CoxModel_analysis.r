
# ==========================
# Cox Model - OPTIMIZED with Diagnostics
# ==========================

# 0. Load essential libraries
if (!require(pacman)) install.packages("pacman")
pacman::p_load(survival, survminer, dplyr, forestmodel, ggplot2, stringr, gridExtra, timeROC, boot)

# Set seed for reproducibility
set.seed(42)

# ========================== 
# Create date-wise directory
# ==========================
analysis_date <- format(Sys.Date(), "%Y%m%d")
output_dir <- paste0("delete_Analysis_", analysis_date)
if (!dir.exists(output_dir)) dir.create(output_dir)

# Create subdirectories
subdirs <- c("UniVar", "MultiVar")
for (sub in subdirs) {
  sub_path <- file.path(output_dir, sub)
  if (!dir.exists(sub_path)) dir.create(sub_path)
}

# Function to create full output paths
output_path <- function(filename, subdir = NULL) {
  if (!is.null(subdir) && subdir %in% subdirs) {
    file.path(output_dir, subdir, filename)
  } else {
    file.path(output_dir, filename)
  }
}

#====================================
# 0. Data Preparation for Analysis
#====================================

# Load data with error handling
tryCatch({
  tcga_data <- readRDS("1.BRCA.RDS")
  limma_res <- readRDS("2.limma_res_COAD.RDS")
}, error = function(e) {
  stop("Error loading RDS files: ", e$message)
})

# Define gene list
my_genes <- c("RRM2","MYBL2","TPX2","DES","KIF18B","SYNM","SLC7A11",
              "TNS4","FOLR1","ABCG2","SORBS1","MASP1","CAV1","MYL9",
              "FHL1","RXRG","SKA3","CRYAB","HJURP","OTX1","IQGAP3","RHPN1",
              "TCEAL2","SIM2","PYGM","NRG1")

# Clean clinical data - COMBINING STAGE III & IV
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
    # Combined stage categorization - III & IV into Advanced
    tumor_stage_combined = case_when(
      str_detect(tumor_stage, "^Stage I$|^Stage IA$|^Stage IB$|^Stage IC$") ~ "I",
      str_detect(tumor_stage, "^Stage II$|^Stage IIA$|^Stage IIB$|^Stage IIC$") ~ "II", 
      str_detect(tumor_stage, "^Stage III$|^Stage IIIA$|^Stage IIIB$|^Stage IIIC$|^Stage IV$|^Stage IVA$|^Stage IVB$|^Stage IVC$") ~ "Advanced",
      TRUE ~ NA_character_
    ),
    age = as.numeric(age),
    gender = factor(gender),
    tumor_stage_combined = factor(tumor_stage_combined, levels = c("I", "II", "Advanced")),
    # Keep original for comparison
    tumor_stage_original = case_when(
      str_detect(tumor_stage, "^Stage I$|^Stage IA$|^Stage IB$|^Stage IC$") ~ "I",
      str_detect(tumor_stage, "^Stage II$|^Stage IIA$|^Stage IIB$|^Stage IIC$") ~ "II", 
      str_detect(tumor_stage, "^Stage III$|^Stage IIIA$|^Stage IIIB$|^Stage IIIC$") ~ "III",
      str_detect(tumor_stage, "^Stage IV$|^Stage IVA$|^Stage IVB$|^Stage IVC$") ~ "IV",
      TRUE ~ NA_character_
    ),
    tumor_stage_original = factor(tumor_stage_original, levels = c("I", "II", "III", "IV"))
  ) %>%
  filter(complete.cases(age, gender, tumor_stage_combined))

# DIAGNOSTIC: Check stage distribution
message("Original stage distribution:")
print(table(clinical_clean$tumor_stage_original, useNA = "always"))
message("Combined stage distribution (III+IV as Advanced):")
print(table(clinical_clean$tumor_stage_combined, useNA = "always"))

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
  clinical_clean[common_patients, c("os_time", "os_status", "age", "gender", "tumor_stage_combined", "tumor_stage_original")],
  gene_expression_data_t[common_patients, , drop = FALSE]
)

# Use combined stage as main variable
final_df$tumor_stage <- final_df$tumor_stage_combined

# ==========================
# 1. Enhanced Data validation
# ==========================
required_cols <- c("os_time", "os_status", "tumor_stage")
if (!all(required_cols %in% names(final_df))) {
  stop("Missing required columns: ", paste(setdiff(required_cols, names(final_df)), collapse = ", "))
}

final_df <- final_df %>%
  mutate(
    os_status = as.numeric(os_status),
    tumor_stage = factor(tumor_stage, levels = c("I", "II", "Advanced")),
    tumor_stage_original = factor(tumor_stage_original, levels = c("I", "II", "III", "IV"))
  )

# Remove patients with extreme survival times
final_df <- final_df %>% 
  filter(os_time > 0 & os_time <= quantile(os_time, 0.95, na.rm = TRUE))

# DIAGNOSTIC: Check final dataset
message("Final dataset summary:")
message("Number of patients: ", nrow(final_df))
message("Combined stage distribution (III+IV as Advanced):")
print(table(final_df$tumor_stage))
message("Original stage distribution:")
print(table(final_df$tumor_stage_original))
message("Event rate: ", round(mean(final_df$os_status), 3))

# ==========================
# AGE CATEGORIES AND GENDER INTEGRATION
# ==========================
message("\n=== CREATING THREE AGE CATEGORIES AND GENDER INTEGRATION ===")

# Convert age from days to years
final_df$age_years <- final_df$age / 365.25

# Create three clinically meaningful age categories
final_df$age_category <- cut(final_df$age_years,
                             breaks = c(0, 40, 60, 120),
                             labels = c("Younger (<40 years)", "Middle (40-59 years)", "Older (≥60 years)"),
                             include.lowest = TRUE)

message("Three age categories created:")
print(table(final_df$age_category))

# Check gender distribution
if ("gender" %in% names(final_df)) {
  message("Gender distribution:")
  print(table(final_df$gender))
}

# ==========================
# 2. Define covariates - UPDATED WITH AGE CATEGORY
# ==========================
gene_signature <- intersect(my_genes, names(final_df))
covariates <- c("tumor_stage", "age_category", "gender", gene_signature)
covariates <- covariates[covariates %in% names(final_df)]

message("Final covariates (with age categories): ", paste(covariates, collapse = ", "))

# ==========================
# 3. UNIVARIATE ANALYSIS FOR CLINICAL VARIABLES
# ==========================
message("\n=== UNIVARIATE ANALYSIS (WITH AGE CATEGORIES) ===")

# Univariate Cox models for each clinical variable
univariate_models <- list()

# Age Categories
if ("age_category" %in% names(final_df)) {
  uni_age_cat <- coxph(Surv(os_time, os_status) ~ age_category, data = final_df)
  univariate_models[["age_category"]] <- uni_age_cat
  age_cat_summary <- summary(uni_age_cat)
  message("Age Categories model fitted")
  for(i in 1:nrow(age_cat_summary$coefficients)) {
    cat_name <- rownames(age_cat_summary$coefficients)[i]
    hr <- round(exp(age_cat_summary$coefficients[i,1]), 3)
    pval <- round(age_cat_summary$coefficients[i,5], 4)
    message("  ", cat_name, ": HR = ", hr, ", p = ", pval)
  }
}

# Gender
if ("gender" %in% names(final_df)) {
  gender_levels <- levels(final_df$gender)
  if (length(gender_levels) >= 2) {
    uni_gender <- coxph(Surv(os_time, os_status) ~ gender, data = final_df)
    univariate_models[["gender"]] <- uni_gender
    
    gender_summary <- summary(uni_gender)
    hr <- round(exp(coef(uni_gender)), 3)
    p_value <- round(gender_summary$coefficients[1,5], 4)
    
    message("Gender - HR: ", hr, ", p-value: ", p_value)
    message("  Reference: ", gender_levels[1])
    message("  Comparison: ", gender_levels[2], " vs ", gender_levels[1])
  } else {
    message("Gender has only 1 level (", paste(gender_levels, collapse = ", "), "), skipping univariate analysis")
  }
}

# Tumor Stage (Combined: I, II, Advanced)
if ("tumor_stage" %in% names(final_df)) {
  uni_stage <- coxph(Surv(os_time, os_status) ~ tumor_stage, data = final_df)
  univariate_models[["tumor_stage"]] <- uni_stage
  message("Tumor Stage (I/II/Advanced) model fitted")
  
  stage_summary <- summary(uni_stage)
  message("Stage II vs I - HR: ", round(exp(stage_summary$coefficients[1,1]), 3),
          ", p: ", round(stage_summary$coefficients[1,5], 4))
  message("Stage Advanced vs I - HR: ", round(exp(stage_summary$coefficients[2,1]), 3),
          ", p: ", round(stage_summary$coefficients[2,5], 4))
}

# ==========================
# 4. COMPARISON: Original vs Combined Stage
# ==========================
message("\n=== STAGE GROUPING COMPARISON ===")

if ("tumor_stage_original" %in% names(final_df)) {
  uni_stage_original <- coxph(Surv(os_time, os_status) ~ tumor_stage_original, data = final_df)
  original_summary <- summary(uni_stage_original)
  
  message("Original stage model (I/II/III/IV):")
  for(i in 1:nrow(original_summary$coefficients)) {
    stage_name <- rownames(original_summary$coefficients)[i]
    hr <- round(exp(original_summary$coefficients[i,1]), 3)
    pval <- round(original_summary$coefficients[i,5], 4)
    message("  ", stage_name, ": HR = ", hr, ", p = ", pval)
  }
}

# ==========================
# 5. CLINICAL VARIABLES FOREST PLOT
# ==========================
message("\n=== CREATING CLINICAL VARIABLES FOREST PLOT ===")

clinical_vars <- c("age_category", "tumor_stage", "gender")
clinical_vars <- clinical_vars[clinical_vars %in% names(final_df)]

message("Clinical variables for forest plot: ", paste(clinical_vars, collapse = ", "))

if (length(clinical_vars) > 0) {
  formula_clinical <- as.formula(paste("Surv(os_time, os_status) ~", 
                                       paste(clinical_vars, collapse = " + ")))
  cox_clinical <- coxph(formula_clinical, data = final_df)
  
  tryCatch({
    forest_clinical <- forest_model(cox_clinical)
    ggsave(output_path("forest_plot_clinical.pdf", "UniVar"), plot = forest_clinical, 
           width = 8, height = 6, dpi = 1200)
    message("Clinical variables forest plot saved")
  }, error = function(e) {
    warning("Clinical forest plot failed: ", e$message)
  })
}

# ==========================
# 6. KAPLAN-MEIER PLOTS FOR CLINICAL VARIABLES
# ==========================
message("\n=== CREATING KAPLAN-MEIER PLOTS ===")

# Create consistent color palettes
binary_palette <- c("#E41A1C", "#377EB8")
stage_palette <- c("#1B9E77", "#D95F02", "#7570B3")
stage_palette_original <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")
age_category_palette <- c("#66C2A5", "#FC8D62", "#8DA0CB")

# KM plot for Age Categories
km_age_cat <- survfit(Surv(os_time, os_status) ~ age_category, data = final_df)
km_plot_age_cat <- ggsurvplot(km_age_cat,
                              data = final_df,
                              pval = TRUE,
                              pval.coord = c(0.1, 0.1),
                              conf.int = FALSE,
                              risk.table = TRUE,
                              risk.table.height = 0.25,
                              legend.labs = levels(final_df$age_category),
                              legend.title = "Age Category",
                              xlab = "Time (Days)",
                              ylab = "Overall Survival Probability",
                              title = "Survival by Age Category",
                              palette = age_category_palette,
                              ggtheme = theme_bw())

# KM plot for Gender
if ("gender" %in% names(final_df) && length(levels(final_df$gender)) >= 2) {
  km_gender <- survfit(Surv(os_time, os_status) ~ gender, data = final_df)
  km_plot_gender <- ggsurvplot(km_gender,
                               data = final_df,
                               pval = TRUE,
                               pval.coord = c(0.1, 0.1),
                               conf.int = FALSE,
                               risk.table = TRUE,
                               risk.table.height = 0.25,
                               legend.labs = levels(final_df$gender),
                               legend.title = "Gender",
                               xlab = "Time (Days)",
                               ylab = "Overall Survival Probability", 
                               title = "Survival by Gender",
                               palette = binary_palette,
                               ggtheme = theme_bw())
}

# KM plot for Combined Tumor Stage
km_stage <- survfit(Surv(os_time, os_status) ~ tumor_stage, data = final_df)
km_plot_stage <- ggsurvplot(km_stage,
                            data = final_df,
                            pval = TRUE,
                            pval.coord = c(0.1, 0.1),
                            conf.int = FALSE,
                            risk.table = TRUE,
                            risk.table.height = 0.25,
                            legend.labs = levels(final_df$tumor_stage),
                            legend.title = "Tumor Stage",
                            xlab = "Time (Days)",
                            ylab = "Overall Survival Probability",
                            title = "Survival by Tumor Stage (I/II/Advanced)",
                            palette = stage_palette,
                            ggtheme = theme_bw())

# Save KM plots
ggsave(output_path("km_plot_age_categories.pdf", "UniVar"), plot = km_plot_age_cat$plot, width = 8, height = 6, dpi = 1200)
if (exists("km_plot_gender")) {
  ggsave(output_path("km_plot_gender.pdf", "UniVar"), plot = km_plot_gender$plot, width = 7, height = 6, dpi = 1200)
}
ggsave(output_path("km_plot_stage_combined.pdf", "UniVar"), plot = km_plot_stage$plot, width = 7, height = 6, dpi = 1200)

# ==========================
# 7. MULTIVARIATE COX MODEL
# ==========================
message("\n=== MULTIVARIATE ANALYSIS ===")

# Remove variables with insufficient levels
valid_covariates <- c()
for (covar in covariates) {
  if (covar %in% c("age_years", "risk_score", "linear_predictor")) {
    next
  } else if (covar == "age_category") {
    valid_covariates <- c(valid_covariates, covar)
  } else if (length(unique(final_df[[covar]])) >= 2) {
    valid_covariates <- c(valid_covariates, covar)
  } else {
    message("Removing ", covar, " from multivariate model - insufficient levels")
  }
}

if (length(valid_covariates) == 0) {
  stop("No valid covariates with sufficient levels for multivariate analysis")
}

message("Final covariates for multivariate model: ", paste(valid_covariates, collapse = ", "))

# Build multivariate model
formula <- as.formula(paste("Surv(os_time, os_status) ~", paste(valid_covariates, collapse = " + ")))
cox_model <- tryCatch({
  coxph(formula, data = final_df, x = TRUE, y = TRUE)
}, error = function(e) {
  stop("Multivariate Cox model failed: ", e$message)
})

# Check proportional hazards
ph_test <- cox.zph(cox_model)
print(ph_test)

ph_violators <- rownames(ph_test$table)[ph_test$table[, "p"] < 0.05]
message("PH violators (p < 0.05): ", paste(ph_violators, collapse = ", "))

# Handle PH violations
use_stratification <- FALSE
if ("tumor_stage" %in% ph_violators) {
  message("Substantial PH violation detected for tumor_stage. Stratifying...")
  covariates_no_stage <- setdiff(valid_covariates, "tumor_stage")
  formula_strat <- as.formula(paste("Surv(os_time, os_status) ~",
                                    paste(covariates_no_stage, collapse = " + "),
                                    "+ strata(tumor_stage)"))
  cox_model <- coxph(formula_strat, data = final_df, x = TRUE, y = TRUE)
  use_stratification <- TRUE
} else {
  message("No substantial PH violations requiring stratification")
}

# ==========================
# 8. RISK SCORE CALCULATION AND STRATIFICATION
# ==========================
message("\n=== RISK STRATIFICATION ===")

# Calculate risk scores
final_df$risk_score <- predict(cox_model, type = "risk")
final_df$linear_predictor <- predict(cox_model, type = "lp")

# Create risk groups using median cutpoint
risk_quantiles <- quantile(final_df$risk_score, probs = c(0, 0.5, 1), na.rm = TRUE)

if (length(unique(risk_quantiles)) >= 2) {
  final_df$risk_group <- cut(final_df$risk_score,
                             breaks = risk_quantiles,
                             labels = c("Low Risk", "High Risk"),
                             include.lowest = TRUE)
  message("Risk groups created")
}

message("Risk group distribution:")
print(table(final_df$risk_group))

# ==========================
# 9. CLINICAL CHARACTERISTICS BY RISK GROUP
# ==========================
message("\n=== CLINICAL CHARACTERISTICS BY RISK GROUP ===")

# Function to create summary table
create_risk_group_summary <- function(data) {
  summary_list <- list()
  
  # Age Categories
  if ("age_category" %in% names(data)) {
    age_cat_summary <- data %>%
      group_by(risk_group, age_category) %>%
      tally() %>%
      group_by(risk_group) %>%
      mutate(percent = round(n / sum(n) * 100, 1))
    summary_list$age_category <- age_cat_summary
  }
  
  # Gender
  if ("gender" %in% names(data)) {
    gender_summary <- data %>%
      group_by(risk_group, gender) %>%
      tally() %>%
      group_by(risk_group) %>%
      mutate(percent = round(n / sum(n) * 100, 1))
    summary_list$gender <- gender_summary
  }
  
  # Tumor Stage
  if ("tumor_stage" %in% names(data)) {
    stage_summary <- data %>%
      group_by(risk_group, tumor_stage) %>%
      tally() %>%
      group_by(risk_group) %>%
      mutate(percent = round(n / sum(n) * 100, 1))
    summary_list$stage <- stage_summary
  }
  
  return(summary_list)
}

# Create and save risk group summary
risk_summary <- create_risk_group_summary(final_df)
if (!is.null(risk_summary$age_category)) {
  write.csv(risk_summary$age_category, output_path("risk_group_age_categories_summary.csv", "MultiVar"))
}
if (!is.null(risk_summary$gender)) {
  write.csv(risk_summary$gender, output_path("risk_group_gender_summary.csv", "MultiVar"))
}
if (!is.null(risk_summary$stage)) {
  write.csv(risk_summary$stage, output_path("risk_group_stage_summary.csv", "MultiVar"))
}

# ==========================
# 10. RISK GROUP PLOTS
# ==========================
message("\n=== CREATING RISK GROUP PLOTS ===")

# Kaplan-Meier plot for risk groups
km_risk <- survfit(Surv(os_time, os_status) ~ risk_group, data = final_df)
km_plot_risk <- ggsurvplot(km_risk,
                           data = final_df,
                           pval = TRUE,
                           pval.coord = c(0.1, 0.1),
                           conf.int = FALSE,
                           risk.table = TRUE,
                           risk.table.height = 0.25,
                           legend.labs = levels(final_df$risk_group),
                           legend.title = "Risk Group",
                           xlab = "Time (Days)",
                           ylab = "Overall Survival Probability",
                           title = "Survival by Risk Group",
                           subtitle = "Clinical + Gene Signature (Stages III+IV as Advanced)",
                           palette = c("blue", "red"),
                           ggtheme = theme_bw())

ggsave(output_path("km_plot_risk_groups.pdf", "MultiVar"), plot = km_plot_risk$plot, 
       width = 7, height = 6, dpi = 1200)

# Age category distribution by risk group
age_cat_risk_plot <- ggplot(final_df, aes(x = age_category, fill = risk_group)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("blue", "red")) +
  labs(title = "Age Category Distribution by Risk Group",
       x = "Age Category", y = "Proportion", fill = "Risk Group") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(output_path("age_three_categories_risk_groups.pdf", "MultiVar"), 
       plot = age_cat_risk_plot, width = 8, height = 6, dpi = 1200)

# ==========================
# 11. BOOTSTRAP VALIDATION OF C-INDEX
# ==========================
message("\n=== BOOTSTRAP VALIDATION OF C-INDEX ===")

# Function to calculate C-index
calc_c_index <- function(data, indices) {
  bootstrap_data <- data[indices, ]
  
  tryCatch({
    formula <- as.formula(paste("Surv(os_time, os_status) ~", paste(valid_covariates, collapse = " + ")))
    cox_boot <- coxph(formula, data = bootstrap_data)
    
    # Calculate C-index on bootstrap sample
    c_apparent <- summary(cox_boot)$concordance[1]
    
    # Calculate C-index on original sample
    predictions <- predict(cox_boot, newdata = data)
    c_test <- survConcordance(Surv(os_time, os_status) ~ predictions, data = data)$concordance
    
    return(c(c_apparent, c_test))
  }, error = function(e) {
    return(c(NA, NA))
  })
}

# Perform bootstrap validation
set.seed(42)
n_boot <- 1000
boot_results <- boot::boot(final_df, calc_c_index, R = n_boot)

# Calculate optimism
optimism <- boot_results$t[,1] - boot_results$t[,2]
mean_optimism <- mean(optimism, na.rm = TRUE)

# Get original C-index
original_c_index <- summary(cox_model)$concordance[1]
corrected_c_index <- original_c_index - mean_optimism

message("Bootstrap Validation Results (n=", n_boot, "):")
message("Original C-index: ", round(original_c_index, 3))
message("Mean Optimism: ", round(mean_optimism, 3))
message("Optimism-corrected C-index: ", round(corrected_c_index, 3))

# Save bootstrap results
bootstrap_validation <- data.frame(
  Metric = c("Original_C_index", "Bootstrap_samples", "Mean_Optimism", "Corrected_C_index"),
  Value = c(original_c_index, n_boot, mean_optimism, corrected_c_index)
)
write.csv(bootstrap_validation, output_path("bootstrap_validation.csv", "MultiVar"), row.names = FALSE)

# ==========================
# 12. ROC CURVE ANALYSIS
# ==========================
message("\n=== ROC CURVE ANALYSIS ===")

# Define time points for time-dependent ROC
if (sum(final_df$os_status) > 0) {
  times <- quantile(final_df$os_time[final_df$os_status == 1], 
                    probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
} else {
  times <- quantile(final_df$os_time, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
}

message("Time-dependent ROC evaluation at: ", paste(round(times), "days"))

# Time-dependent ROC
tryCatch({
  roc_data <- timeROC(T = final_df$os_time,
                      delta = final_df$os_status,
                      marker = final_df$linear_predictor,
                      cause = 1,
                      times = times,
                      ROC = TRUE)
  
  # Combined ROC plot
  pdf(output_path("time_dependent_roc_combined.pdf", "MultiVar"), width = 8, height = 8)
  plot(roc_data, time = times[1], title = FALSE, col = "#1B9E77", lwd = 3)
  plot(roc_data, time = times[2], col = "#D95F02", lwd = 3, add = TRUE)
  plot(roc_data, time = times[3], col = "#7570B3", lwd = 3, add = TRUE)
  
  abline(a = 0, b = 1, col = "gray50", lty = 2, lwd = 1)
  
  legend("bottomright", 
         legend = c(paste0(round(times[1]), " days (AUC = ", round(roc_data$AUC[1], 3), ")"),
                    paste0(round(times[2]), " days (AUC = ", round(roc_data$AUC[2], 3), ")"),
                    paste0(round(times[3]), " days (AUC = ", round(roc_data$AUC[3], 3), ")"),
                    "Reference"),
         col = c("#1B9E77", "#D95F02", "#7570B3", "gray50"),
         lty = c(1, 1, 1, 2),
         lwd = c(3, 3, 3, 1),
         bty = "n",
         cex = 0.9)
  
  title(main = "Time-Dependent ROC Curves\nClinical + Gene Signature Model",
        xlab = "1 - Specificity (False Positive Rate)",
        ylab = "Sensitivity (True Positive Rate)")
  
  grid(col = "gray80", lty = 3)
  dev.off()
  
  message("Time-dependent AUC values:")
  for(i in 1:length(times)) {
    message("  ", round(times[i]), " days: AUC = ", round(roc_data$AUC[i], 3))
  }
  
}, error = function(e) {
  warning("Time-dependent ROC calculation failed: ", e$message)
})

# ==========================
# 13. FOREST PLOTS
# ==========================
message("\n=== CREATING FOREST PLOTS ===")

# Clinical-only forest plot
clinical_vars <- c("tumor_stage", "age_category", "gender")
clinical_vars <- clinical_vars[clinical_vars %in% names(final_df)]

if (length(clinical_vars) > 0) {
  formula_clinical <- as.formula(paste("Surv(os_time, os_status) ~", 
                                       paste(clinical_vars, collapse = " + ")))
  cox_clinical <- coxph(formula_clinical, data = final_df, x = TRUE)
  
  forest_clinical <- forest_model(cox_clinical)
  ggsave(output_path("forest_plot_clinical_only.pdf", "MultiVar"), plot = forest_clinical, 
         width = 8, height = 6, dpi = 1200)
  message("Clinical-only forest plot saved")
}

# Full model forest plot
if (use_stratification) {
  formula_for_plot <- as.formula(paste("Surv(os_time, os_status) ~", 
                                       paste(valid_covariates, collapse = " + ")))
  cox_model_for_plot <- coxph(formula_for_plot, data = final_df, x = TRUE)
} else {
  cox_model_for_plot <- cox_model
}

forest_full <- forest_model(cox_model_for_plot)
ggsave(output_path("forest_plot_full_model.pdf","MultiVar"), plot = forest_full, 
       width = 10, height = 12, dpi = 1200)
message("Full model forest plot saved")

# ==========================
# 14. GLOBAL LOG-RANK P-VALUE CALCULATION
# ==========================
message("\n=== GLOBAL LOG-RANK P-VALUES ===")

clinical_vars <- c("age_category", "tumor_stage", "gender", "risk_group")

for (var in clinical_vars) {
  if (var %in% names(final_df)) {
    formula <- as.formula(paste("Surv(os_time, os_status) ~", var))
    logrank_test <- survdiff(formula, data = final_df)
    global_p <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
    message("Global log-rank P-value (", var, "): ", signif(global_p, 4))
  }
}

# ==========================
# 15. FINAL SUMMARY AND SAVE
# ==========================
# Save complete analysis objects
saveRDS(list(
  cox_model = cox_model,
  final_df = final_df,
  c_index_multi = original_c_index,
  c_index_corrected = corrected_c_index,
  roc_data = if(exists("roc_data")) roc_data else NULL,
  risk_summary = risk_summary
), output_path("analysis_results.rds","MultiVar"))

# Save session info
writeLines(capture.output(sessionInfo()), output_path("session_info.txt","MultiVar"))

message("\n=== ANALYSIS COMPLETED SUCCESSFULLY ===")
message("All results saved in: ", output_dir)
message("\nKEY OUTPUTS GENERATED:")
message("AGE CATEGORIES:")
message("  - Three age categories: <40, 40-59, ≥60 years")
message("  - km_plot_age_categories.pdf: KM plot with three age groups")
message("STAGE GROUPING:")
message("  - Stages III and IV combined into 'Advanced' stage")
message("  - km_plot_stage_combined.pdf: KM plot with I/II/Advanced")
message("RISK GROUP PLOTS:")
message("  - km_plot_risk_groups.pdf: KM plot for risk groups")
message("  - age_three_categories_risk_groups.pdf: Age category distribution by risk")
message("PERFORMANCE METRICS:")
message("  - bootstrap_validation.csv: C-index validation results")
message("  - time_dependent_roc_combined.pdf: ROC curves")
message("FOREST PLOTS:")
message("  - forest_plot_clinical_only.pdf: Clinical variables only")
message("  - forest_plot_full_model.pdf: Clinical + genes")
