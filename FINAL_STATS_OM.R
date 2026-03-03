# ==============================================================================
# ULTIMATE LMM SCRIPT: EEG LMM Analysis (Baseline, Longitudinal, Clinical)
# ==============================================================================

#### - STEP 1. Load Libraries & Setup ####
libraries <- c('pacman', 'MHTdiscrete', 'dplyr', 'readr', 'readxl', 'tidyr', 'rstudioapi',
               'ggplot2', 'stringr', 'BayesFactor', 'emmeans', 'boot', 'boot.pval', 
               'olsrr', 'sjPlot', 'car', 'scales', 'ggeffects', 'lme4', 'lmerTest', 'broom.mixed', "purrr",
               "ggdist", "glue")



do.call(pacman::p_load, as.list(libraries))
options(warn = -1, scipen = 999)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(123)  # For reproducibility


#### - STEP 2. Load and Prepare the Data ####
eeg_data  <- read_csv("../data/Data/Data Tables/eeg_data_table_091225.csv")
metatable <- read_csv("../data/Data/Metatable/metatable_081225.csv") %>%
  mutate(across(c(PANSS_pos, PANSS_neg, PANSS_total, PANSS_gen), ~ ifelse(.x == 0, NA, .x)))

# Genc's Delta Dataframe
genc_df <- read_csv("/Users/genchasanaj/Library/CloudStorage/GoogleDrive-gencxhasanaj@gmail.com/My Drive/LMU Klinikum/Projects/OscillationInMotion/data/Data/Data Tables/GencTable.csv") %>%
  rename_with(~ make.names(.x)) %>%
  mutate(
    Condition = factor(
      case_when(
        Condition == "eyes-closed" ~ "EC",
        Condition == "eyes-open"   ~ "EO",
        TRUE ~ as.character(Condition)
      )
    ),
    sex = factor(
      case_when(
        sex == "m" ~ "Male",
        sex == "f" ~ "Female",
        TRUE ~ as.character(sex)
      )
    )
  )



# Merge and Clean Master Data
dataframe <- eeg_data %>%
  left_join(metatable, by = c("ID", "Session")) %>%
  mutate(
    Session = as.factor(Session),
    Condition = dplyr::case_when(
      Condition == "eyes-closed" ~ "EC",
      Condition == "eyes-open"   ~ "EO",
      TRUE                       ~ as.character(Condition)),
      sex = factor(
        case_when(
          sex == "m" ~ "Male",
          sex == "f" ~ "Female",
          TRUE ~ as.character(sex)
        )
    )
  ) %>%
  filter((excluded_ec == 0 | excluded_eo == 0), excluded == 0) %>%
  rename_with(~ make.names(.x)) %>%
  filter(grepl("BHC|SCZ", ID))

# Subsets for specific analyses
dataframe_baseline <- dataframe %>% filter(Session == "V1")

dataframe_deniz <- dataframe %>% filter(Group_SSD == "SSD") %>%
  mutate(cpz_eq_bl = if_else(Session == "V1", CPZ_eq, NA_real_),
         age_bl = ifelse(Session == "V1", age, NA_real_),
         sex_bl = ifelse(Session == "V1", sex, NA_real_)) %>%
  group_by(ID) %>%
  fill(cpz_eq_bl, .direction = "downup") %>% 
  fill(age_bl, .direction = "downup") %>%
  fill(sex_bl, .direction = "downup") %>%
  ungroup() %>%
  filter(grepl("SCZ", ID))


#### - STEP 3. The Ultimate LMM Pipeline ####
run_full_lmm_pipeline <- function(data, outcomes, predictors, predictor_name = "Main",
                                  covariates = NULL, random_effects = "(1|ID)", 
                                  loop_vars = c("Region", "Periodic_param"),
                                  use_reml = TRUE,                   
                                  interaction_term = NULL,           
                                  interaction_type = "categorical", 
                                  condition_col = "Condition",      
                                  continuous_var = NULL,             
                                  id_col = "ID",                     
                                  save_dir = "../results/BaselineLmmOutputs") {
  
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  loop_combinations <- data %>% 
    dplyr::select(dplyr::all_of(loop_vars)) %>% 
    dplyr::distinct() %>%
    tidyr::drop_na() 
  
  all_main_results <- list()
  all_posthoc_results <- list()
  all_models_list <- list() 
  
  for (current_outcome in outcomes) {
    message(paste("\n--- Processing:", current_outcome, "| Predictor:", predictor_name, "---"))
    
    for (i in 1:nrow(loop_combinations)) {
      current_combo <- loop_combinations[i, , drop = FALSE]
      key_elements <- c(as.character(unname(unlist(current_combo))), current_outcome, predictor_name)
      model_key <- paste(key_elements, collapse = "__")
      
      subset_data <- data %>% dplyr::inner_join(current_combo, by = loop_vars)
      
      # FIX: Added backticks around outcome for formula safety
      formula_str <- paste0("`", current_outcome, "` ~ ", predictors)
      if (!is.null(covariates)) formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
      formula_str <- paste(formula_str, "+", random_effects)
      
      model_fit <- tryCatch({
        lmerTest::lmer(as.formula(formula_str), data = subset_data, REML = use_reml)
      }, error = function(e) {
        message(paste("  -> SKIPPED MODEL:", model_key, "| Error:", e$message))
        return(NULL)
      })
      
      if (!is.null(model_fit)) {
        all_models_list[[model_key]] <- model_fit 
        n_obs <- nobs(model_fit)
        n_id <- as.numeric(summary(model_fit)$ngrps[id_col])
        
        # 4. MAIN EFFECTS
        tidy_res <- tryCatch({
          broom.mixed::tidy(model_fit, conf.int = TRUE, conf.method = "profile")
        }, error = function(e) {
          broom.mixed::tidy(model_fit, conf.int = TRUE, conf.method = "Wald")
        })
        
        tidy_res <- tidy_res %>%
          dplyr::filter(effect == "fixed") %>% 
          dplyr::mutate(Outcome = current_outcome, Predictor_Set = predictor_name,
                        N_unique_ID = n_id, N_observations = n_obs)
        
        for (var in loop_vars) tidy_res[[var]] <- current_combo[[var]][1]
        
        all_main_results[[model_key]] <- tidy_res %>%
          dplyr::select(dplyr::all_of(loop_vars), Outcome, Predictor_Set, Term = term, 
                        N_unique_ID, N_observations, estimate, std_error = std.error, 
                        df, t_val = statistic, p_value = p.value, low_CI = conf.low, high_CI = conf.high)
        
        # 5. POST-HOCS
        if (!is.null(interaction_term)) {
          model_data <- model.frame(model_fit) 
          
          tryCatch({
            # PATH A: CATEGORICAL
            if (interaction_type == "categorical" && grepl("\\|", interaction_term)) {
              lhs_var <- trimws(strsplit(interaction_term, "\\|")[[1]][1])
              rhs_var <- trimws(strsplit(interaction_term, "\\|")[[1]][2])
              
              emm <- emmeans::emmeans(model_fit, as.formula(paste("~", interaction_term)))
              tidy_ph <- broom::tidy(pairs(emm, adjust = "none"), conf.int = TRUE)
              
              stats_df <- model_data %>%
                dplyr::group_by(.data[[lhs_var]], .data[[rhs_var]]) %>%
                dplyr::summarise(
                  n = dplyr::n_distinct(.data[[id_col]]),
                  mean_val = mean(.data[[current_outcome]], na.rm = TRUE),
                  .groups = 'drop'
                ) %>%
                dplyr::mutate(
                  formatted_n = paste0(.data[[lhs_var]], " (", n, ")"),
                  formatted_mean = paste0(.data[[lhs_var]], " (", round(mean_val, 3), ")")
                ) %>%
                dplyr::group_by(.data[[rhs_var]]) %>%
                dplyr::summarise(
                  # FIX: Added backticks to column names with spaces
                  `Group (Count)` = paste(formatted_n, collapse = " vs "),
                  `Group (Mean)` = paste(formatted_mean, collapse = " vs "),
                  .groups = "drop"
                )
              
              tidy_ph <- tidy_ph %>%
                dplyr::left_join(stats_df, by = rhs_var) %>%
                dplyr::mutate(
                  Outcome = current_outcome, Predictor_Set = predictor_name,
                  Significance = dplyr::case_when(
                    p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", 
                    p.value < 0.05 ~ "*", p.value < 0.1 ~ "marginal", TRUE ~ "ns"
                  )
                )
              
              for (var in loop_vars) tidy_ph[[var]] <- current_combo[[var]][1]
              all_posthoc_results[[model_key]] <- tidy_ph %>%
                dplyr::relocate(dplyr::all_of(loop_vars), Outcome, Predictor_Set)
              
              # PATH B: CONTINUOUS
            } else if (interaction_type == "continuous" && !is.null(continuous_var)) {
              tr <- emmeans::emtrends(model_fit, specs = stats::as.formula(paste0("~ ", condition_col)), var = continuous_var)
              tr_df <- as.data.frame(summary(tr, infer = c(TRUE, TRUE)))
              
              candidates <- c("trend", continuous_var, paste0(continuous_var, ".trend"))
              slope_col <- candidates[candidates %in% names(tr_df)][1]
              
              counts_df <- model_data %>%
                dplyr::group_by(!!rlang::sym(condition_col)) %>%
                dplyr::summarise(N_unique_ID = dplyr::n_distinct(!!rlang::sym(id_col)),
                                 N_observations = dplyr::n(), .groups = "drop") %>%
                dplyr::rename(Condition = !!rlang::sym(condition_col))
              
              tidy_ph <- tr_df %>%
                dplyr::transmute(Condition = .data[[condition_col]], Slope_Estimate = .data[[slope_col]],
                                 Slope_SE = .data[["SE"]], df = .data[["df"]], t_val = .data[["t.ratio"]],
                                 p_value = .data[["p.value"]], low_CI = .data[["lower.CL"]], high_CI = .data[["upper.CL"]]) %>%
                dplyr::left_join(counts_df, by = "Condition") %>%
                dplyr::mutate(Outcome = current_outcome, Predictor_Set = predictor_name, Trend_Variable = continuous_var,
                              Significance = dplyr::case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
                                                              p_value < 0.05 ~ "*", p_value < 0.1 ~ "marginal", TRUE ~ "ns"))
              
              for (var in loop_vars) tidy_ph[[var]] <- current_combo[[var]][1]
              all_posthoc_results[[model_key]] <- tidy_ph %>%
                dplyr::select(dplyr::all_of(loop_vars), Outcome, Predictor_Set, Condition, Trend_Variable, 
                              N_unique_ID, N_observations, Slope_Estimate, Slope_SE, df, t_val, p_value, low_CI, high_CI, Significance)
            }
          }, error = function(e) {
            message("POSTHOC FAILED for ", model_key, " | ", e$message)
          })
        }
      }
    }
  }
  
  master_results <- if(length(all_main_results) > 0) dplyr::bind_rows(all_main_results) else NULL
  master_posthocs <- if(length(all_posthoc_results) > 0) dplyr::bind_rows(all_posthoc_results) else NULL
  
  return(list(main_table = master_results, posthoc_table = master_posthocs, models = all_models_list))
}


#### - STEP 4. Helper Function for Standardization ####
filter_eeg_metrics <- function(df) {
  if(is.null(df)) return(NULL)
  df %>%
    filter(
      (Outcome == "Theta.4.8Hz"      & Periodic_param == "PW") |
        (Outcome == "Alpha.8.12Hz"     & Periodic_param == "CF") |
     #   (Outcome == "Exponent_3.40Hz"  & Periodic_param == "PW") |
        (Outcome == "Exponent_30.48Hz" & Periodic_param == "PW")
    ) %>%
    mutate(to_corr1 = if_else(Region == "Global", "Global", "Regional"))
}


# ==============================================================================
# ANALYSIS 1: NINA (Baseline HC vs SSD - Categorical)
# ==============================================================================
message("\n>>> STARTING NINA ANALYSIS <<<")
Results_Nina <- list()
Results_Nina$raw_output <- run_full_lmm_pipeline(
  data = dataframe_baseline, 
  outcomes = c("Theta.4.8Hz", "Alpha.8.12Hz", "Exponent_3.40Hz", "Exponent_1.70Hz", "Exponent_30.48Hz"), 
  predictors = "Group_SSD * Condition", predictor_name = "Group_SSD",
  covariates = c("age", "sex"), random_effects = "(1|ID)", 
  loop_vars = c("Region", "Periodic_param"), 
  use_reml = TRUE,
  interaction_term = "Group_SSD | Condition", interaction_type = "categorical",
  save_dir = "../results/Nina_Outputs"
)
saveRDS(Results_Nina, file = "../results/Nina_Outputs/Results_Nina_Master.rds")




# ==============================================================================
# ANALYSIS 2: DENIZ (Longitudinal SSD - Categorical)
# ==============================================================================
message("\n>>> STARTING DENIZ ANALYSIS <<<")
Results_Deniz <- list()
Results_Deniz$raw_output <- run_full_lmm_pipeline(
  data = dataframe_deniz, 
  outcomes = c("Theta.4.8Hz", "Alpha.8.12Hz", "Exponent_3.40Hz", "Exponent_1.70Hz", "Exponent_30.48Hz"), 
  predictors = "Session * Condition", predictor_name = "Session",
  covariates = c("age_bl", "sex_bl", "cpz_eq_bl"), random_effects = "(1|ID)", 
  loop_vars = c("Region", "Periodic_param"), 
  use_reml = TRUE,
  interaction_term = "Session | Condition", interaction_type = "categorical",
  save_dir = "../results/Deniz_Outputs"
)
saveRDS(Results_Deniz, file = "../results/Deniz_Outputs/Results_Deniz_Master.rds")


# ==============================================================================
# ANALYSIS 3: BERKHAN (Baseline Clinical Predictors - Continuous Simple Slopes)
# ==============================================================================
message("\n>>> STARTING BERKHAN ANALYSIS <<<")
pred_vars <- c("BACS_composite_z", "PANSS_total", "PANSS_pos", "PANSS_neg", "PANSS_gen", "FROGS_total")
Results_Berkhan <- list(main_table = data.frame(), posthoc_table = data.frame(), models = list())

for (p in pred_vars) {
  formula_pred <- paste0(p, " * Condition")
  temp_out <- run_full_lmm_pipeline(
    data = dataframe_baseline, 
    outcomes = c("Theta.4.8Hz", "Alpha.8.12Hz", "Exponent_3.40Hz", "Exponent_1.70Hz", "Exponent_30.48Hz"), 
    predictors = formula_pred, predictor_name = p,
    covariates = c("age", "sex", "CPZ_eq"), random_effects = "(1|ID)", 
    loop_vars = c("Region", "Periodic_param"), 
    use_reml = TRUE,
    interaction_term = "Condition", interaction_type = "continuous", 
    condition_col = "Condition", continuous_var = p, 
    save_dir = "../results/Berkhan_Outputs"
  )
  Results_Berkhan$main_table <- bind_rows(Results_Berkhan$main_table, temp_out$main_table)
  Results_Berkhan$posthoc_table <- bind_rows(Results_Berkhan$posthoc_table, temp_out$posthoc_table)
  Results_Berkhan$models <- c(Results_Berkhan$models, temp_out$models)
}
saveRDS(Results_Berkhan, file = "../results/Berkhan_Outputs/Results_Berkhan_Master.rds")


# ==============================================================================
# ANALYSIS 4: GENC (Delta Clinical Predictors - Continuous Simple Slopes)
# ==============================================================================
message("\n>>> STARTING GENC ANALYSIS <<<")
Results_Genc <- list(main_table = data.frame(), posthoc_table = data.frame(), models = list())

for (p in pred_vars) {
  formula_pred <- paste0(p, " * Condition")
  temp_out <- run_full_lmm_pipeline(
    data = genc_df, 
    outcomes = c("Theta.4.8Hz","Alpha.8.12Hz", "Exponent_3.40Hz", "Exponent_1.70Hz", "Exponent_30.48Hz"), 
    predictors = formula_pred, predictor_name = p,
    covariates = c("age", "sex", "CPZ_eq_baseline"), random_effects = "(1|ID)", 
    loop_vars = c("Region", "Periodic_param"), 
    use_reml = TRUE,
    interaction_term = "Condition", interaction_type = "continuous", 
    condition_col = "Condition", continuous_var = p, 
    save_dir = "../results/Genc_Outputs"
  )
  Results_Genc$main_table <- bind_rows(Results_Genc$main_table, temp_out$main_table)
  Results_Genc$posthoc_table <- bind_rows(Results_Genc$posthoc_table, temp_out$posthoc_table)
  Results_Genc$models <- c(Results_Genc$models, temp_out$models)
}
saveRDS(Results_Genc, file = "../results/Genc_Outputs/Results_Genc_Master.rds")

message("\n✅ ALL ANALYSES COMPLETED AND SAVED SUCCESSFULLY!")








# ==============================================================================
# 0. HELPER: STANDARDIZE COLUMN NAMES
# ==============================================================================
standardize_lmm_cols <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if ("p.value" %in% names(df) && !"p_value" %in% names(df)) df <- df %>% rename(p_value = p.value)
  if ("Estimate" %in% names(df) && !"estimate" %in% names(df)) df <- df %>% rename(estimate = Estimate)
  if ("std.error" %in% names(df) && !"std_error" %in% names(df)) df <- df %>% rename(std_error = std.error)
  if ("Std. Error" %in% names(df) && !"std_error" %in% names(df)) df <- df %>% rename(std_error = `Std. Error`)
  if ("statistic" %in% names(df) && !"t_val" %in% names(df)) df <- df %>% rename(t_val = statistic)
  if ("t value" %in% names(df) && !"t_val" %in% names(df)) df <- df %>% rename(t_val = `t value`)
  if ("conf.low" %in% names(df) && !"low_CI" %in% names(df)) df <- df %>% rename(low_CI = conf.low)
  if ("conf.high" %in% names(df) && !"high_CI" %in% names(df)) df <- df %>% rename(high_CI = conf.high)
  return(df)
}

# ==============================================================================
# 1. MASTER FORMATTING: MAIN EFFECTS
# ==============================================================================
format_publication_main <- function(df, is_clinical = FALSE, is_regional = FALSE) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- standardize_lmm_cols(df)
  
  # Filter Region and Target Metrics (Mapping 30-48 strictly)
  pub_df <- df %>%
    filter(if(is_regional) Region != "Global" else Region == "Global") %>%
    filter(
      (Periodic_param == "PW" & Outcome == "Theta.4.8Hz") |
        (Periodic_param == "CF" & Outcome == "Alpha.8.12Hz") |
        (Periodic_param == "PW" & (Outcome == "Exponent_30.48Hz"))
    )
  
  if (is_clinical) {
    pub_df <- pub_df %>% 
      filter(Predictor_Set %in% c("BACS_composite_z", "PANSS_neg", "PANSS_pos", "PANSS_gen")) %>%
      mutate(Clinical = if_else(Predictor_Set == "BACS_composite_z", "Cognition", "Symptoms"))
  }
  
  if (nrow(pub_df) == 0) return(NULL)
  
  # FDR Correction
  if (!is_regional) {
    group_vars <- if(is_clinical) c("Term", "Predictor_Set") else "Term"
    pub_df <- pub_df %>% group_by(across(all_of(group_vars))) %>% 
      mutate(p_fdr = p.adjust(p_value, method = "fdr")) %>% ungroup()
  } else {
    pub_df <- pub_df %>% mutate(p_fdr = NA_real_)
  }
  
  pub_df <- pub_df %>%
    mutate(
      Outcome_Clean = case_when(
        Outcome == "Theta.4.8Hz" ~ "Theta (4-8 Hz)",
        Outcome == "Alpha.8.12Hz" ~ "iAPF (8-12 Hz)",
        Outcome == "Exponent.30.48Hz" ~ "Exponent (30-48 Hz)",
        TRUE ~ Outcome
      ),
      Specparam = if_else(grepl("Exponent", Outcome), "Exponent", 
                          if_else(Periodic_param == "PW", "Power", "Center frequency")),
      Term_Clean = case_when(
        Term == "(Intercept)" ~ "Intercept",
        grepl("^Group_SSD", Term) & !grepl(":", Term) ~ "Group (SSD)",
        grepl("^Condition", Term) & !grepl(":", Term) ~ "Condition (EO)",
        grepl("age", Term, ignore.case = TRUE) ~ "Age",
        grepl("sex", Term, ignore.case = TRUE) ~ "Sex (Male)",
        grepl("CPZ", Term, ignore.case = TRUE) ~ "CPZ Equivalent",
        TRUE ~ str_replace_all(Term, "_", " ")
      ),
      Term_Clean = str_replace_all(Term_Clean, c("Group SSDSSD" = "Group (SSD)", "ConditionEO" = "Condition (EO)", "SessionV3" = "Session (V3)", ":" = " \u00D7 ")),
      `N (ID/Obs)` = paste0(N_unique_ID, " / ", N_observations),
      `Estimate (95% CI)` = sprintf("%.3f (%.3f, %.3f)", estimate, low_CI, high_CI),
      `p (uncorr.)` = if_else(p_value < 0.001, "< .001", sprintf("%.3f", p_value)),
      `sig_stars` = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**", p_value < 0.05 ~ "*", TRUE ~ ""),
      `p (FDR)` = if_else(p_fdr < 0.001, "< .001", sprintf("%.3f", p_fdr)),
      `fdr_stars` = case_when(p_fdr < 0.001 ~ "***", p_fdr < 0.01 ~ "**", p_fdr < 0.05 ~ "*", TRUE ~ "")
    )
  
  if (is_regional) {
    pub_df <- pub_df %>%
      mutate(Region_Clean = factor(str_replace_all(Region, "_", " "), 
                                   levels = c("Frontal", "Central", "Parietal", "Left Temporal", "Right Temporal", "Occipital"))) %>%
      arrange(Region_Clean, Outcome)
  }
  
  # Build final selection
  sel_cols <- c()
  if (is_regional) sel_cols <- c(sel_cols, "Region" = "Region_Clean")
  if (is_clinical) sel_cols <- c(sel_cols, "Clinical", "Predictor" = "Predictor_Set")
  sel_cols <- c(sel_cols, "Specparam", "Outcome" = "Outcome_Clean", "Term" = "Term_Clean", "N (ID/Obs)", "Estimate (95% CI)")
  sel_cols <- c(sel_cols, "SE" = "std_error", "t" = "t_val", "df", "p (uncorr.)", " " = "sig_stars")
  if (!is_regional) sel_cols <- c(sel_cols, "p (FDR)", "  " = "fdr_stars")
  
  pub_df %>% select(all_of(sel_cols)) %>%
    mutate(across(any_of(c("SE", "t", "df")), ~ sprintf("%.2f", as.numeric(.))))
}

## HELPER FUNCTION TO CREATE TABLES
format_publication_posthoc <- function(ph_df,
                                       main_df,
                                       is_clinical = FALSE,
                                       is_regional = FALSE) {
  
  if (is.null(ph_df) || nrow(ph_df) == 0) return(NULL)
  
  ph_df   <- standardize_lmm_cols(ph_df)
  main_df <- standardize_lmm_cols(main_df)
  
  # ----------------------------------------------------------
  # 1. KEEP ONLY CONFIRMATORY SPECS (STRICT)
  # ----------------------------------------------------------
  allowed_spec <- tibble::tibble(
    Outcome = c("Theta.4.8Hz", "Alpha.8.12Hz", "Exponent_30.48Hz"),
    Periodic_param = c("PW", "CF", "PW")
  )
  
  main_df <- main_df %>%
    semi_join(allowed_spec, by = c("Outcome", "Periodic_param"))
  
  ph_df <- ph_df %>%
    semi_join(allowed_spec, by = c("Outcome", "Periodic_param"))
  
  # ----------------------------------------------------------
  # 2. REGION FILTER
  # ----------------------------------------------------------
  if (is_regional) {
    main_df <- main_df %>% filter(Region != "Global")
  } else {
    main_df <- main_df %>% filter(Region == "Global")
  }
  
  # ----------------------------------------------------------
  # 3. CLINICAL FILTER
  # ----------------------------------------------------------
  if (is_clinical) {
    main_df <- main_df %>%
      filter(Predictor_Set %in% c("BACS_composite_z",
                                  "PANSS_neg",
                                  "PANSS_pos",
                                  "PANSS_gen"))
  }
  
  # ----------------------------------------------------------
  # 4. FIND SIGNIFICANT INTERACTIONS (STRICT MATCH)
  # ----------------------------------------------------------
  sig_models <- main_df %>%
    filter(grepl(":", Term) & p_value < 0.05) %>%
    select(Region, Outcome, Periodic_param, Predictor_Set) %>%
    distinct()
  
  if (nrow(sig_models) == 0) return(NULL)
  
  # ----------------------------------------------------------
  # 5. MATCH POSTHOC STRICTLY TO SIGNIFICANT MODELS
  # ----------------------------------------------------------
  join_cols <- c("Region", "Outcome", "Periodic_param", "Predictor_Set")
  
  ph_pub <- ph_df %>%
    semi_join(sig_models, by = join_cols)
  
  if (nrow(ph_pub) == 0) return(NULL)
  
  # ----------------------------------------------------------
  # 6. CLEAN OUTPUT NAMES
  # ----------------------------------------------------------
  ph_pub <- ph_pub %>%
    mutate(
      Outcome = case_when(
        Outcome == "Theta.4.8Hz" ~ "Theta (4-8 Hz)",
        Outcome == "Alpha.8.12Hz" ~ "iAPF (8-12 Hz)",
        Outcome == "Exponent_30.48Hz" ~ "Exponent (30-48 Hz)",
        TRUE ~ Outcome
      ),
      `p value` = ifelse(p_value < 0.001,
                         "< .001",
                         sprintf("%.3f", p_value)),
      Sig. = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  # ----------------------------------------------------------
  # 7. CONTINUOUS (SIMPLE SLOPES)
  # ----------------------------------------------------------
  if (is_clinical) {
    
    ph_pub <- ph_pub %>%
      mutate(
        Clinical = if_else(Predictor_Set == "BACS_composite_z",
                           "Cognition",
                           "Symptoms"),
        `N (ID/Obs)` = paste0(N_unique_ID, " / ", N_observations),
        `Slope (95% CI)` =
          sprintf("%.3f (%.3f, %.3f)",
                  Slope_Estimate,
                  low_CI,
                  high_CI)
      )
    
    sel <- c()
    if (is_regional) sel <- c(sel, "Region")
    
    sel <- c(sel,
             "Clinical",
             "Predictor" = "Trend_Variable",
             "Outcome",
             "Condition",
             "N (ID/Obs)",
             "Slope (95% CI)",
             "SE" = "Slope_SE",
             "t" = "t_val",
             "df",
             "p value",
             "Sig.")
    
    return(
      ph_pub %>%
        select(all_of(sel)) %>%
        mutate(across(c(SE, t, df),
                      ~ sprintf("%.2f", as.numeric(.))))
    )
  }
  
  # ----------------------------------------------------------
  # 8. CATEGORICAL
  # ----------------------------------------------------------
  if (!"SE" %in% names(ph_pub) && "std_error" %in% names(ph_pub)) {
    ph_pub <- ph_pub %>% rename(SE = std_error)
  }
  
  ph_pub <- ph_pub %>%
    mutate(
      `Estimate (95% CI)` =
        sprintf("%.3f (%.3f, %.3f)",
                estimate,
                low_CI,
                high_CI)
    )
  
  sel <- c()
  if (is_regional) sel <- c(sel, "Region")
  if ("Condition" %in% names(ph_pub)) sel <- c(sel, "Condition")
  
  sel <- c(sel,
           "Outcome",
           "Contrast" = "contrast",
           "Group (Count)",
           "Group (Mean)",
           "Estimate (95% CI)",
           "SE",
           "t" = "t_val",
           "df",
           "p value",
           "Sig.")
  
  return(
    ph_pub %>%
      select(all_of(sel)) %>%
      mutate(across(c(SE, t, df),
                    ~ sprintf("%.2f", as.numeric(.))))
  )
}



# ==============================================================================
# 3. GENERATE ALL PUBLICATION TABLES
# ==============================================================================
# Add Δ prefix to delta predictors in Genc tables (but NOT to covariates / intercept / condition)
add_delta_triangle <- function(df, posthoc = FALSE) {
  
  if (is.null(df) || nrow(df) == 0) return(df)
  
  if (!posthoc) {
    
    # ---- MAIN TABLES ----
    df <- df %>%
      dplyr::mutate(
        
        Term = dplyr::if_else(
          stringr::str_detect(Term,
                              "^(Intercept|Age|Sex \\(Male\\)|CPZ Equivalent|Condition \\(EO\\))$") |
            stringr::str_detect(Term, "^Condition"),
          Term,
          if_else(stringr::str_starts(Term, "\u0394"),
                  Term,
                  paste0("\u0394", Term))
        ),
        
        Outcome = dplyr::if_else(
          stringr::str_starts(Outcome, "\u0394"),
          Outcome,
          paste0("\u0394", Outcome)
        )
      )
    
  } else {
    
    # ---- POSTHOC TABLES ----
    df <- df %>%
      dplyr::mutate(
        
        Predictor = dplyr::if_else(
          stringr::str_starts(Predictor, "\u0394"),
          Predictor,
          paste0("\u0394", Predictor)
        ),
        
        Outcome = dplyr::if_else(
          stringr::str_starts(Outcome, "\u0394"),
          Outcome,
          paste0("\u0394", Outcome)
        )
      )
  }
  
  return(df)
}
# --- A. NINA (HC vs SSD) ---
Table_Nina_Confirmatory        <- format_publication_main(Results_Nina$raw_output$main_table, is_clinical = FALSE, is_regional = FALSE)
Table_Nina_Regional_Exploratory <- format_publication_main(Results_Nina$raw_output$main_table, is_clinical = FALSE, is_regional = TRUE)

Table_Nina_PostHoc_Global      <- format_publication_posthoc(Results_Nina$raw_output$posthoc_table, Results_Nina$raw_output$main_table, is_clinical = FALSE, is_regional = FALSE)
Table_Nina_PostHoc_Regional    <- format_publication_posthoc(Results_Nina$raw_output$posthoc_table, Results_Nina$raw_output$main_table, is_clinical = FALSE, is_regional = TRUE)



# --- B. DENIZ (Session V1 vs V3) ---
Table_Deniz_Confirmatory        <- format_publication_main(Results_Deniz$raw_output$main_table, is_clinical = FALSE, is_regional = FALSE)
Table_Deniz_Regional_Exploratory <- format_publication_main(Results_Deniz$raw_output$main_table, is_clinical = FALSE, is_regional = TRUE)

Table_Deniz_PostHoc_Global      <- format_publication_posthoc(Results_Deniz$raw_output$posthoc_table, Results_Deniz$raw_output$main_table, is_clinical = FALSE, is_regional = FALSE)
Table_Deniz_PostHoc_Regional    <- format_publication_posthoc(Results_Deniz$raw_output$posthoc_table, Results_Deniz$raw_output$main_table, is_clinical = FALSE, is_regional = TRUE)

# --- C. BERKHAN (Baseline Clinical) ---
Table_Berkhan_Confirmatory        <- format_publication_main(Results_Berkhan$main_table, is_clinical = TRUE, is_regional = FALSE)
Table_Berkhan_Regional_Exploratory <- format_publication_main(Results_Berkhan$main_table, is_clinical = TRUE, is_regional = TRUE)


Table_Berkhan_PostHoc_Global      <- format_publication_posthoc(Results_Berkhan$posthoc_table, Results_Berkhan$main_table, is_clinical = TRUE, is_regional = FALSE)
Table_Berkhan_PostHoc_Regional    <- format_publication_posthoc(Results_Berkhan$posthoc_table, Results_Berkhan$main_table, is_clinical = TRUE, is_regional = TRUE)

# --- D. GENC (Delta Clinical) ---
Table_Genc_Confirmatory        <- format_publication_main(Results_Genc$main_table, is_clinical = TRUE, is_regional = FALSE)%>%add_delta_triangle()
Table_Genc_Regional_Exploratory <- format_publication_main(Results_Genc$main_table, is_clinical = TRUE, is_regional = TRUE)%>%add_delta_triangle()


Table_Genc_PostHoc_Global      <- format_publication_posthoc(Results_Genc$posthoc_table, Results_Genc$main_table, is_clinical = TRUE, is_regional = FALSE)%>%add_delta_triangle(posthoc = TRUE)
Table_Genc_PostHoc_Regional    <- format_publication_posthoc(Results_Genc$posthoc_table, Results_Genc$main_table, is_clinical = TRUE, is_regional = TRUE)%>%add_delta_triangle(posthoc = TRUE)





# ------------------------------------------------------------
#  Define output folder
# ------------------------------------------------------------
output_dir <- "/Users/genchasanaj/Library/CloudStorage/GoogleDrive-gencxhasanaj@gmail.com/My Drive/LMU Klinikum/Projects/OscillationInMotion/Publication/Tables"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------------------------------------------
# 2. Put all tables into a named list
# ------------------------------------------------------------
all_tables <- list(
  Table_Nina_Confirmatory = Table_Nina_Confirmatory,
  Table_Nina_Regional_Exploratory = Table_Nina_Regional_Exploratory,
  Table_Nina_PostHoc_Global = Table_Nina_PostHoc_Global,
  Table_Nina_PostHoc_Regional = Table_Nina_PostHoc_Regional,
  
  Table_Deniz_Confirmatory = Table_Deniz_Confirmatory,
  Table_Deniz_Regional_Exploratory = Table_Deniz_Regional_Exploratory,
  Table_Deniz_PostHoc_Global = Table_Deniz_PostHoc_Global,
  Table_Deniz_PostHoc_Regional = Table_Deniz_PostHoc_Regional,
  
  Table_Berkhan_Confirmatory = Table_Berkhan_Confirmatory,
  Table_Berkhan_Regional_Exploratory = Table_Berkhan_Regional_Exploratory,
  Table_Berkhan_PostHoc_Global = Table_Berkhan_PostHoc_Global,
  Table_Berkhan_PostHoc_Regional = Table_Berkhan_PostHoc_Regional,
  
  Table_Genc_Confirmatory = Table_Genc_Confirmatory,
  Table_Genc_Regional_Exploratory = Table_Genc_Regional_Exploratory,
  Table_Genc_PostHoc_Global = Table_Genc_PostHoc_Global,
  Table_Genc_PostHoc_Regional = Table_Genc_PostHoc_Regional
)

# ------------------------------------------------------------
# 3. Save each non-null table
# ------------------------------------------------------------
iwalk(all_tables, function(tbl, name) {
  if (!is.null(tbl)) {
    write_csv(tbl, file.path(output_dir, paste0(name, ".csv")))
  }
})




