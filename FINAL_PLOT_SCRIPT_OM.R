#### ===========================================================================
#### ULTIMATE EEG VISUALIZATION PIPELINE (NINA, DENIZ, BERKHAN, GENC)
#### ===========================================================================

#### - STEP 1. Load Libraries & Setup ####
libraries <- c('pacman', 'dplyr', 'tidyr', 'ggplot2', 'stringr', 'emmeans', 
               'lme4', 'lmerTest', 'broom.mixed', 'purrr', 'ggdist', 'glue')

do.call(pacman::p_load, as.list(libraries))
options(warn = -1, scipen = 999)
set.seed(123)

# Create Output Directories
dir.create("../results/BioRender_Assets/Maps", recursive = TRUE, showWarnings = FALSE)
dir.create("../results/BioRender_Assets/Scatters", recursive = TRUE, showWarnings = FALSE)

#### - STEP 2. Load the RDS files ####
message("Loading RDS Files...")
Results_Nina    <- readRDS("../results/Nina_Outputs/Results_Nina_Master.rds")
Results_Deniz   <- readRDS("../results/Deniz_Outputs/Results_Deniz_Master.rds")
Results_Berkhan <- readRDS("../results/Berkhan_Outputs/Results_Berkhan_Master.rds")
Results_Genc    <- readRDS("../results/Genc_Outputs/Results_Genc_Master.rds")

#### - STEP 3. Prep Stats Tables ####
message("Preparing Statistical Dataframes...")

prep_stats <- function(df, is_clinical = FALSE) {
  raw_df <- if("raw_output" %in% names(df)) df$raw_output$main_table else df$main_table
  
  df_clean <- raw_df %>% 
    filter(
      (Outcome == "Exponent_30.48Hz" & Periodic_param == "PW") |
        (Outcome == "Alpha.8.12Hz"  & Periodic_param == "CF") |
        (Outcome == "Theta.4.8Hz"   & Periodic_param == "PW")
    ) %>%
    mutate(Analysis = if_else(Region == "Global", "Confirmatory", "Regional Exploratory"))
  
  if(is_clinical) {
    df_clean <- df_clean %>%
      filter(Term %in% c("BACS_composite_z", "PANSS_neg", "PANSS_pos", "PANSS_gen")) %>%
      mutate(Clinical = if_else(Term == "BACS_composite_z", "Cognition", "Symptoms")) %>%
      group_by(Region, Clinical, Term)
  } else {
    df_clean <- df_clean %>% group_by(Region, Term)
  }
  
  df_clean %>% 
    mutate(p_fdr = p.adjust(p_value, method = "fdr")) %>% 
    ungroup() %>%
    mutate(
      p_report = case_when(
        Analysis == "Confirmatory" ~ paste0("p = ", sprintf("%.3f", p_value), "; q = ", sprintf("%.3f", p_fdr)),
        TRUE ~ paste0("p = ", sprintf("%.3f", p_value))
      )
    )
}

nina_baseline_stats    <- prep_stats(Results_Nina, is_clinical = FALSE)
deniz_baseline_stats   <- prep_stats(Results_Deniz, is_clinical = FALSE)
berkhan_baseline_stats <- prep_stats(Results_Berkhan, is_clinical = TRUE)
genc_baseline_stats    <- prep_stats(Results_Genc, is_clinical = TRUE)

# Helper: Fetch Raw Model Data safely
get_model_data <- function(results_obj, region, param, outcome, predictor) {
  key <- paste(region, param, outcome, predictor, sep = "__")
  model_list <- if("models" %in% names(results_obj)) results_obj$models else results_obj$raw_output$models
  if (key %in% names(model_list)) return(model_list[[key]]@frame)
  key_bw <- paste(region, "BW", outcome, predictor, sep = "__")
  if (key_bw %in% names(model_list)) return(model_list[[key_bw]]@frame)
  return(NULL)
}

#### - STEP 4. Map Geometry & Plotting Function ####
create_arc <- function(start, end, radius = 10, n = 50) {
  theta <- seq(start, end, length.out = n) * pi / 180
  data.frame(x = radius * cos(theta), y = radius * sin(theta))
}
brain_geometry <- bind_rows(
  rbind(data.frame(x = 3, y = 4), data.frame(x = -3, y = 4), create_arc(135, 45)) %>% mutate(Region = "Frontal"),
  data.frame(x = c(-3, 3, 3, -3), y = c(4, 4, -0.5, -0.5), Region = "Central"),
  data.frame(x = c(-3, 3, 3, -3), y = c(-0.5, -0.5, -5, -5), Region = "Parietal"),
  rbind(data.frame(x = -3, y = -5), data.frame(x = 3, y = -5), create_arc(-45, -135)) %>% mutate(Region = "Occipital"),
  rbind(data.frame(x = -3, y = -5), data.frame(x = -3, y = 4), create_arc(135, 225)) %>% mutate(Region = "Left_Temporal"),
  rbind(data.frame(x = 3, y = 4), data.frame(x = 3, y = -5), create_arc(-45, 45)) %>% mutate(Region = "Right_Temporal")
)
label_coords <- data.frame(
  Region = c("Frontal", "Central", "Parietal", "Left_Temporal", "Right_Temporal", "Occipital"),
  Label = c("Frontal", "Central", "Parietal", "Left\nTemporal", "Right\nTemporal", "Occipital"),
  x = c(0, 0, 0, -6.5, 6.5, 0), y = c(7.5, 1.75, -2.75, -0.5, -0.5, -7.5)
)

plot_eeg_map <- function(stats_df, plot_title = "", legend_title = paste0("Estimate (B)"),  #expression("Estimate (" * italic(B) * ")")
                         low_color = "#2166AC", mid_color = "white", high_color = "#B2182B", 
                         show_labels = TRUE, label_color = "black", star_color = "black") {
  
  if(!"sig_star" %in% names(stats_df)) stats_df <- stats_df %>% mutate(sig_star = if_else(p_value < 0.05, "*", ""))
  map_data <- brain_geometry %>% left_join(stats_df, by = "Region")
  map_labels <- label_coords %>% left_join(stats_df, by = "Region")
  max_est <- max(abs(map_data$estimate), na.rm = TRUE)
  if(is.infinite(max_est) || max_est == 0) max_est <- 1 
  
  p_map <- ggplot() +
    geom_polygon(data = map_data, aes(x = x, y = y, group = Region, fill = estimate), color = "grey30", size = 1.5) +
    annotate("polygon", x = c(-1.5, 0, 1.5), y = c(9.9, 11.9, 9.9), fill = "white", color = "grey30", size = 1.5) +
    annotate("polygon", x = c(-9.8, -11, -11, -9.8), y = c(2, 1, -2, -3), fill = "white", color = "grey30", size = 1.5) +
    annotate("polygon", x = c(9.8, 11, 11, 9.8), y = c(2, 1, -2, -3), fill = "white", color = "grey30", size = 1.5) +
    geom_text(data = map_labels, aes(x = x, y = y - 2, label = sig_star), size = 20, fontface = "bold", color = star_color) +
    scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color, midpoint = 0, limits = c(-max_est, max_est), na.value = "grey90", name = legend_title, guide = guide_colorbar(barheight = unit(3, "in"), barwidth = unit(0.3, "in"), frame.colour = "black", ticks.colour = "black")) +
    coord_fixed(ratio = 1) + theme_void(base_size = 16) +
    theme(legend.position = "right", legend.title = element_text(face = "plain", size = 18, margin = margin(b = 10)),
          legend.text = element_text(size = 14), plot.title = element_text(face = "bold", size = 20, hjust = 0.5, margin = margin(b = 20)),
          plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_rect(fill = "transparent", color = NA))
  
  if(show_labels) p_map <- p_map + geom_text(data = map_labels, aes(x = x, y = y, label = Label), size = 5, fontface = "bold", color = label_color, lineheight = 0.9)
  return(p_map)
}

#### - STEP 5. Scatter & Regression Plotting Function ####
plot_eeg_scatter <- function(raw_data, 
                             plot_type = c("cat_main", "cat_int", "cont_main", "cont_int"),
                             x_var, 
                             y_var, 
                             condition_var = "Condition", 
                             stats_df = NULL, 
                             plot_title = "", 
                             x_label = NULL,
                             y_label = "",
                             color_palette = c("#3667E9", "#FC766D"),
                             show_legend = FALSE,
                             paired_id_var = NULL) {
  
  plot_type <- match.arg(plot_type)
  
  # 1. Clean data safely
  plot_data <- raw_data %>% filter(!is.na(.data[[y_var]]), !is.na(.data[[x_var]]))
  if (plot_type %in% c("cat_int", "cont_int")) plot_data <- plot_data %>% filter(!is.na(.data[[condition_var]]))
  
  # SAFEGUARD
  if(nrow(plot_data) == 0) {
    message("  -> Skipping plot: Missing data for variables.")
    return(NULL) 
  }
  
  # Calculate dynamic Y limits for annotations
  y_max <- max(plot_data[[y_var]], na.rm = TRUE)
  y_min <- min(plot_data[[y_var]], na.rm = TRUE)
  y_range <- y_max - y_min
  x_mid <- mean(range(as.numeric(plot_data[[x_var]]), na.rm=TRUE)) 
  x_left <- min(as.numeric(plot_data[[x_var]]), na.rm=TRUE) + diff(range(as.numeric(plot_data[[x_var]]), na.rm=TRUE)) * 0.05
  
  # Prepare categorical significance markers
  anno_df <- NULL
  if (!is.null(stats_df)) {
    anno_df <- stats_df %>%
      mutate(sig_star = case_when(
        p_value < 0.001 ~ "***", p_value < 0.01 ~ "**", p_value < 0.05 ~ "*", TRUE ~ ""
      )) %>% filter(sig_star != "")
  }
  
  # BUILD BASE PLOT
  p <- ggplot(plot_data) +
    scale_fill_manual(values = color_palette) +
    scale_color_manual(values = color_palette) +
    theme_classic(base_size = 16) +
    theme(
      strip.background = element_rect(fill = "white", color = "black", size = 1), 
      strip.text = element_text(face = "bold", size = 16, margin = margin(t=8, b=8)),
      axis.line = element_line(color = "black", size = 1), axis.ticks = element_line(color = "black", size = 1), 
      axis.title.y = element_text(face = "bold", size = 20, margin = margin(r = 15)), 
      axis.title.x = element_text(face = "bold", size = 20, margin = margin(t = 15)), 
      axis.text.x = element_text(face = "bold", size = 18, color = "black"),          
      axis.text.y = element_text(face = "bold", size = 18, color = "black"),          
      plot.title = element_text(face = "bold", size = 20, hjust = 0, margin = margin(b = 20)), 
      plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_rect(fill = "transparent", color = NA),
      legend.position = ifelse(show_legend, "right", "none"), legend.title = element_blank(), legend.text = element_text(face = "bold", size = 14)
    ) +
    labs(title = plot_title, y = y_label, x = ifelse(is.null(x_label), x_var, x_label))
  
  if (plot_type %in% c("cat_main", "cat_int")) {
    # --- CATEGORICAL: RAINCLOUD PLOTS ---
    p <- p + ggdist::stat_halfeye(aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]], color = .data[[x_var]]), adjust = 0.6, width = 0.6, .width = 0, justification = -0.3, point_colour = NA, alpha = 0.6)
    
    # Check if paired lines are requested and the ID variable exists
    if(!is.null(paired_id_var) && paired_id_var %in% names(plot_data)) {
      # Draw dotted paired lines first so they are behind points
      p <- p + geom_line(aes(x = .data[[x_var]], y = .data[[y_var]], group = .data[[paired_id_var]]), 
                         color = "grey50", linetype = "dotted", alpha = 0.6, size = 0.8)
      # Add points without jitter so they align perfectly with the lines
      p <- p + geom_point(aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[x_var]]), size = 2.5, alpha = 0.8)
    } else {
      # Standard jittered points for un-paired data
      p <- p + geom_point(aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[x_var]]), position = position_jitter(width = 0.1, height = 0, seed = 123), size = 2.5, alpha = 0.6)
    }
    
    # Add Boxplot
    p <- p + geom_boxplot(aes(x = .data[[x_var]], y = .data[[y_var]]), width = 0.25, outlier.shape = NA, alpha = 0, color = "black", fatten = 3)
    if (plot_type == "cat_int") p <- p + facet_wrap(as.formula(paste("~", condition_var)))
    
    # Categorical Annotations
    if (!is.null(anno_df) && nrow(anno_df) > 0) {
      p <- p + geom_segment(data = anno_df, aes(x = 1, xend = 2, y = y_max + y_range*0.06, yend = y_max + y_range*0.06), inherit.aes = FALSE, color = "black", size = 1) +
        geom_text(data = anno_df, aes(x = 1.5, y = y_max + y_range*0.12, label = sig_star), inherit.aes = FALSE, size = 12, fontface = "bold", color = "black")
    }
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.20)))
    
  } else {
    # --- CONTINUOUS: SCATTER & REGRESSION PLOTS ---
    if (plot_type == "cont_main") {
      p <- p + geom_point(aes(x = .data[[x_var]], y = .data[[y_var]]), color = "darkgray", alpha = 0.5, size = 3) +
        geom_smooth(aes(x = .data[[x_var]], y = .data[[y_var]]), method = "lm", color = color_palette[1], fill = color_palette[1], alpha = 0.2, size = 1.5)
      
      if (!is.null(stats_df)) {
        anno_text <- if("p_report" %in% names(stats_df)) stats_df$p_report[1] else paste0("p = ", sprintf("%.3f", stats_df$p_value[1]))
        p <- p + annotate("text", x = x_mid, y = Inf, label = anno_text, hjust = 0.5, vjust = 1.3, fontface = "plain", size = 6, color = "black")
      }
      
    } else if (plot_type == "cont_int") {
      p <- p + geom_point(aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[condition_var]]), alpha = 0.3, size = 3) +
        geom_smooth(aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[condition_var]], fill = .data[[condition_var]]), method = "lm", alpha = 0.2, size = 1.5)
      
      if (!is.null(stats_df) && "p_value" %in% names(stats_df)) {
        ec_val <- stats_df$p_value[stats_df[[condition_var]] %in% c("EC", "Eyes Closed")]
        eo_val <- stats_df$p_value[stats_df[[condition_var]] %in% c("EO", "Eyes Open")]
        ec_text <- paste0("Eyes Closed: p = ", if(length(ec_val) > 0 && !is.na(ec_val[1])) sprintf("%.3f", ec_val[1]) else "N/A")
        eo_text <- paste0("Eyes Open: p = ", if(length(eo_val) > 0 && !is.na(eo_val[1])) sprintf("%.3f", eo_val[1]) else "N/A")
        
        p <- p + annotate("text", x = -Inf, y = Inf, label = paste0(ec_text, "\n", eo_text), hjust = -0.05, vjust = 1.5, fontface = "plain", size = 6, color = "black", lineheight = 1.2)
      }
    }
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) 
  }
  return(p)
}


# ==============================================================================
# THE MASTER LOOP FUNCTION 
# ==============================================================================
clean_term_name <- function(x) {
  case_when(
    x == "BACS_composite_z" ~ "BACS Composite Score (z)",
    x == "PANSS_gen"         ~ "PANSS General",
    x == "PANSS_pos"         ~ "PANSS Positive",
    x == "PANSS_neg"         ~ "PANSS Negative",
    x == "Group_SSDSSD"      ~ "Group (SSD)",
    x == "SessionV3"         ~ "Session", 
    TRUE ~ str_replace_all(x, "_", " ")
  )
}



#### - STEP 3. The Master Loop Function ####
# THE MASTER LOOP FUNCTION (WITH SYMBOL PREFIXES & CUSTOM LABELS)
# ==============================================================================
clean_term_name <- function(x) {
  case_when(
    x == "BACS_composite_z" ~ "BACS Composite Score (z)",
    x == "PANSS_gen"         ~ "PANSS General",
    x == "PANSS_pos"         ~ "PANSS Positive",
    x == "PANSS_neg"         ~ "PANSS Negative",
    x == "Group_SSDSSD"      ~ "Group (SSD)",
    x == "SessionV3"         ~ "Session", 
    TRUE ~ str_replace_all(x, "_", " ")
  )
}

generate_all_plots <- function(name, results_obj, stats_table, is_clinical, predictors, main_term, 
                               cat_contrast = "HC - SSD", 
                               colors = c("#3667E9", "#FC766D"),
                               level_labels = NULL,         
                               title_prefix = "",           
                               x_prefix = "",               
                               y_prefix = "",               
                               paired_id_var = NULL,       
                               map_dims = c(8, 7), 
                               main_dims = c(8, 7), 
                               int_dims = c(8, 7)) {
  
  path_main <- file.path("../results/BioRender_Assets/Scatters", name, "Main")
  path_int  <- file.path("../results/BioRender_Assets/Scatters", name, "Int")
  path_maps <- file.path("../results/BioRender_Assets/Maps", name)
  dir.create(path_main, recursive = TRUE, showWarnings = FALSE)
  dir.create(path_int, recursive = TRUE, showWarnings = FALSE)
  dir.create(path_maps, recursive = TRUE, showWarnings = FALSE)
  
  regions <- c("Global", "Frontal", "Central", "Parietal", "Left_Temporal", "Right_Temporal", "Occipital")
  outcomes <- data.frame(
    outcome = c("Exponent_30.48Hz", "Alpha.8.12Hz", "Theta.4.8Hz"), param = c("PW", "CF", "PW"),
    title = c("Aperiodic Exponent (30-48 Hz)", "Alpha CF (8-12 Hz)", "Theta Power (4-8 Hz)"),
    ylabel = c("Exponent", "Frequency (Hz)", "Power")
  )
  
  ph_table <- if("raw_output" %in% names(results_obj)) results_obj$raw_output$posthoc_table else results_obj$posthoc_table
  
  for(pred in predictors) {
    term_filter <- if(is_clinical) pred else main_term
    pretty_pred_name <- clean_term_name(pred)
    
    for(i in 1:nrow(outcomes)) {
      o_name <- outcomes$outcome[i]; p_name <- outcomes$param[i]
      title_base <- paste0(title_prefix, outcomes$title[i])
      y_lab <- paste0(y_prefix, outcomes$ylabel[i])
      x_lab <- paste0(x_prefix, pretty_pred_name)
      
      # --- A. MAPS ---
      st_map <- stats_table %>% filter(Outcome == o_name, Periodic_param == p_name, Term == term_filter)
      if(nrow(st_map) > 0) {
        p_map <- plot_eeg_map(st_map, plot_title = paste(title_base, "\nEffect of", pretty_pred_name), star_color = "white", low_color = colors[1], high_color = colors[2])
        ggsave(file.path(path_maps, paste0("Map_", name, "_", o_name, "_", pred, ".png")), p_map, width = map_dims[1], height = map_dims[2], dpi = 300, bg = "transparent")
      }
      
      # --- B. SCATTERS ---
      for(reg in regions) {
        raw_df <- get_model_data(results_obj, reg, p_name, o_name, pred)
        if(is.null(raw_df) || nrow(raw_df) == 0) next 
        
        # 1. Main Scatters
        st_main <- stats_table %>% filter(Region == reg, Outcome == o_name, Term == term_filter)
        if (nrow(st_main) > 0) {
          type_main <- if(is_clinical) "cont_main" else "cat_main"
          raw_df_main <- raw_df 
          if(!is_clinical) {
            if (!is.null(level_labels) && pred %in% names(raw_df_main)) raw_df_main[[pred]] <- factor(as.character(raw_df_main[[pred]]), levels = names(level_labels), labels = level_labels)
            else if (pred %in% names(raw_df_main)) raw_df_main[[pred]] <- factor(raw_df_main[[pred]])
          }
          
          p_main <- plot_eeg_scatter(raw_df_main, type_main, x_var = pred, y_var = o_name, stats_df = st_main, 
                                     plot_title = paste(reg, title_base), x_label = x_lab, y_label = y_lab, 
                                     color_palette = if(is_clinical) c(colors[2], colors[2]) else colors,
                                     paired_id_var = paired_id_var) # Passed here
          
          if(!is.null(p_main)) ggsave(file.path(path_main, paste0(reg, "_", o_name, "_", pred, ".png")), p_main, width = main_dims[1], height = main_dims[2], dpi = 300, bg = "transparent")
        }
        
        # 2. Interaction Scatters
        if(!is.null(ph_table) && "Condition" %in% names(raw_df)) {
          raw_df_int <- raw_df %>% mutate(Condition = case_when(Condition %in% c("EC", "Eyes Closed") ~ "Eyes Closed", Condition %in% c("EO", "Eyes Open") ~ "Eyes Open", TRUE ~ as.character(Condition))) %>% filter(!is.na(Condition))
          if(nrow(raw_df_int) > 0) {
            st_ph <- if(!is_clinical) ph_table %>% filter(Region == reg, Outcome == o_name, contrast == cat_contrast) else ph_table %>% filter(Region == reg, Outcome == o_name, Trend_Variable == pred)
            if(nrow(st_ph) > 0) {
              st_ph_int <- st_ph %>% mutate(Condition = case_when(Condition %in% c("EC", "Eyes Closed") ~ "Eyes Closed", Condition %in% c("EO", "Eyes Open") ~ "Eyes Open", TRUE ~ as.character(Condition)))
              if("p.value" %in% names(st_ph_int)) st_ph_int <- st_ph_int %>% rename(p_value = p.value)
              
              if(!is_clinical && !is.null(level_labels) && pred %in% names(raw_df_int)) {
                raw_df_int[[pred]] <- factor(as.character(raw_df_int[[pred]]), levels = names(level_labels), labels = level_labels)
              }
              
              p_int <- plot_eeg_scatter(raw_df_int, if(is_clinical) "cont_int" else "cat_int", x_var = pred, y_var = o_name, condition_var = "Condition", 
                                        stats_df = st_ph_int, plot_title = paste(reg, title_base, "by Condition"), x_label = x_lab, y_label = y_lab, 
                                        show_legend = TRUE, color_palette = colors,
                                        paired_id_var = paired_id_var) # Passed here
              
              if(!is.null(p_int)) ggsave(file.path(path_int, paste0(reg, "_", o_name, "_", pred, ".png")), p_int, width = int_dims[1], height = int_dims[2], dpi = 300, bg = "transparent")
            }
          }
        }
      } 
    } 
  } 
  message(paste("Successfully generated all assets for:", name))
}

message("Generating All Scatter Plots...")

#Nina
generate_all_plots(
  name = "Nina", 
  results_obj = Results_Nina, 
  stats_table = nina_baseline_stats, 
  is_clinical = FALSE, 
  predictors = "Group_SSD", 
  main_term = "Group_SSDSSD", 
  cat_contrast = "HC - SSD", 
  colors = c("#3667E9", "#FC766D"),
  #level_labels = NULL, 
  map_dims = c(7, 7), main_dims = c(8, 7), int_dims = c(9, 7) 
)

# Deniz
generate_all_plots(
  name = "Deniz", 
  results_obj = Results_Deniz, 
  stats_table = deniz_baseline_stats, 
  is_clinical = FALSE, 
  predictors = "Session", 
  main_term = "SessionV3", 
  cat_contrast = "V1 - V3", 
  colors = c("#762A83", "#1B7837"),
  level_labels = c("V1" = "Pre", "V3" = "Post"), 
  paired_id_var = "ID", # <--- DRAWS DOTTED LINES BETWEEN PRE AND POST
  map_dims = c(8, 7), main_dims = c(8, 7), int_dims = c(9, 7) 
)

#Berkhan
generate_all_plots(
  name = "Berkhan", 
  results_obj = Results_Berkhan, 
  stats_table = berkhan_baseline_stats, 
  is_clinical = TRUE, 
  predictors = c("BACS_composite_z", "PANSS_neg", "PANSS_pos", "PANSS_gen"), 
  main_term = "", 
  cat_contrast = NULL, 
  colors = c("#FEA104", "#028180"),
  level_labels = NULL, 
  map_dims = c(8, 7), main_dims = c(8, 7), int_dims = c(9, 7) 
)


#Genc
generate_all_plots(
  name = "Genc", 
  results_obj = Results_Genc, 
  stats_table = genc_baseline_stats, 
  is_clinical = TRUE, 
  predictors = c("BACS_composite_z", "PANSS_neg", "PANSS_pos", "PANSS_gen"), 
  main_term = "", 
  cat_contrast = NULL, 
  colors = c("#FEA104", "#028180"),
  level_labels = NULL, 
  #add delta symbol as prefix
  y_prefix = "\u0394",
  x_prefix = "\u0394",
  map_dims = c(8, 7), main_dims = c(8, 7), int_dims = c(9, 7) 
)
# 

message("DONE! Check the BioRender_Assets/Scatters folder.")
