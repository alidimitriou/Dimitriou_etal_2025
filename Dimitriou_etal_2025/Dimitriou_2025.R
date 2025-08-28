# ---------------------------------------------------------------
## Habitat Use Models
# ---------------------------------------------------------------
# Load required libraries
packages <- c(
  "bayesplot", "brms", "cmdstanr", "dplyr", "ggdist", "ggplot2",
  "iNEXT", "kableExtra", "knitr", "lubridate", "patchwork",
  "posterior", "purrr", "readr", "sf", "spdep",
  "tidybayes", "tidyr", "tidyverse"
)
invisible(lapply(packages, library, character.only = TRUE))

# Load data
model_data <- read.csv("~/Desktop/Dimitriou_etal_2025/data/model_data.csv")

# Prepare output directory
output_dir <- "~/Desktop/Dimitriou_etal_2025/brms_output"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Define function to run and save BRMS model
options(brms.backend = "cmdstanr", brms.parallel_chains = 4)
run_and_save_brms_model <- function(response_var, species_name) {
  model <- brm(
    formula = as.formula(
      paste(response_var, "~",
            "mean_ndvi_scaled + canopy_cover_scaled +",
            "rec_month_100_scaled * prox_to_trail_decay +",
            "(1 | Deployment.Location.ID) + offset(log_effort)")
    ),
    data = model_data,
    family = negbinomial(),
    prior = c(
      prior(normal(0, 2), class = "b"),
      prior(normal(0, 2), class = "Intercept"),
      prior(cauchy(0, 2), class = "sd")
    ),
    chains = 4,
    iter = 10000,
    warmup = 4000,
    seed = 123,
    control = list(adapt_delta = 0.999, max_treedepth = 20),
    save_pars = save_pars(all = TRUE),
    init = 0
  )
  
  # Save model
  model_path <- file.path(output_dir, paste0(species_name, "_model.rds"))
  saveRDS(model, model_path)
  cat(paste0("Model saved: ", model_path, "\n"))
  
  # Save model summary
  summary_path <- file.path(output_dir, paste0(species_name, "_summary.txt"))
  sink(summary_path)
  print(summary(model))
  sink()
  cat(paste0("Summary saved: ", summary_path, "\n"))
  
  # Save effect size plot
  plot_and_save_effect_sizes(model, species_name)
  
  gc()
  return(model)
}

# Species response variable mapping
species_models <- list(
  "Hare" = "lepus_detections",
  "Black bear" = "ursus_detections",
  "Deer" = "odocoileus_detections",
  "Coyote" = "canis_detections",
  "Bobcat" = "lynx_detections",
  "Marten" = "martes_detections",
  "Marmot" = "marmota_detections",
  "Rare Carnivores" = "rare_detections"
)
species_models <- species_models[order(unlist(species_models))]

# Run and save each species model
for (species_name in names(species_models)) {
  response_var <- species_models[[species_name]]
  cat(paste0("\nRunning model for: ", species_name, "...\n"))
  run_and_save_brms_model(response_var, species_name)
}

# Initialize empty dataframe
summary_df <- data.frame()

# Loop through each saved model
for (species_name in names(species_models)) {
  model_path <- file.path(output_dir, paste0(species_name, "_model.rds"))
  
  if (!file.exists(model_path)) {
    warning(paste("Model file not found for", species_name))
    next
  }
  
  model <- readRDS(model_path)
  draws <- as_draws_df(model)
  
  # Get all population-level (fixed effect) terms
  params <- grep("^b_", names(draws), value = TRUE)
  ess_vals <- ess_bulk(draws)
  
  for (param in params) {
    vec <- draws[[param]]
    
    summary_df <- bind_rows(summary_df, tibble(
      Species = species_name,
      Parameter = param,
      Estimate = round(mean(vec), 2),
      Lower_95 = round(quantile(vec, 0.025), 2),
      Upper_95 = round(quantile(vec, 0.975), 2),
      Lower_80 = round(quantile(vec, 0.10), 2),
      Upper_80 = round(quantile(vec, 0.90), 2),
      Rhat = round(rhat(model)[param], 2)
    ))
  }
  
}

# Write to csv
write.csv(summary_df, file.path(output_dir, "posterior_estimates_by_species.csv"), row.names = FALSE)
print(summary_df)


# Save random effect standard deviation results
# Initialize empty dataframe
re_sd_df <- data.frame()

for (species_name in names(species_models)) {
  model_path <- file.path(output_dir, paste0(species_name, "_model.rds"))
  
  if (!file.exists(model_path)) {
    warning(paste("Model file not found for", species_name))
    next
  }
  
  model <- readRDS(model_path)
  sum_obj <- summary(model)
  
  if (!is.null(sum_obj$random$Deployment.Location.ID)) {
    re_row <- sum_obj$random$Deployment.Location.ID["sd(Intercept)", ]
    
    re_sd_df <- bind_rows(re_sd_df, tibble(
      Species = species_name,
      SD_Estimate = round(re_row["Estimate"], 2),
      SD_Lower_95 = round(re_row["l-95% CI"], 2),
      SD_Upper_95 = round(re_row["u-95% CI"], 2),
      SD_Rhat = round(re_row["Rhat"], 2)
    ))
  } else {
    warning(paste("Random effect not found in summary for", species_name))
  }
}

# Save and print
write.csv(re_sd_df, file.path(output_dir, "random_intercept_sds_by_species.csv"), row.names = FALSE)
print(re_sd_df)

# ---------------------------------------------------------------
## Effect size plots
# ---------------------------------------------------------------
# Variable labels
var_labels <- c(
  "b_mean_ndvi_scaled" = "NDVI",
  "b_prox_to_trail_decay" = "Proximity to Trail",
  "b_elevation_scaled" = "Elevation",
  "b_canopy_cover_scaled" = "Canopy Cover",
  "b_rec_month_100_scaled" = "Park-Month Recreation",
  "b_rec_month_100_scaled:prox_to_trail_decay" = "Recreation × Trail",
  "b_dist_to_trail_decay" = "Proximity to Trail", 
  "b_rec_month_100_scaled:dist_to_trail_decay" = "Recreation × Trail"
)

# Prepare data for plotting
plot_df <- summary_df %>%
  filter(Parameter %in% names(var_labels)) %>%
  rename(
    median = Estimate,
    lower95 = Lower_95,
    upper95 = Upper_95,
    lower80 = Lower_80,
    upper80 = Upper_80
  ) %>%
  # Filter: keep Recreation always, others only if 80% CI does NOT cross 0
  filter(
    Parameter == "b_rec_month_100_scaled" |
      (lower80 > 0 | upper80 < 0)
  ) %>%
  mutate(
    Variable = factor(var_labels[Parameter], levels = rev(unique(var_labels))),
    is_recreation = ifelse(Parameter == "b_rec_month_100_scaled", "Recreation", "Other")
  )

# Plot
combined_effects_plot <- ggplot(plot_df, aes(y = Variable)) +
  geom_segment(aes(x = lower80, xend = upper80, yend = Variable, color = is_recreation), linewidth = 3) +
  geom_segment(aes(x = lower95, xend = upper95, yend = Variable, color = is_recreation), linewidth = 1.2) +
  geom_point(aes(x = median, color = is_recreation), size = 4.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  facet_wrap(~Species, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("Recreation" = "#2E8B57", "Other" = "#084594")) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.5, "lines"),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    legend.position = "none"
  ) +
  labs(x = "Posterior Estimate", y = NULL)

# Save
ggsave(filename = file.path(output_dir, "Recreation_and_Significant_Other_Effects.png"),
       plot = combined_effects_plot, width = 12, height = 14, dpi = 300)


# ---------------------------------------------------------------   
# Autocorrelation tests
# ---------------------------------------------------------------
# Define species and detection columns
species_name_map <- list(
  "Hare" = "lepus_detections",
  "Black bear" = "ursus_detections",
  "Deer" = "odocoileus_detections",
  "Coyote" = "canis_detections",
  "Bobcat" = "lynx_detections",
  "Marten" = "martes_detections",
  "Marmot" = "marmota_detections",
  "Rare Carnivores" = "rare_detections")

parks <- c("Garibaldi", "Joffre Lakes")

moran_results <- data.frame(
  Species = character(),
  Park = character(),
  Month = character(),  
  Moran_I = numeric(),
  Expectation = numeric(),
  Variance = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through species
for (species_name in names(species_name_map)) {
  
  response_var <- species_name_map[[species_name]]
  cat("\nUsing detection column:", response_var, "for species:", species_name, "\n")
  
  model_path <- paste0("~/Desktop/Dimitriou_etal_2025/brms_output/", species_name, "_model.rds")
  if (!file.exists(model_path)) {
    cat("Model file not found for", species_name, "- Skipping.\n")
    next
  }
  
  # Load model and re join start_date and spatial info using row_index
  model <- readRDS(model_path)
  model$data$row_index <- 1:nrow(model$data)
  model_data$row_index <- 1:nrow(model_data)
  
  species_data <- left_join(
    model$data,
    model_data %>% select(row_index, start_date, latitude, longitude, park),
    by = "row_index"
  )
  
  if (!("park" %in% colnames(species_data))) {
    stop(paste("Missing 'park' column in model data for species:", species_name))
  }
  
  # Compute residuals and assign months
  predicted_values <- fitted(model, summary = TRUE)[, "Estimate"]
  
  species_data <- species_data %>%
    mutate(
      Residuals = !!sym(response_var) - predicted_values,
      Month = floor_date(as.Date(start_date), "month")
    )
  
  model_data_summary <- species_data %>%
    select(Month, park, Deployment.Location.ID, latitude, longitude, Residuals)
  
  # Loop through parks
  for (park_name in parks) {
    
    available_months <- unique(model_data_summary$Month[model_data_summary$park == park_name])
    if (length(available_months) == 0) {
      cat("No available months found for", species_name, "in", park_name, "- Skipping.\n")
      next
    }
    
    cat("Available months for", species_name, "in", park_name, ":", available_months, "\n")
    
    for (selected_month in available_months) {
      
      cat("\nRunning Moran's I for", species_name, "in", park_name, "for Month:", selected_month, "...\n")
      
      spatial_residuals <- model_data_summary %>%
        filter(park == park_name, Month == selected_month) %>%
        select(Deployment.Location.ID, latitude, longitude, Residuals) %>%
        filter(!is.na(Residuals))
      
      if (nrow(spatial_residuals) == 0) {
        cat("No residuals found for", species_name, "in", park_name, "Month:", selected_month, "- Skipping.\n")
        next
      }
      
      spatial_residuals_sf <- st_as_sf(spatial_residuals, coords = c("longitude", "latitude"), crs = 4326)
      num_stations <- nrow(spatial_residuals_sf)
      
      if (num_stations < 4) {
        cat("Not enough spatial locations for", species_name, "in", park_name, "Month:", selected_month, "- Skipping.\n")
        next
      }
      
      k_neighbors <- min(4, num_stations - 1)
      cat("Using", k_neighbors, "nearest neighbors\n")
      
      coords <- st_coordinates(spatial_residuals_sf)
      nb <- knearneigh(coords, k = k_neighbors) %>% knn2nb()
      listw <- nb2listw(nb, style = "W")
      
      moran_test <- moran.test(spatial_residuals$Residuals, listw)
      
      moran_results <- rbind(moran_results, data.frame(
        Species = species_name,
        Park = park_name,
        Month = as.character(selected_month),
        Moran_I = moran_test$estimate[1],
        Expectation = moran_test$estimate[2],
        Variance = moran_test$estimate[3],
        P_Value = moran_test$p.value
      ))
      
      print(moran_test)
    }
  }
}

# Save results
write.csv(moran_results, file.path(output_dir, "moran_results.csv"), row.names = FALSE)


cat("\nMoran's I monthly analysis complete. Results saved to:", output_dir, "\n")
print(moran_results)

# Number of significant Moran's I results (p < 0.05)
sum(moran_results$P_Value < 0.05, na.rm = TRUE)

# Species to test
species_list <- c("Black bear", "Rare Carnivores", "Marmot", "Marten", "Bobcat", "Hare", "Coyote", "Deer")

# Function to compute Durbin-Watson
compute_dw <- function(residuals) {
  sum_diff_sq <- sum(diff(residuals)^2)
  sum_sq <- sum(residuals^2)
  dw_stat <- sum_diff_sq / sum_sq
  return(dw_stat)
}

# Store results
dw_results <- data.frame(
  Species = character(),
  DW = numeric(),
  stringsAsFactors = FALSE
)

# Loop through species
for (species in species_list) {
  model_path <- file.path(output_dir, paste0(species, "_model.rds"))
  
  if (!file.exists(model_path)) {
    warning(paste("Model not found for", species))
    next
  }
  
  cat("Running manual Durbin-Watson for:", species, "\n")
  
  model <- readRDS(model_path)
  res <- tryCatch({
    residuals(model)[, 1]
  }, error = function(e) {
    warning(paste("Could not extract residuals for", species))
    return(NULL)
  })
  
  if (is.null(res)) next
  
  dw_stat <- compute_dw(res)
  
  dw_results <- rbind(dw_results, data.frame(
    Species = species,
    DW = round(dw_stat, 3)
  ))
}

# Save and print results
dir.create(file.path(output_dir, "diagnostics"), recursive = TRUE, showWarnings = FALSE)
write.csv(dw_results, file.path(output_dir, "diagnostics", "durbin_watson_results.csv"), row.names = FALSE)

print(dw_results)

# ---------------------------------------------------------------
# Model Diagnostics
# ---------------------------------------------------------------
# Create output folders
trace_dir <- "~/Desktop/Dimitriou_etal_2025/brms_output/diagnostics/trace_plots"
pp_dir <- "~/Desktop/Dimitriou_etal_2025/brms_output/diagnostics/pp_checks"
dir.create(trace_dir, showWarnings = FALSE)
dir.create(pp_dir, showWarnings = FALSE)

# Loop through models
for (f in model_files) {
  model <- readRDS(f)
  species_name <- gsub("_model.rds", "", basename(f))
  cat("\nSaving diagnostics for:", species_name, "\n")
  
# Save posterior predictive check (mean) 
  pp_plot <- pp_check(model, type = "stat", stat = "mean", ndraws = 1000) +
    ggtitle(paste("Posterior Predictive Check:", species_name)) +
    theme_bw(base_size = 12)
  
  ggsave(
    filename = file.path(pp_dir, paste0(species_name, "_ppcheck.png")),
    plot = pp_plot,
    width = 6,
    height = 4,
    dpi = 300
  )
  

# Save Trace Plot for Fixed Effects
  
  posterior_array <- as.array(model)
  b_params <- grep("^b_", dimnames(posterior_array)[[3]], value = TRUE)
  
  if (length(b_params) > 0) {
    trace_plot <- mcmc_trace(posterior_array, pars = b_params) +
      ggtitle(paste("Trace Plots:", species_name)) +
      theme_bw(base_size = 12)
    
    ggsave(
      filename = file.path(trace_dir, paste0(species_name, "_trace_plot.png")),
      plot = trace_plot,
      width = 10,
      height = 6,
      dpi = 300
    )
  } else {
    cat("No fixed effects found to plot.\n")
  }
}

# Prepare empty lists
pp_plots <- list()
trace_plots <- list()

# Collect plots
for (f in model_files) {
  model <- readRDS(f)
  species_name <- gsub("_model.rds", "", basename(f))
  cat("Preparing plots for:", species_name, "\n")
  
# Posterior predictive plot
  pp_plot <- pp_check(model, type = "stat", stat = "mean", ndraws = 100) +
    ggtitle(species_name) +
    theme_bw(base_size = 10)
  pp_plots[[species_name]] <- pp_plot
  
# Trace plot 
  posterior_array <- as.array(model)
  b_params <- grep("^b_", dimnames(posterior_array)[[3]], value = TRUE)
  if (length(b_params) > 0) {
    trace_plot <- mcmc_trace(posterior_array, pars = b_params) +
      ggtitle(species_name) +
      theme_bw(base_size = 10)
    trace_plots[[species_name]] <- trace_plot
  }
}

# Combine plots into one multi-panel layout
pp_combined <- wrap_plots(pp_plots, ncol = 3)
trace_combined <- wrap_plots(trace_plots, ncol = 1)

# Save combined plots
ggsave("~/Desktop/Dimitriou_etal_2025/brms_output/diagnostics/Combined_PP_Checks.pdf", pp_combined, width = 12, height = 9)
ggsave("~/Desktop/Dimitriou_etal_2025/brms_output/diagnostics/Combined_Trace_Plots.pdf", trace_combined, width = 10, height = 15)

# Extract Bayesian p-values 

# Initialize results list
bayes_pvals <- data.frame(
  Species = character(),
  Bayes_p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through models
for (f in model_files) {
  model <- readRDS(f)
  species_name <- gsub("_model.rds", "", basename(f))
  
  y_var <- as.character(model$formula$formula[[2]])
  y_obs <- model$data[[y_var]]
  
  # Generate posterior predictive draws
  y_rep <- posterior_predict(model, ndraws = 1000)
  
  # Test statistic: mean
  T_obs <- mean(y_obs)
  T_rep <- apply(y_rep, 1, mean)
  
  # Bayesian p-value
  p_bayes <- mean(T_rep >= T_obs)
  
  # Store in results
  bayes_pvals <- rbind(
    bayes_pvals,
    data.frame(Species = species_name, Bayes_p_value = round(p_bayes, 3))
  )
}

# View table
print(bayes_pvals)


# ---------------------------------------------------------------
## Community analysis 
# ---------------------------------------------------------------
   
# Define species to exclude
exclude_species <- c(
  "Unknown.species", "Homo.sapiens", "Bird.spp.", "Dendragapus.fuliginosus", "Tamias.striatus",
  "Mus.musculus", "Canis.familiaris", "Lagopus.leucura", "Felis.catus", "Neotamias.townsendii",
  "Phenacomys.intermedius", "Tamiasciurus.hudsonicus", "Glaucomys.sabrinus", "Myodes.gapperi",
  "Ochotona.princeps", "Peromyscus.maniculatus", "Tamiasciurus.douglasii", "Neotamias.minimus",
  "Neotamias.amoenus", "Sciurus.carolinensis", "Bonasa.umbellus", "Dendragapus.obscurus",
  "Falcipennis.canadensis", "Mustela.erminea"
)

# Generate and save accumulation curves
plot_species_accumulation_q <- function(data_path, park_name, q = 0, y_label = "Species Richness", y_limit = c(0, 17)) {
  park_obs <- read.csv(data_path) %>%
    mutate(Date = as.Date(Date, format = "%Y-%m-%d"))
  
  periods <- list(
    "2020" = interval(dmy("06-08-2020"), dmy("03-11-2020")),
    "2021" = interval(dmy("22-06-2021"), dmy("19-09-2021")),
    "2022" = interval(dmy("10-05-2022"), dmy("07-08-2022"))
  )
  
  year_colours <- c("2020" = "#E63946", "2021" = "#457B9D", "2022" = "#2A9D8F")
  plots <- list()
  
  for (year in names(periods)) {
    park_period <- park_obs %>%
      filter(Date %within% periods[[year]]) %>%
      group_by(Date) %>%
      summarise(
        Days = 1,
        across(where(is.numeric) & !any_of(c("Effort", "Deployment.Location.ID", exclude_species)), ~ sum(. > 0, na.rm = TRUE))
      )
    
    species_cols <- setdiff(names(park_period), c("Date", "Days"))
    if (length(species_cols) == 0) next
    
    inc_dat <- park_period %>% mutate(across(all_of(species_cols), ~ as.integer(. > 0)))
    
    total_days <- sum(park_period$Days)
    species_counts <- colSums(inc_dat[, species_cols, drop = FALSE])
    
    project_level <- list(project_level = c(total_days, species_counts))
    names(project_level$project_level)[1] <- "SamplingUnits"
    
    if (total_days > 0 && sum(species_counts) > 0) {
      out <- iNEXT(project_level, q = q, datatype = "incidence_freq", knots = 40, se = TRUE, conf = 0.95, nboot = 10)
      plot_data <- fortify(out, type = 1) %>%
        mutate(LineType = ifelse(Method == "Extrapolation", "dashed", "solid"))
      
      plots[[year]] <- ggplot(plot_data, aes(x = x, y = y)) +
        geom_line(aes(linetype = LineType), colour = year_colours[year], linewidth = 1.2) +
        geom_ribbon(aes(ymin = y.lwr, ymax = y.upr), fill = year_colours[year], alpha = 0.2) +
        scale_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed")) +
        geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
        geom_vline(xintercept = 0, colour = "black", linewidth = 0.5) +
        theme_minimal(base_size = 16) +
        theme(
          panel.grid = element_blank(),
          legend.position = "none",
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          plot.title = element_text(size = 18, face = "bold")
        ) +
        labs(
          title = paste0(park_name, " - ", year),
          x = "Sampling Effort (Days)",
          y = y_label
        ) +
        coord_cartesian(ylim = y_limit)
    }
  }
  
  return(plots)
}

# Generate accumulation curves from processed data
joffre_richness <- plot_species_accumulation_q(
  "~/Desktop/Dimitriou_etal_2025/data/JOFF_30min_Independent_daily_observations.csv",
  "Joffre Lakes", q = 0, y_label = "Species Richness", y_limit = c(0, 15)
)

garibaldi_richness <- plot_species_accumulation_q(
  "~/Desktop/Dimitriou_etal_2025/data/GARI_30min_Independent_daily_observations.csv",
  "Garibaldi", q = 0, y_label = "Species Richness", y_limit = c(0, 15)
)

joffre_shannon <- plot_species_accumulation_q(
  "~/Desktop/Dimitriou_etal_2025/data/JOFF_30min_Independent_daily_observations.csv",
  "Joffre Lakes", q = 1, y_label = "Shannon Diversity", y_limit = c(0, 8)
)

garibaldi_shannon <- plot_species_accumulation_q(
  "~/Desktop/Dimitriou_etal_2025/data/GARI_30min_Independent_daily_observations.csv",
  "Garibaldi", q = 1, y_label = "Shannon Diversity", y_limit = c(0, 8)
)


accumulation_grid <- (
  joffre_richness[["2020"]] | joffre_richness[["2021"]] | joffre_richness[["2022"]]
) /
  (
    garibaldi_richness[["2020"]] | garibaldi_richness[["2021"]] | garibaldi_richness[["2022"]]
  ) /
  (
    joffre_shannon[["2020"]] | joffre_shannon[["2021"]] | joffre_shannon[["2022"]]
  ) /
  (
    garibaldi_shannon[["2020"]] | garibaldi_shannon[["2021"]] | garibaldi_shannon[["2022"]]
  )
# Save plots
ggsave(
  filename = file.path("~/Desktop/Dimitriou_etal_2025/diversity_output", "species_accumulation_4x3_grid.png"),
  plot = accumulation_grid,
  width = 18,
  height = 12,
  dpi = 300
)

  
process_park_data_periods <- function(data_path, park_name) {
  park_obs <- read.csv(data_path) %>%
    mutate(Date = as.Date(Date))
  
  periods <- list(
    "2020" = interval(dmy("06-08-2020"), dmy("03-11-2020")),
    "2021" = interval(dmy("22-06-2021"), dmy("19-09-2021")),
    "2022" = interval(dmy("10-05-2022"), dmy("07-08-2022"))
  )
  
  richness_estimates <- shannon_estimates <- data.frame()
  
  for (year in names(periods)) {
    park_period <- park_obs %>% filter(Date %within% periods[[year]])
    if (nrow(park_period) == 0) next
    
    park_period <- park_period %>%
      select(-any_of(exclude_species)) %>%
      group_by(Date) %>%
      summarise(
        TotalEffort = sum(Effort, na.rm = TRUE),
        across(where(is.numeric) & !any_of("Effort"), ~ sum(. > 0, na.rm = TRUE))
      )
    
    species_cols <- setdiff(names(park_period), c("Date", "TotalEffort"))
    if (length(species_cols) == 0) next
    
    inc_dat <- park_period %>% mutate(across(all_of(species_cols), ~ as.integer(. > 0)))
    
    project_level <- list(project_level = c(sum(park_period$TotalEffort), colSums(inc_dat[, species_cols])))
    names(project_level$project_level)[1] <- "SamplingUnits"
    
    if (sum(project_level$project_level[-1]) > 0) {
      out <- iNEXT(project_level, q = c(0, 1), datatype = "incidence_freq", knots = 40, se = TRUE, conf = 0.95, nboot = 50)
      asy_estimates <- out$AsyEst %>%
        filter(Diversity %in% c("Species richness", "Shannon diversity")) %>%
        mutate(Year = year, Park = park_name) %>%
        select(Year, Park, Diversity, Estimator, LCL, UCL)
      
      richness_estimates <- bind_rows(richness_estimates, filter(asy_estimates, Diversity == "Species richness"))
      shannon_estimates <- bind_rows(shannon_estimates, filter(asy_estimates, Diversity == "Shannon diversity"))
    }
  }
  
  list(richness = richness_estimates, shannon = shannon_estimates)
}

# Run Diversity Processing for Each Park

joffre <- process_park_data_periods("~/Desktop/Dimitriou_etal_2025/data/JOFF_30min_Independent_daily_observations.csv", "Joffre Lakes")
garibaldi <- process_park_data_periods("~/Desktop/Dimitriou_etal_2025/data/GARI_30min_Independent_daily_observations.csv", "Garibaldi")

richness_estimates <- bind_rows(joffre$richness, garibaldi$richness)
shannon_estimates <- bind_rows(joffre$shannon, garibaldi$shannon)


# Calculate Evenness
custom_colors <- c("Joffre Lakes" = "#E63946", "Garibaldi" = "#457B9D")

evenness_estimates <- richness_estimates %>%
  rename(Richness = Estimator, LCL_R = LCL, UCL_R = UCL) %>%
  left_join(
    shannon_estimates %>%
      rename(Shannon_Hill = Estimator, LCL_S = LCL, UCL_S = UCL),
    by = c("Year", "Park")
  ) %>%
  mutate(
    Shannon = log(Shannon_Hill),
    Evenness = Shannon / log(Richness),
    SE_Shannon = (log(UCL_S) - log(LCL_S)) / (2 * 1.96),
    SE_Richness = (UCL_R - LCL_R) / (2 * 1.96),
    SE_Evenness = sqrt((1 / log(Richness))^2 * SE_Shannon^2 + (Shannon / (Richness * log(Richness)^2))^2 * SE_Richness^2),
    LCL_Evenness = Evenness - 1.96 * SE_Evenness,
    UCL_Evenness = Evenness + 1.96 * SE_Evenness
  )


# Plotting Functions and Outputs
plot_metrics <- function(data, y_label, show_legend = FALSE) {
  ggplot(data, aes(x = factor(Year), y = Estimator, color = Park)) +
    geom_point(size = 4, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.2, position = position_dodge(0.5)) +
    scale_color_manual(values = custom_colors) +
    theme_minimal(base_size = 20) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.position = if (show_legend) "right" else "none",
      plot.title = element_blank()
    ) +
    labs(x = "Year", y = y_label, color = "Park")
}

# Apply to each metric
richness_plot <- plot_metrics(richness_estimates, "Species Richness", show_legend = FALSE)
shannon_plot  <- plot_metrics(shannon_estimates, "Shannon Diversity", show_legend = TRUE)

# Evenness plot
evenness_plot <- ggplot(evenness_estimates, aes(x = factor(Year), y = Evenness, color = Park)) +
  geom_point(size = 4, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = LCL_Evenness, ymax = UCL_Evenness), width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.position = "none",
    plot.title = element_blank()
  ) +
  labs(x = "Year", y = "Pielou's Evenness", color = "Park")

# Combine plots
diversity_grid <- richness_plot / shannon_plot / evenness_plot

# Save plot
ggsave(filename = file.path(output_dir, "diversity_metrics_grid.png"),
       plot = diversity_grid, width = 10, height = 12, dpi = 300)

# Save estimates
dir.create("~/Desktop/Dimitriou_etal_2025/diversity_output", recursive = TRUE, showWarnings = FALSE)
saveRDS(richness_estimates, "~/Desktop/Dimitriou_etal_2025/diversity_output/richness_estimates.rds")
saveRDS(shannon_estimates, "~/Desktop/Dimitriou_etal_2025/diversity_output/shannon_estimates.rds")
saveRDS(evenness_estimates, "~/Desktop/Dimitriou_etal_2025/diversity_output/evenness_estimates.rds")
 

# ---------------------------------------------------------------
# Bootstrap hypothesis test
# ---------------------------------------------------------------

# Bootstrap simulation function
simulate_bootstrap <- function(mean, se, n = 10000) {
  rnorm(n, mean = mean, sd = se)
}

# Function to run bootstrapped comparison
compare_groups_bootstrap <- function(data, group_var, metric, level1, level2, n_boot) {
  sim1_list <- data %>% filter(!!sym(group_var) == level1) %>% pull(Simulated)
  sim2_list <- data %>% filter(!!sym(group_var) == level2) %>% pull(Simulated)

  sim1 <- sim1_list[[1]]
  sim2 <- sim2_list[[1]]
  
  est1 <- data %>% filter(!!sym(group_var) == level1) %>% pull(.data[[metric]])
  se1  <- data %>% filter(!!sym(group_var) == level1) %>% pull(.data[[paste0("SE_", metric)]])
  est2 <- data %>% filter(!!sym(group_var) == level2) %>% pull(.data[[metric]])
  se2  <- data %>% filter(!!sym(group_var) == level2) %>% pull(.data[[paste0("SE_", metric)]])
  
  diff_sim <- sim2 - sim1
  
  # Continuity-corrected two-sided p-value
  prop_ge_0 <- (sum(diff_sim >= 0) + 1) / (n_boot + 1)
  prop_le_0 <- (sum(diff_sim <= 0) + 1) / (n_boot + 1)
  boot_p_value <- 2 * min(prop_ge_0, prop_le_0)
  
  data.frame(
    Metric = metric,
    Comparison = paste(level1, "vs", level2),
    Mean_Diff = mean(diff_sim),
    CI_Lower = quantile(diff_sim, 0.025),
    CI_Upper = quantile(diff_sim, 0.975),
    Boot_P_Value = boot_p_value,
    Estimate_1 = est1,
    SE_1 = se1,
    Estimate_2 = est2,
    SE_2 = se2,
    N_Boot = n_boot
  )
}

# Run Bootstrap Comparisons
n_boot <- 10000
metrics <- c("Richness", "Shannon", "Evenness")
all_results <- list()

for (metric in metrics) {
  se_col <- paste0("SE_", metric)

  sim_estimates <- evenness_estimates %>%
    mutate(Simulated = map2(.data[[metric]], .data[[se_col]], ~ simulate_bootstrap(.x, .y, n_boot)))
  
  # Between-Year Comparisons (Within Each Park)
  for (park in unique(sim_estimates$Park)) {
    park_data <- sim_estimates %>% filter(Park == park)
    years <- unique(park_data$Year)
    
    for (i in 1:(length(years) - 1)) {
      for (j in (i + 1):length(years)) {
        res <- compare_groups_bootstrap(park_data, "Year", metric, years[i], years[j], n_boot)
        res$Park <- park
        res$Comparison_Type <- "Between-Year"
        all_results <- append(all_results, list(res))
      }
    }
  }
  
  # Between-Park Comparisons (Within Each Year)
  for (year in unique(sim_estimates$Year)) {
    year_data <- sim_estimates %>% filter(Year == year)
    parks <- unique(year_data$Park)
    
    if (length(parks) == 2) {
      res <- compare_groups_bootstrap(year_data, "Park", metric, parks[1], parks[2], n_boot)
      res$Year <- year
      res$Comparison_Type <- "Between-Park"
      all_results <- append(all_results, list(res))
    }
  }
}

# Combine and Correct Results (Benjamini-Hochberg)
results_df <- bind_rows(all_results)
results_df$Boot_BH_P_Value <- p.adjust(results_df$Boot_P_Value, method = "BH")

# Round numeric values
results_df <- results_df %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Save Results
write.csv(results_df, "~/Desktop/Dimitriou_etal_2025/diversity_output/diversity_bootstrap_BH_results_.csv", row.names = FALSE)
