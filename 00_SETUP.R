### 00_SETUP.R ###

#### 1. Analysis metadata ####

analysis_date <- "2026.07.17"
message("Running AR-attractor-producer analysis: ", analysis_date)


#### 2. Packages ####

required_packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "forcats", "lubridate", "tibble", "purrr",
  "glmmTMB", "lme4", "lmerTest", "emmeans", "broom.mixed", "performance", "DHARMa",
  "splines", "mgcv", "vegan", "ggplot2", "patchwork", "ggeffects", "here"
)

missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages)) stop("Install missing packages: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages(invisible(lapply(required_packages, library, character.only = TRUE)))


#### 3. Project paths ####

project_dir <- here::here()
data_raw_dir <- here::here("data", "raw")
data_interim_dir <- here::here("data", "interim")
data_processed_dir <- here::here("data", "processed")
outputs_root_dir <- here::here("outputs")
outputs_dir <- file.path(outputs_root_dir, paste0("Analysis_", analysis_date))
figures_dir <- file.path(outputs_dir, "figures")
tables_dir <- file.path(outputs_dir, "tables")
models_dir <- file.path(outputs_dir, "model_objects")
diagnostics_dir <- file.path(outputs_dir, "diagnostics")
exploration_dir <- file.path(outputs_dir, "exploration")
docs_dir <- here::here("docs")

dirs <- c(
  data_raw_dir, data_interim_dir, data_processed_dir, outputs_root_dir,
  outputs_dir, figures_dir, tables_dir, models_dir, diagnostics_dir,
  exploration_dir, docs_dir
)

invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
message("Outputs will be saved to: ", outputs_dir)


#### 4. Raw data paths ####

abundance_raw_path <- file.path(data_raw_dir, "2026.07.01_ArtificialReefs_MASTER.csv")
size_raw_path <- file.path(data_raw_dir, "2025.11.13_FishSize_MASTER.csv")


#### 5. Shared levels and lookups ####

reef_type_levels <- c("Natural", "Artificial")
period_levels <- c("Pre", "Post")
pair_levels <- c("Aow Mao", "No Name", "Sattakut")
feeding_guild_levels <- c("Herbivore", "Invertivore", "Mesopredator", "HTLP")
life_stage_levels <- c("juvenile", "adult")

site_lookup <- tibble::tribble(
  ~site,              ~type,        ~pair,
  "Aow Mao",          "Natural",    "Aow Mao",
  "Aow Mao Wreck",    "Artificial", "Aow Mao",
  "No Name Pinnacle", "Natural",    "No Name",
  "No Name Wreck",    "Artificial", "No Name",
  "Hin Pee Wee",      "Natural",    "Sattakut",
  "Sattakut",         "Artificial", "Sattakut"
)

deployment_lookup <- tibble::tribble(
  ~pair,      ~deployment_date,
  "Aow Mao",  as.Date("2023-09-01"),
  "No Name",  as.Date("2023-09-01"),
  "Sattakut", as.Date(NA)
)


#### 6. Colours and theme ####

reef_cols <- c("Artificial" = "#253494", "Natural" = "#66BFA6")

feeding_guild_cols <- c(
  "Herbivore" = "#66c2a4",
  "Invertivore" = "#41b6c4",
  "Mesopredator" = "#2c7fb8",
  "HTLP" = "#253494"
)

life_stage_cols <- c("juvenile" = "#66BFA6", "adult" = "#253494")

theme_clean <- theme_minimal(base_family = "Arial") +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_blank()
  )


#### 7. Helper functions ####

format_p <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ formatC(p, format = "f", digits = 3)
  )
}

model_export <- function(model, model_name, output_dir = tables_dir, sigfigs = 3) {
  fx_out <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    transmute(
      Effect = term,
      Estimate = signif(estimate, sigfigs),
      SE = signif(std.error, sigfigs),
      CI = paste0("[", signif(conf.low, sigfigs), ", ", signif(conf.high, sigfigs), "]"),
      p_value = format_p(p.value)
    )
  
  readr::write_csv(fx_out, file.path(output_dir, paste0(model_name, "_summarytable.csv")))
  invisible(fx_out)
}

save_model_summary <- function(model, model_name, output_dir = diagnostics_dir) {
  capture.output(summary(model), file = file.path(output_dir, paste0(model_name, "_summary_", analysis_date, ".txt")))
  invisible(NULL)
}

save_model_rds <- function(model, model_name, output_dir = models_dir) {
  saveRDS(model, file.path(output_dir, paste0(model_name, "_", analysis_date, ".rds")))
  invisible(NULL)
}


#### 8. Size-processing helpers and session record ####

source("size-fixes/A_matlookup.R")
source("size-fixes/B_size-bin-fix.R")

writeLines(capture.output(sessionInfo()), file.path(diagnostics_dir, paste0("session_info_", analysis_date, ".txt")))
message("Setup complete.")