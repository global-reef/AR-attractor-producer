### 00_SETUP.R ###

#### 1. Analysis metadata ####

analysis_date <- "2026.07.01"

message("Running AR-attractor-producer analysis: ", analysis_date)


#### 2. Packages ####

required_packages <- c(
  # data handling
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "forcats",
  "lubridate",
  "tibble",
  "purrr",
  
  # modelling
  "glmmTMB",
  "lme4",
  "lmerTest",
  "emmeans",
  "broom.mixed",
  "performance",
  "DHARMa",
  "splines",
  "mgcv",
  
  # diversity / community metrics
  "vegan",
  
  # plotting
  "ggplot2",
  "patchwork",
  "ggeffects",
  
  # project helpers
  "here"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  stop(
    "These packages are missing. Install them before running the analysis:\n",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  lapply(required_packages, library, character.only = TRUE)
})


#### 3. Project paths ####
#### 3. Project paths ####

project_dir <- here::here()

data_raw_dir       <- here::here("data", "raw")
data_interim_dir   <- here::here("data", "interim")
data_processed_dir <- here::here("data", "processed")

# Create a new dated analysis folder for each analysis version.
outputs_root_dir <- here::here("outputs")
outputs_dir <- file.path(
  outputs_root_dir,
  paste0("Analysis_", analysis_date)
)

figures_dir     <- file.path(outputs_dir, "figures")
tables_dir      <- file.path(outputs_dir, "tables")
models_dir      <- file.path(outputs_dir, "model_objects")
diagnostics_dir <- file.path(outputs_dir, "diagnostics")

docs_dir <- here::here("docs")

dir.create(data_raw_dir,       recursive = TRUE, showWarnings = FALSE)
dir.create(data_interim_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(data_processed_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(outputs_root_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(outputs_dir,      recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir,      recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir,       recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir,       recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostics_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(docs_dir,         recursive = TRUE, showWarnings = FALSE)

message("Outputs will be saved to: ", outputs_dir)
#### 4. Raw data file paths ####
abundance_raw_path <- file.path(
  data_raw_dir,
  "2026.07.01_ArtificialReefs_MASTER.csv"
)

size_raw_path <- file.path(
  data_raw_dir,
  "2025.11.13_FishSize_MASTER.csv"
)


#### 5. Shared factor levels ####

reef_type_levels <- c("Artificial", "Natural")

period_levels <- c("Pre", "Post")

feeding_guild_levels <- c(
  "Herbivore",
  "Invertivore",
  "Mesopredator",
  "HTLP"
)

life_stage_levels <- c("juvenile", "adult")


#### 6. Site lookup ####

site_lookup <- tibble::tribble(
  ~site,               ~type,         ~pair,
  "Aow Mao",           "Natural",     "Aow Mao",
  "Aow Mao Wreck",     "Artificial",  "Aow Mao",
  "No Name Pinnacle",  "Natural",     "No Name",
  "No Name Wreck",     "Artificial",  "No Name",
  "Hin Pee Wee",       "Natural",     "Sattakut",
  "Sattakut",          "Artificial",  "Sattakut"
)


#### 7. Deployment dates ####

deployment_lookup <- tibble::tribble(
  ~pair,       ~deployment_date,
  "Aow Mao",   as.Date("2023-09-01"),
  "No Name",   as.Date("2023-09-01"),
  "Sattakut",  as.Date(NA)
)


#### 8. Colour palettes ####

reef_cols <- c(
  "Artificial" = "#253494",
  "Natural"    = "#66BFA6"
)

feeding_guild_cols <- c(
  "Herbivore"       = "#66c2a4",
  "Invertivore" = "#41b6c4",
  "Mesopredator"= "#2c7fb8",
  "HTLP"         = "#253494"
)

life_stage_cols <- c(
  "juvenile" = "#66BFA6",
  "adult"    = "#253494"
)


#### 9. Plot theme ####

theme_clean <- ggplot2::theme_minimal(base_family = "Arial") +
  ggplot2::theme(
    legend.position = "right",
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = "white", colour = NA),
    plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
    plot.title       = ggplot2::element_blank()
  )


#### 10. Helper functions ####

format_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3))
  )
}


model_export <- function(model, model_name, output_dir = tables_dir, sigfigs = 3) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  fx <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  
  fx_out <- fx %>%
    dplyr::transmute(
      Effect   = term,
      Estimate = signif(estimate, sigfigs),
      SE       = signif(std.error, sigfigs),
      CI       = paste0(
        "[",
        signif(conf.low, sigfigs),
        ", ",
        signif(conf.high, sigfigs),
        "]"
      ),
      p_value = format_p(p.value)
    )
  
  readr::write_csv(
    fx_out,
    file.path(output_dir, paste0(model_name, "_summarytable.csv"))
  )
  
  invisible(fx_out)
}


save_model_summary <- function(model, model_name, output_dir = diagnostics_dir) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  capture.output(
    summary(model),
    file = file.path(output_dir, paste0(model_name, "_summary_", analysis_date, ".txt"))
  )
  
  invisible(NULL)
}


save_model_rds <- function(model, model_name, output_dir = models_dir) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  saveRDS(
    model,
    file = file.path(output_dir, paste0(model_name, "_", analysis_date, ".rds"))
  )
  
  invisible(NULL)
}


#### 11. Session info ####

writeLines(
  capture.output(sessionInfo()),
  con = file.path(diagnostics_dir, paste0("session_info_", analysis_date, ".txt"))
)
source("size-fixes/A_matlookup.R")
source("size-fixes/B_size-bin-fix.R")

message("Setup complete.")


