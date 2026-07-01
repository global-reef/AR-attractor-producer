## 02_EFFORT.R ####

library(dplyr)
library(tidyr)

summarise_effort <- function(fish_long, fish_counts, fish_size) {
  
  ## Abundance survey effort ##
  abundance_effort <- fish_long %>%
    distinct(site, pair, type, period, Date, survey_id) %>%
    group_by(pair, type, period) %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      first_survey = min(Date),
      last_survey = max(Date),
      .groups = "drop"
    )
  
  abundance_effort_site <- fish_long %>%
    distinct(site, pair, type, period, Date, survey_id) %>%
    group_by(site, pair, type, period) %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      first_survey = min(Date),
      last_survey = max(Date),
      .groups = "drop"
    )
  
  abundance_dataset_summary <- fish_long %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      n_sites = n_distinct(site),
      n_pairs = n_distinct(pair),
      n_species = n_distinct(Species),
      total_individuals = sum(Count, na.rm = TRUE),
      first_survey = min(Date),
      last_survey = max(Date)
    )
  
  abundance_observer_summary <- fish_long %>%
    distinct(survey_id, site, pair, type, period, Date,
             n_observer_rows, n_observers, researchers) %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      surveys_with_multiple_observer_rows = sum(n_observer_rows > 1, na.rm = TRUE),
      median_observer_rows = median(n_observer_rows, na.rm = TRUE),
      max_observer_rows = max(n_observer_rows, na.rm = TRUE),
      median_observers = median(n_observers, na.rm = TRUE),
      max_observers = max(n_observers, na.rm = TRUE)
    )
  
  abundance_feeding_guild_totals <- fish_counts %>%
    group_by(pair, type, period, feeding_guild) %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      total_individuals = sum(Count, na.rm = TRUE),
      mean_count_per_survey = mean(Count, na.rm = TRUE),
      sd_count_per_survey = sd(Count, na.rm = TRUE),
      .groups = "drop"
    )
  
  ## Size survey effort ##
  size_effort <- fish_size %>%
    distinct(site, pair, type, Date, survey_id) %>%
    group_by(pair, type) %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      first_survey = min(Date),
      last_survey = max(Date),
      .groups = "drop"
    )
  
  size_effort_site <- fish_size %>%
    distinct(site, pair, type, Date, survey_id) %>%
    group_by(site, pair, type) %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      first_survey = min(Date),
      last_survey = max(Date),
      .groups = "drop"
    )
  
  size_dataset_summary <- fish_size %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      n_sites = n_distinct(site),
      n_species = n_distinct(Species),
      total_stage_individuals = sum(stage_Count, na.rm = TRUE),
      first_survey = min(Date),
      last_survey = max(Date)
    )
  
  size_life_stage_totals <- fish_size %>%
    group_by(pair, type, life_stage) %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      total_individuals = sum(stage_Count, na.rm = TRUE),
      mean_count_per_survey = mean(stage_Count, na.rm = TRUE),
      sd_count_per_survey = sd(stage_Count, na.rm = TRUE),
      .groups = "drop"
    )
  
  size_species_stage_totals <- fish_size %>%
    group_by(Species, Sci_Name, feeding_guild, life_stage) %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      total_individuals = sum(stage_Count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Species, life_stage)
  
  ## Save publication-useful tables ##
  write.csv(
    abundance_effort,
    file.path(tables_dir, paste0("abundance_effort_pair_type_period_", analysis_date, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    abundance_effort_site,
    file.path(tables_dir, paste0("abundance_effort_site_", analysis_date, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    abundance_dataset_summary,
    file.path(tables_dir, paste0("abundance_dataset_summary_", analysis_date, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    abundance_observer_summary,
    file.path(tables_dir, paste0("abundance_observer_summary_", analysis_date, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    abundance_feeding_guild_totals,
    file.path(tables_dir, paste0("abundance_feeding_guild_totals_", analysis_date, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    size_effort,
    file.path(tables_dir, paste0("size_effort_pair_type_", analysis_date, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    size_effort_site,
    file.path(tables_dir, paste0("size_effort_site_", analysis_date, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    size_dataset_summary,
    file.path(tables_dir, paste0("size_dataset_summary_", analysis_date, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    size_life_stage_totals,
    file.path(tables_dir, paste0("size_life_stage_totals_", analysis_date, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    size_species_stage_totals,
    file.path(tables_dir, paste0("size_species_stage_totals_", analysis_date, ".csv")),
    row.names = FALSE
  )
  
  ## Return summaries ##
  effort_summaries <- list(
    abundance_effort = abundance_effort,
    abundance_effort_site = abundance_effort_site,
    abundance_dataset_summary = abundance_dataset_summary,
    abundance_observer_summary = abundance_observer_summary,
    abundance_feeding_guild_totals = abundance_feeding_guild_totals,
    size_effort = size_effort,
    size_effort_site = size_effort_site,
    size_dataset_summary = size_dataset_summary,
    size_life_stage_totals = size_life_stage_totals,
    size_species_stage_totals = size_species_stage_totals
  )
  
  message("Saved effort and dataset summary tables.")
  
  return(effort_summaries)
}

## Run summaries ##
effort_summaries <- summarise_effort(
  fish_long = fish_long,
  fish_counts = fish_counts,
  fish_size = fish_size
)

## Check key summaries ##
effort_summaries$abundance_effort
effort_summaries$abundance_dataset_summary
effort_summaries$size_effort
effort_summaries$size_dataset_summary