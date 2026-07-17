### 02_EFFORT.R ###

summarise_effort <- function(fish_long, fish_counts, fish_size) {
  abundance_surveys <- fish_long %>%
    distinct(survey_id, site, pair, type, period, Date, n_observer_rows, n_observers, researchers)
  
  size_surveys <- fish_size %>%
    distinct(survey_id, site, pair, type, Date)
  
  size_stage_survey <- fish_size %>%
    group_by(survey_id, site, pair, type, Date, life_stage) %>%
    summarise(stage_count = sum(stage_Count), .groups = "drop")
  
  abundance_effort <- abundance_surveys %>%
    group_by(pair, type, period) %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      first_survey = min(Date),
      last_survey = max(Date),
      .groups = "drop"
    )
  
  abundance_effort_site <- abundance_surveys %>%
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
      total_individuals = sum(Count),
      first_survey = min(Date),
      last_survey = max(Date)
    )
  
  abundance_observer_summary <- abundance_surveys %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      surveys_with_multiple_observer_rows = sum(n_observer_rows > 1),
      median_observer_rows = median(n_observer_rows),
      max_observer_rows = max(n_observer_rows),
      median_observers = median(n_observers),
      max_observers = max(n_observers)
    )
  
  abundance_feeding_guild_totals <- fish_counts %>%
    group_by(pair, type, period, feeding_guild) %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      total_individuals = sum(Count),
      mean_count_per_survey = mean(Count),
      sd_count_per_survey = sd(Count),
      .groups = "drop"
    )
  
  abundance_effort_year <- abundance_surveys %>%
    mutate(year = format(Date, "%Y")) %>%
    count(year, type, name = "n_surveys") %>%
    pivot_wider(names_from = type, values_from = n_surveys, values_fill = 0)
  
  abundance_effort_2025 <- abundance_surveys %>%
    filter(format(Date, "%Y") == "2025") %>%
    count(pair, type, name = "n_surveys") %>%
    pivot_wider(names_from = type, values_from = n_surveys, values_fill = 0)
  
  size_effort <- size_surveys %>%
    group_by(pair, type) %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      first_survey = min(Date),
      last_survey = max(Date),
      .groups = "drop"
    )
  
  size_effort_site <- size_surveys %>%
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
      n_pairs = n_distinct(pair),
      n_species = n_distinct(Species),
      total_stage_individuals = sum(stage_Count),
      first_survey = min(Date),
      last_survey = max(Date)
    )
  
  size_life_stage_totals <- size_stage_survey %>%
    group_by(pair, type, life_stage) %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      total_individuals = sum(stage_count),
      mean_count_per_survey = mean(stage_count),
      sd_count_per_survey = sd(stage_count),
      .groups = "drop"
    )
  
  size_species_stage_totals <- fish_size %>%
    group_by(Species, Sci_Name, feeding_guild, life_stage) %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      total_individuals = sum(stage_Count),
      .groups = "drop"
    ) %>%
    arrange(Species, life_stage)
  
  effort_summaries <- list(
    abundance_effort = abundance_effort,
    abundance_effort_site = abundance_effort_site,
    abundance_dataset_summary = abundance_dataset_summary,
    abundance_observer_summary = abundance_observer_summary,
    abundance_feeding_guild_totals = abundance_feeding_guild_totals,
    abundance_effort_year = abundance_effort_year,
    abundance_effort_2025 = abundance_effort_2025,
    size_effort = size_effort,
    size_effort_site = size_effort_site,
    size_dataset_summary = size_dataset_summary,
    size_life_stage_totals = size_life_stage_totals,
    size_species_stage_totals = size_species_stage_totals
  )
  
  purrr::iwalk(
    effort_summaries,
    ~ readr::write_csv(.x, file.path(tables_dir, paste0(.y, "_", analysis_date, ".csv")))
  )
  
  message("Saved effort and dataset summaries to: ", tables_dir)
  effort_summaries
}


#### Run summaries ####

effort_summaries <- summarise_effort(fish_long, fish_counts, fish_size)

effort_summaries$abundance_dataset_summary
effort_summaries$abundance_effort
effort_summaries$abundance_effort_2025
effort_summaries$size_dataset_summary
effort_summaries$size_effort