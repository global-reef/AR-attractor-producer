## 03_EXPLORE.R ####
## Formal pre-model exploration following Zuur et al. (2010).
## Runs inside a function so temporary objects are not retained.
## By default, results are printed only and nothing is saved.

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(vegan)

SAVE_EXPLORATION <- TRUE # or FALSE

run_exploration <- function(
    fish_long,
    fish_counts,
    fish_size,
    deployment_lookup,
    save_outputs = SAVE_EXPLORATION
) {
  
  #### Helpers ####
  
  explore_dir <- file.path(outputs_dir, "exploration")
  
  if (save_outputs) {
    dir.create(explore_dir, recursive = TRUE, showWarnings = FALSE)
    
    report_path <- file.path(
      explore_dir,
      paste0("exploration_report_", analysis_date, ".txt")
    )
    
    sink(report_path, split = TRUE)
    on.exit(sink(), add = TRUE)
  }
  
  heading <- function(x) {
    cat("\n\n", paste0("### ", x, " ###"), "\n", sep = "")
  }
  
  print_tbl <- function(x, n = 30) {
    print(x, n = n, width = Inf)
    invisible(x)
  }
  
  flag <- function(ok, pass, fail) {
    if (isTRUE(ok)) {
      cat("PASS: ", pass, "\n", sep = "")
    } else {
      cat("FLAG: ", fail, "\n", sep = "")
    }
    invisible(ok)
  }
  
  safe_ratio <- function(x, y) {
    ifelse(is.na(y) | y == 0, NA_real_, x / y)
  }
  
  brown_forsythe <- function(data, response, group_vars) {
    y <- data[[response]]
    g <- interaction(data[group_vars], drop = TRUE, lex.order = TRUE)
    
    keep <- complete.cases(y, g)
    y <- y[keep]
    g <- droplevels(g[keep])
    
    if (length(unique(g)) < 2 || any(table(g) < 2)) {
      return(NA_real_)
    }
    
    group_median <- ave(y, g, FUN = median)
    absolute_deviation <- abs(y - group_median)
    
    anova(lm(absolute_deviation ~ g))[["Pr(>F)"]][1]
  }
  
  design_check <- function(data, formula) {
    mm <- model.matrix(formula, data = data)
    qr_mm <- qr(mm)
    
    tibble(
      formula = paste(deparse(formula), collapse = ""),
      n_rows = nrow(mm),
      n_columns = ncol(mm),
      rank = qr_mm$rank,
      rank_deficient = qr_mm$rank < ncol(mm),
      condition_number = kappa(mm, exact = FALSE)
    )
  }
  
  #### Core survey objects ####
  
  survey_dates <- fish_long %>%
    distinct(
      survey_id, site, pair, type, period,
      Date, date_num, date_s, t_since
    ) %>%
    arrange(pair, type, Date)
  
  total_counts <- fish_counts %>%
    group_by(
      survey_id, site, pair, type, period,
      Date, date_s, t_since
    ) %>%
    summarise(Count = sum(Count), .groups = "drop")
  
  #### 1. Structural integrity ####
  
  heading("1. Structural integrity")
  
  survey_structure <- fish_long %>%
    distinct(survey_id, site, pair, type, period, Date) %>%
    count(survey_id, name = "n_metadata_combinations")
  
  missing_guilds <- fish_long %>%
    filter(is.na(feeding_guild)) %>%
    count(Species, sort = TRUE)
  
  invalid_counts <- fish_long %>%
    filter(is.na(Count) | Count < 0 | Count != floor(Count))
  
  species_duplicates <- fish_long %>%
    count(survey_id, Species, name = "n_rows") %>%
    filter(n_rows > 1)
  
  flag(
    all(survey_structure$n_metadata_combinations == 1),
    "Each survey_id maps to one site/date/type combination.",
    "Some survey_id values map to multiple metadata combinations."
  )
  
  flag(
    nrow(missing_guilds) == 0,
    "All species have feeding-guild assignments.",
    paste(nrow(missing_guilds), "species have missing feeding-guild assignments.")
  )
  
  flag(
    nrow(invalid_counts) == 0,
    "Counts are non-negative integers with no missing values.",
    paste(nrow(invalid_counts), "invalid cleaned count rows detected.")
  )
  
  flag(
    nrow(species_duplicates) == 0,
    "There is one cleaned row per survey and species.",
    paste(nrow(species_duplicates), "duplicated survey-species combinations detected.")
  )
  
  #### 2. Explicit zeros and count reconciliation ####
  
  heading("2. Explicit zeros and count reconciliation")
  
  survey_species_grid <- tidyr::crossing(
    survey_id = unique(fish_long$survey_id),
    Species = levels(fish_long$Species)
  )
  
  missing_survey_species <- survey_species_grid %>%
    anti_join(
      fish_long %>% distinct(survey_id, Species),
      by = c("survey_id", "Species")
    )
  
  count_reconciliation <- fish_long %>%
    group_by(survey_id, feeding_guild) %>%
    summarise(species_sum = sum(Count), .groups = "drop") %>%
    full_join(
      fish_counts %>%
        select(survey_id, feeding_guild, guild_count = Count),
      by = c("survey_id", "feeding_guild")
    ) %>%
    mutate(difference = species_sum - guild_count)
  
  flag(
    nrow(missing_survey_species) == 0,
    "Every survey-species combination is represented explicitly.",
    paste(
      nrow(missing_survey_species),
      "survey-species combinations are absent rather than explicit zeroes."
    )
  )
  
  flag(
    all(count_reconciliation$difference == 0, na.rm = TRUE) &&
      !anyNA(count_reconciliation$difference),
    "Species counts reconcile exactly with feeding-guild counts.",
    "Species and feeding-guild totals do not reconcile."
  )
  
  #### 3. Observer aggregation ####
  
  heading("3. Observer aggregation")
  
  observer_summary <- fish_long %>%
    distinct(
      survey_id, n_observer_rows,
      n_observers, researchers
    ) %>%
    summarise(
      n_surveys = n(),
      multiple_row_surveys = sum(n_observer_rows > 1),
      multiple_observer_surveys = sum(n_observers > 1),
      repeated_same_observer_surveys =
        sum(n_observer_rows > n_observers),
      median_observers = median(n_observers),
      maximum_observers = max(n_observers)
    )
  
  print_tbl(observer_summary)
  
  flag(
    observer_summary$repeated_same_observer_surveys == 0,
    "No surveys contain repeated rows from the same observer.",
    paste(
      observer_summary$repeated_same_observer_surveys,
      "surveys may contain repeated rows from the same observer."
    )
  )
  
  cat(
    "NOTE: Observer disagreement cannot be reconstructed after observer-level\n",
    "records have been averaged. It must be checked before aggregation in\n",
    "01A_CLEAN_abund.R.\n",
    sep = ""
  )
  
  #### 4. Temporal support and AR-NR matching ####
  
  heading("4. Temporal support and AR-NR matching")
  
  temporal_support <- survey_dates %>%
    group_by(pair, type, period) %>%
    summarise(
      n_surveys = n(),
      first_date = min(Date),
      last_date = max(Date),
      date_span_days = as.numeric(last_date - first_date),
      median_gap_days = if (n() > 1) median(diff(Date)) else NA_real_,
      maximum_gap_days = if (n() > 1) max(diff(Date)) else NA_real_,
      minimum_date_s = min(date_s),
      maximum_date_s = max(date_s),
      .groups = "drop"
    )
  
  print_tbl(temporal_support)
  
  nearest_date_match <- survey_dates %>%
    select(pair, type, Date) %>%
    group_by(pair) %>%
    group_modify(~{
      ar <- .x %>% filter(type == "Artificial")
      nr <- .x %>% filter(type == "Natural")
      
      if (nrow(ar) == 0 || nrow(nr) == 0) {
        return(tibble())
      }
      
      bind_rows(
        tibble(
          focal_type = "Artificial",
          focal_date = ar$Date,
          nearest_gap_days = vapply(
            ar$Date,
            function(x) min(abs(as.numeric(x - nr$Date))),
            numeric(1)
          )
        ),
        tibble(
          focal_type = "Natural",
          focal_date = nr$Date,
          nearest_gap_days = vapply(
            nr$Date,
            function(x) min(abs(as.numeric(x - ar$Date))),
            numeric(1)
          )
        )
      )
    }) %>%
    ungroup()
  
  match_summary <- nearest_date_match %>%
    group_by(pair, focal_type) %>%
    summarise(
      n_dates = n(),
      median_nearest_gap_days = median(nearest_gap_days),
      maximum_nearest_gap_days = max(nearest_gap_days),
      percent_within_7_days = mean(nearest_gap_days <= 7) * 100,
      percent_within_30_days = mean(nearest_gap_days <= 30) * 100,
      .groups = "drop"
    )
  
  print_tbl(match_summary)
  
  #### 5. Predictor collinearity and design support ####
  
  heading("5. Predictor collinearity and design support")
  
  temporal_collinearity <- survey_dates %>%
    group_by(pair) %>%
    summarise(
      n_surveys = n(),
      cor_date_t_since = if (
        sd(date_s) > 0 && sd(t_since) > 0
      ) cor(date_s, t_since) else NA_real_,
      cor_date_period = if (
        n_distinct(period) > 1
      ) cor(date_s, as.numeric(period == "Post")) else NA_real_,
      .groups = "drop"
    )
  
  print_tbl(temporal_collinearity)
  
  design_diagnostics <- bind_rows(
    design_check(
      fish_counts,
      ~ type * pair * date_s + feeding_guild
    ),
    design_check(
      fish_counts,
      ~ type * (date_s + period) + pair + feeding_guild
    ),
    design_check(
      fish_counts,
      ~ type * pair * date_s + period + feeding_guild
    )
  )
  
  print_tbl(design_diagnostics)
  
  flag(
    !any(design_diagnostics$rank_deficient),
    "Candidate fixed-effect design matrices are full rank.",
    "At least one candidate fixed-effect design matrix is rank deficient."
  )
  
  flag(
    all(design_diagnostics$condition_number < 30),
    "Candidate fixed-effect design matrices have acceptable condition numbers.",
    "At least one candidate design has a condition number above 30, indicating possible collinearity."
  )
  
  cat(
    "Interpretation guide: |r| > 0.7 is concerning; condition numbers above\n",
    "approximately 30 indicate potentially unstable fixed-effect estimation.\n",
    sep = ""
  )
  
  #### 6. Mean-variance relationship, zeros, and dispersion ####
  
  heading("6. Mean-variance relationship, zeros, and dispersion")
  
  count_distribution <- fish_counts %>%
    group_by(pair, type, feeding_guild) %>%
    summarise(
      n = n(),
      mean_count = mean(Count),
      variance = var(Count),
      variance_to_mean = safe_ratio(variance, mean_count),
      median_count = median(Count),
      maximum_count = max(Count),
      zero_percent = mean(Count == 0) * 100,
      .groups = "drop"
    )
  
  print_tbl(count_distribution)
  
  overdispersed_groups <- count_distribution %>%
    filter(variance_to_mean > 1.5)
  
  flag(
    nrow(overdispersed_groups) == 0,
    "No strong raw overdispersion was detected.",
    paste(
      nrow(overdispersed_groups),
      "pair/type/guild groups have variance-to-mean ratios above 1.5."
    )
  )
  
  high_zero_groups <- count_distribution %>%
    filter(zero_percent >= 50)
  
  flag(
    nrow(high_zero_groups) == 0,
    "No pair/type/guild group contains at least 50% zeroes.",
    paste(
      nrow(high_zero_groups),
      "pair/type/guild groups contain at least 50% zeroes."
    )
  )
  
  #### 7. Variance heterogeneity before modelling ####
  
  heading("7. Raw variance heterogeneity")
  
  total_variance <- total_counts %>%
    mutate(log_count = log1p(Count)) %>%
    group_by(pair, type) %>%
    summarise(
      n = n(),
      raw_variance = var(Count),
      log_variance = var(log_count),
      .groups = "drop"
    )
  
  guild_variance <- fish_counts %>%
    mutate(log_count = log1p(Count)) %>%
    group_by(pair, type, feeding_guild) %>%
    summarise(
      n = n(),
      raw_variance = var(Count),
      log_variance = var(log_count),
      .groups = "drop"
    )
  
  print_tbl(total_variance)
  print_tbl(guild_variance)
  
  total_bf_p <- brown_forsythe(
    total_counts %>% mutate(log_count = log1p(Count)),
    response = "log_count",
    group_vars = c("pair", "type")
  )
  
  guild_bf_p <- brown_forsythe(
    fish_counts %>% mutate(log_count = log1p(Count)),
    response = "log_count",
    group_vars = c("pair", "type", "feeding_guild")
  )
  
  cat(
    "Brown-Forsythe p-value, log total abundance: ",
    format.pval(total_bf_p, digits = 3), "\n",
    "Brown-Forsythe p-value, log guild abundance: ",
    format.pval(guild_bf_p, digits = 3), "\n",
    sep = ""
  )
  
  cat(
    "NOTE: Equal raw variance is not an assumption of a negative-binomial\n",
    "GLMM. These checks describe heterogeneity only. Formal variance adequacy,\n",
    "dispersion, zero inflation, and residual patterns must be assessed after\n",
    "model fitting using DHARMa and performance diagnostics.\n",
    sep = ""
  )
  
  #### 8. Extreme observations ####
  
  heading("8. Extreme observations")
  
  extreme_surveys <- total_counts %>%
    group_by(pair, type) %>%
    mutate(
      q1 = quantile(Count, 0.25),
      q3 = quantile(Count, 0.75),
      iqr = q3 - q1,
      extreme = Count > q3 + 3 * iqr
    ) %>%
    ungroup() %>%
    filter(extreme) %>%
    select(survey_id, site, pair, type, Date, Count)
  
  flag(
    nrow(extreme_surveys) == 0,
    "No surveys exceed the upper 3-IQR threshold.",
    paste(nrow(extreme_surveys), "extreme surveys require biological verification.")
  )
  
  if (nrow(extreme_surveys) > 0) {
    print_tbl(extreme_surveys)
  }
  
  #### 9. Random-effect and interaction support ####
  
  heading("9. Random-effect and interaction support")
  
  random_effect_support <- bind_rows(
    survey_dates %>%
      count(site, name = "n") %>%
      summarise(
        grouping_factor = "site",
        levels = n(),
        minimum_n = min(n),
        median_n = median(n),
        maximum_n = max(n)
      ),
    survey_dates %>%
      count(pair, name = "n") %>%
      summarise(
        grouping_factor = "pair",
        levels = n(),
        minimum_n = min(n),
        median_n = median(n),
        maximum_n = max(n)
      ),
    fish_counts %>%
      count(survey_id, name = "n") %>%
      summarise(
        grouping_factor = "survey_id",
        levels = n(),
        minimum_n = min(n),
        median_n = median(n),
        maximum_n = max(n)
      )
  )
  
  print_tbl(random_effect_support)
  
  cell_support <- fish_counts %>%
    filter(!(pair == "Sattakut" & period == "Pre")) %>%
    count(pair, type, period, feeding_guild, name = "n")
  
  weak_cells <- cell_support %>% filter(n < 5)
  
  flag(
    nrow(weak_cells) == 0,
    "All pair/type/period/guild cells contain at least five observations.",
    paste(
      nrow(weak_cells),
      "interaction cells contain fewer than five observations."
    )
  )
  
  if (nrow(weak_cells) > 0) {
    print_tbl(weak_cells)
  }
  
  #### 10. Diversity prerequisites and Gaussian variance checks ####
  
  heading("10. Diversity prerequisites")
  
  community_wide <- fish_long %>%
    select(
      survey_id, site, pair, type,
      period, Date, date_s, Species, Count
    ) %>%
    pivot_wider(
      names_from = Species,
      values_from = Count,
      values_fill = 0
    )
  
  species_columns <- setdiff(
    names(community_wide),
    c("survey_id", "site", "pair", "type", "period", "Date", "date_s")
  )
  
  community_matrix <- as.matrix(community_wide[species_columns])
  
  diversity_data <- community_wide %>%
    transmute(
      survey_id, site, pair, type, period, Date, date_s,
      total_abundance = rowSums(community_matrix),
      richness = rowSums(community_matrix > 0),
      shannon = vegan::diversity(community_matrix, index = "shannon")
    )
  
  diversity_summary <- diversity_data %>%
    group_by(pair, type) %>%
    summarise(
      n = n(),
      mean_shannon = mean(shannon),
      variance_shannon = var(shannon),
      minimum_shannon = min(shannon),
      maximum_shannon = max(shannon),
      zero_abundance_surveys = sum(total_abundance == 0),
      .groups = "drop"
    )
  
  print_tbl(diversity_summary)
  
  shannon_bf_p <- brown_forsythe(
    diversity_data,
    response = "shannon",
    group_vars = c("pair", "type")
  )
  
  cat(
    "Brown-Forsythe p-value, Shannon diversity: ",
    format.pval(shannon_bf_p, digits = 3), "\n",
    sep = ""
  )
  
  richness_abundance_cor <- cor(
    diversity_data$richness,
    diversity_data$total_abundance,
    method = "spearman"
  )
  
  cat(
    "Spearman correlation, richness vs total abundance: ",
    round(richness_abundance_cor, 3), "\n",
    sep = ""
  )
  
  #### 11. Size-data support ####
  
  heading("11. Size-data support")
  
  size_support <- fish_size %>%
    group_by(pair, type) %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      first_date = min(Date),
      last_date = max(Date),
      missing_visibility_percent =
        mean(is.na(Visibility_m)) * 100,
      .groups = "drop"
    )
  
  print_tbl(size_support)
 
   size_visibility <- fish_size %>%
    distinct(survey_id, site, pair, type, Date, Visibility_m) %>%
    group_by(pair, type) %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      missing_visibility_surveys =
        n_distinct(survey_id[is.na(Visibility_m)]),
      missing_visibility_percent =
        100 * missing_visibility_surveys / n_surveys,
      .groups = "drop"
    )
  print(size_visibility)
  
  species_stage_eligibility <- fish_size %>%
    group_by(Species, Sci_Name, feeding_guild) %>%
    summarise(
      occupied_surveys =
        n_distinct(survey_id[stage_Count > 0]),
      occupied_dates =
        n_distinct(Date[stage_Count > 0]),
      occupied_reef_types =
        n_distinct(type[stage_Count > 0]),
      juvenile_individuals =
        sum(stage_Count[life_stage == "juvenile"]),
      adult_individuals =
        sum(stage_Count[life_stage == "adult"]),
      juvenile_surveys =
        n_distinct(
          survey_id[
            life_stage == "juvenile" &
              stage_Count > 0
          ]
        ),
      adult_surveys =
        n_distinct(
          survey_id[
            life_stage == "adult" &
              stage_Count > 0
          ]
        ),
      .groups = "drop"
    ) %>%
    mutate(
      eligible =
        occupied_surveys >= 10 &
        occupied_dates >= 5 &
        occupied_reef_types == 2 &
        juvenile_individuals >= 10 &
        adult_individuals >= 10 &
        juvenile_surveys >= 3 &
        adult_surveys >= 3
    ) %>%
    arrange(desc(eligible), desc(occupied_surveys))
  
  cat(
    "Species eligible for formal juvenile/adult models: ",
    sum(species_stage_eligibility$eligible), " of ",
    nrow(species_stage_eligibility), "\n",
    sep = ""
  )
  
  print_tbl(species_stage_eligibility)
  
  #### 12. Timeline display ####
  
  heading("12. Survey timeline")
  
  timeline <- ggplot(
    survey_dates,
    aes(Date, site, shape = type)
  ) +
    geom_point(size = 2, alpha = 0.75) +
    geom_vline(
      data = deployment_lookup %>%
        filter(!is.na(deployment_date)),
      aes(xintercept = deployment_date),
      linetype = 2,
      inherit.aes = FALSE
    ) +
    facet_wrap(~pair, scales = "free_y", ncol = 1) +
    scale_shape_manual(
      values = c("Artificial" = 16, "Natural" = 1)
    ) +
    labs(
      x = "Survey date",
      y = NULL,
      shape = "Reef type"
    ) +
    theme_clean
  
  print(timeline)
  
  if (save_outputs) {
    ggsave(
      file.path(
        explore_dir,
        paste0("survey_timeline_", analysis_date, ".png")
      ),
      timeline,
      width = 10,
      height = 7,
      dpi = 300
    )
    
    cat(
      "\nSaved only:\n",
      "- exploration report\n",
      "- survey timeline\n",
      "to ", explore_dir, "\n",
      sep = ""
    )
  } else {
    cat(
      "\nExploration complete. Nothing was saved locally.\n",
      "Set SAVE_EXPLORATION <- TRUE to save only the report and timeline.\n",
      sep = ""
    )
  }
  
  invisible(NULL)
}


run_exploration(
  fish_long = fish_long,
  fish_counts = fish_counts,
  fish_size = fish_size,
  deployment_lookup = deployment_lookup
)



### 
#### 99. Verify extreme abundance surveys against raw records ####
read.csv(
  abundance_raw_path,
  stringsAsFactors = FALSE,
  strip.white = TRUE
) %>%
  as_tibble() %>%
  mutate(
    Date = as.Date(Date, "%m/%d/%Y"),
    Site = recode(Site, "No Name" = "No Name Pinnacle")
  ) %>%
  semi_join(
    tibble::tribble(
      ~Site,               ~Date,
      "Aow Mao",           as.Date("2024-09-16"),
      "No Name Pinnacle",  as.Date("2024-03-17"),
      "No Name Pinnacle",  as.Date("2024-03-18"),
      "No Name Wreck",     as.Date("2023-10-01"),
      "No Name Wreck",     as.Date("2024-08-01")
    ),
    by = c("Site", "Date")
  ) %>%
  pivot_longer(
    cols = any_of(c(
      "Parrotfish", "Rabbitfish", "Butterflyfish", "Angelfish",
      "Cleaner_Wrasse", "Batfish", "Thicklip", "Red_Breast",
      "Slingjaw", "Sweetlips", "Squirrel.Soldier", "Triggerfish",
      "Porcupine.Puffer", "Ray", "Brown_Stripe_Snapper",
      "Russels_Snapper", "lrg_Snapper", "Eel", "Trevally",
      "Emperorfish", "sml_Grouper", "lrg_Grouper", "Barracuda"
    )),
    names_to = "Species",
    values_to = "Count"
  ) %>%
  mutate(
    Count = as.numeric(Count),
    Count = replace_na(Count, 0)
  ) %>%
  group_by(Site, Date, Researcher, Species) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  group_by(Site, Date, Researcher) %>%
  mutate(observer_total = sum(Count)) %>%
  group_by(Site, Date, Species) %>%
  summarise(
    mean_count = mean(Count),
    min_count = min(Count),
    max_count = max(Count),
    n_observers = n_distinct(Researcher),
    mean_observer_total = mean(observer_total),
    .groups = "drop"
  ) %>%
  group_by(Site, Date) %>%
  mutate(
    contribution_pct = 100 * mean_count / sum(mean_count)
  ) %>%
  filter(mean_count > 0) %>%
  arrange(Site, Date, desc(mean_count)) %>%
  print(n = Inf, width = Inf)
