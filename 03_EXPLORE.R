### 03_EXPLORE.R ###
### Formal pre-model exploration following Zuur et al. (2010) ###

SAVE_EXPLORATION <- TRUE

run_exploration <- function(
    fish_long,
    fish_counts,
    fish_size,
    deployment_lookup,
    save_outputs = SAVE_EXPLORATION
) {
  if (save_outputs) {
    report_path <- file.path(exploration_dir, paste0("exploration_report_", analysis_date, ".txt"))
    sink(report_path, split = TRUE)
    on.exit(sink(), add = TRUE)
  }
  
  heading <- function(x) cat("\n\n### ", x, " ###\n", sep = "")
  print_tbl <- function(x, n = 30) { print(x, n = n, width = Inf); invisible(x) }
  flag <- function(ok, pass, fail) { cat(if (isTRUE(ok)) "PASS: " else "FLAG: ", if (isTRUE(ok)) pass else fail, "\n", sep = ""); invisible(ok) }
  safe_ratio <- function(x, y) ifelse(is.na(y) | y == 0, NA_real_, x / y)
  
  brown_forsythe <- function(data, response, group_vars) {
    y <- data[[response]]
    g <- interaction(data[group_vars], drop = TRUE, lex.order = TRUE)
    keep <- complete.cases(y, g)
    y <- y[keep]
    g <- droplevels(g[keep])
    
    if (length(unique(g)) < 2 || any(table(g) < 2)) return(NA_real_)
    
    deviation <- abs(y - ave(y, g, FUN = median))
    anova(lm(deviation ~ g))[["Pr(>F)"]][1]
  }
  
  design_check <- function(data, formula) {
    mm <- model.matrix(formula, data)
    tibble(
      formula = paste(deparse(formula), collapse = ""),
      n_rows = nrow(mm),
      n_columns = ncol(mm),
      rank = qr(mm)$rank,
      rank_deficient = qr(mm)$rank < ncol(mm),
      condition_number = kappa(mm, exact = FALSE)
    )
  }
  
  survey_dates <- fish_long %>%
    distinct(survey_id, site, pair, type, period, Date, date_num, date_s, t_since) %>%
    arrange(pair, type, Date)
  
  total_counts <- fish_counts %>%
    group_by(survey_id, site, pair, type, period, Date, date_s, t_since) %>%
    summarise(Count = sum(Count), .groups = "drop")
  
  
  #### 1. Structural integrity ####
  
  heading("1. Structural integrity")
  
  survey_structure <- fish_long %>%
    distinct(survey_id, site, pair, type, period, Date) %>%
    count(survey_id, name = "n_metadata_combinations")
  
  missing_guilds <- fish_long %>% filter(is.na(feeding_guild)) %>% count(Species, sort = TRUE)
  invalid_counts <- fish_long %>% filter(is.na(Count) | Count < 0 | Count != floor(Count))
  species_duplicates <- fish_long %>% count(survey_id, Species, name = "n_rows") %>% filter(n_rows > 1)
  
  flag(all(survey_structure$n_metadata_combinations == 1),
       "Each survey_id maps to one site/date/type combination.",
       "Some survey_id values map to multiple metadata combinations.")
  
  flag(nrow(missing_guilds) == 0,
       "All species have feeding-guild assignments.",
       paste(nrow(missing_guilds), "species have missing feeding-guild assignments."))
  
  flag(nrow(invalid_counts) == 0,
       "Counts are non-negative integers with no missing values.",
       paste(nrow(invalid_counts), "invalid count rows detected."))
  
  flag(nrow(species_duplicates) == 0,
       "There is one cleaned row per survey and species.",
       paste(nrow(species_duplicates), "duplicated survey-species combinations detected."))
  
  
  #### 2. Explicit zeros and reconciliation ####
  
  heading("2. Explicit zeros and count reconciliation")
  
  missing_survey_species <- tidyr::crossing(
    survey_id = unique(fish_long$survey_id),
    Species = levels(fish_long$Species)
  ) %>%
    anti_join(fish_long %>% distinct(survey_id, Species), by = c("survey_id", "Species"))
  
  count_reconciliation <- fish_long %>%
    group_by(survey_id, feeding_guild) %>%
    summarise(species_sum = sum(Count), .groups = "drop") %>%
    full_join(
      fish_counts %>% select(survey_id, feeding_guild, guild_count = Count),
      by = c("survey_id", "feeding_guild")
    ) %>%
    mutate(difference = species_sum - guild_count)
  
  flag(nrow(missing_survey_species) == 0,
       "Every survey-species combination is represented explicitly.",
       paste(nrow(missing_survey_species), "survey-species combinations are absent."))
  
  flag(all(count_reconciliation$difference == 0, na.rm = TRUE) && !anyNA(count_reconciliation$difference),
       "Species counts reconcile exactly with feeding-guild counts.",
       "Species and feeding-guild totals do not reconcile.")
  
  
  #### 3. Observer aggregation ####
  
  heading("3. Observer aggregation")
  
  observer_summary <- fish_long %>%
    distinct(survey_id, n_observer_rows, n_observers, researchers) %>%
    summarise(
      n_surveys = n(),
      multiple_row_surveys = sum(n_observer_rows > 1),
      multiple_observer_surveys = sum(n_observers > 1),
      repeated_same_observer_surveys = sum(n_observer_rows > n_observers),
      median_observers = median(n_observers),
      maximum_observers = max(n_observers)
    )
  
  print_tbl(observer_summary)
  
  flag(observer_summary$repeated_same_observer_surveys == 0,
       "No surveys contain repeated rows from the same observer.",
       paste(observer_summary$repeated_same_observer_surveys, "surveys may contain repeated rows from the same observer."))
  
  
  #### 4. Temporal support and reef matching ####
  
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
      ar <- filter(.x, type == "Artificial")
      nr <- filter(.x, type == "Natural")
      if (!nrow(ar) || !nrow(nr)) return(tibble())
      
      bind_rows(
        tibble(
          focal_type = "Artificial",
          focal_date = ar$Date,
          nearest_gap_days = vapply(ar$Date, \(x) min(abs(as.numeric(x - nr$Date))), numeric(1))
        ),
        tibble(
          focal_type = "Natural",
          focal_date = nr$Date,
          nearest_gap_days = vapply(nr$Date, \(x) min(abs(as.numeric(x - ar$Date))), numeric(1))
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
  
  
  #### 5. Predictor support ####
  
  heading("5. Predictor collinearity and design support")
  
  temporal_collinearity <- survey_dates %>%
    group_by(pair) %>%
    summarise(
      n_surveys = n(),
      cor_date_t_since = if (sd(date_s) > 0 && sd(t_since) > 0) cor(date_s, t_since) else NA_real_,
      cor_date_period = if (n_distinct(period) > 1) cor(date_s, as.numeric(period == "Post")) else NA_real_,
      .groups = "drop"
    )
  
  design_diagnostics <- bind_rows(
    design_check(fish_counts, ~ type * pair * date_s + feeding_guild),
    design_check(fish_counts, ~ type * (date_s + period) + pair + feeding_guild),
    design_check(fish_counts, ~ type * pair * date_s + period + feeding_guild)
  )
  
  print_tbl(temporal_collinearity)
  print_tbl(design_diagnostics)
  
  flag(!any(design_diagnostics$rank_deficient),
       "Candidate fixed-effect design matrices are full rank.",
       "At least one candidate fixed-effect design matrix is rank deficient.")
  
  flag(all(design_diagnostics$condition_number < 30),
       "Candidate fixed-effect design matrices have acceptable condition numbers.",
       "At least one candidate design has a condition number above 30.")
  
  
  #### 6. Count distributions ####
  
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
  
  flag(!any(count_distribution$variance_to_mean > 1.5, na.rm = TRUE),
       "No strong raw overdispersion was detected.",
       paste(sum(count_distribution$variance_to_mean > 1.5, na.rm = TRUE), "groups have variance-to-mean ratios above 1.5."))
  
  flag(!any(count_distribution$zero_percent >= 50),
       "No group contains at least 50% zeroes.",
       paste(sum(count_distribution$zero_percent >= 50), "groups contain at least 50% zeroes."))
  
  
  #### 7. Raw variance heterogeneity ####
  
  heading("7. Raw variance heterogeneity")
  
  total_variance <- total_counts %>%
    mutate(log_count = log1p(Count)) %>%
    group_by(pair, type) %>%
    summarise(n = n(), raw_variance = var(Count), log_variance = var(log_count), .groups = "drop")
  
  guild_variance <- fish_counts %>%
    mutate(log_count = log1p(Count)) %>%
    group_by(pair, type, feeding_guild) %>%
    summarise(n = n(), raw_variance = var(Count), log_variance = var(log_count), .groups = "drop")
  
  print_tbl(total_variance)
  print_tbl(guild_variance)
  
  total_bf_p <- brown_forsythe(total_counts %>% mutate(log_count = log1p(Count)), "log_count", c("pair", "type"))
  guild_bf_p <- brown_forsythe(fish_counts %>% mutate(log_count = log1p(Count)), "log_count", c("pair", "type", "feeding_guild"))
  
  cat(
    "Brown-Forsythe p-value, log total abundance: ", format.pval(total_bf_p, digits = 3), "\n",
    "Brown-Forsythe p-value, log guild abundance: ", format.pval(guild_bf_p, digits = 3), "\n",
    sep = ""
  )
  
  
  #### 8. Extreme observations ####
  
  heading("8. Extreme observations")
  
  extreme_surveys <- total_counts %>%
    group_by(pair, type) %>%
    mutate(q1 = quantile(Count, 0.25), q3 = quantile(Count, 0.75), iqr = q3 - q1, extreme = Count > q3 + 3 * iqr) %>%
    ungroup() %>%
    filter(extreme) %>%
    select(survey_id, site, pair, type, Date, Count)
  
  flag(nrow(extreme_surveys) == 0,
       "No surveys exceed the upper 3-IQR threshold.",
       paste(nrow(extreme_surveys), "extreme surveys require biological verification."))
  
  if (nrow(extreme_surveys)) print_tbl(extreme_surveys)
  
  
  #### 9. Random effects and cell support ####
  
  heading("9. Random-effect and interaction support")
  
  random_effect_support <- bind_rows(
    survey_dates %>% count(site, name = "n") %>% summarise(grouping_factor = "site", levels = n(), minimum_n = min(n), median_n = median(n), maximum_n = max(n)),
    survey_dates %>% count(pair, name = "n") %>% summarise(grouping_factor = "pair", levels = n(), minimum_n = min(n), median_n = median(n), maximum_n = max(n)),
    fish_counts %>% count(survey_id, name = "n") %>% summarise(grouping_factor = "survey_id", levels = n(), minimum_n = min(n), median_n = median(n), maximum_n = max(n))
  )
  
  print_tbl(random_effect_support)
  
  weak_cells <- fish_counts %>%
    filter(!(pair == "Sattakut" & period == "Pre")) %>%
    count(pair, type, period, feeding_guild, name = "n") %>%
    filter(n < 5)
  
  flag(nrow(weak_cells) == 0,
       "All pair/type/period/guild cells contain at least five observations.",
       paste(nrow(weak_cells), "interaction cells contain fewer than five observations."))
  
  if (nrow(weak_cells)) print_tbl(weak_cells)
  
  
  #### 10. Diversity prerequisites ####
  
  heading("10. Diversity prerequisites")
  
  community_wide <- fish_long %>%
    select(survey_id, site, pair, type, period, Date, date_s, Species, Count) %>%
    pivot_wider(names_from = Species, values_from = Count, values_fill = 0)
  
  metadata_cols <- c("survey_id", "site", "pair", "type", "period", "Date", "date_s")
  community_matrix <- as.matrix(community_wide[setdiff(names(community_wide), metadata_cols)])
  
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
  
  cat(
    "Brown-Forsythe p-value, Shannon diversity: ",
    format.pval(brown_forsythe(diversity_data, "shannon", c("pair", "type")), digits = 3),
    "\nSpearman correlation, richness vs total abundance: ",
    round(cor(diversity_data$richness, diversity_data$total_abundance, method = "spearman"), 3),
    "\n",
    sep = ""
  )
  
  
  #### 11. Size-data support ####
  
  heading("11. Size-data support")
  
  size_support <- fish_size %>%
    distinct(survey_id, site, pair, type, Date, Visibility_m) %>%
    group_by(pair, type) %>%
    summarise(
      n_surveys = n_distinct(survey_id),
      first_date = min(Date),
      last_date = max(Date),
      missing_visibility_surveys = n_distinct(survey_id[is.na(Visibility_m)]),
      missing_visibility_percent = 100 * missing_visibility_surveys / n_surveys,
      .groups = "drop"
    )
  
  print_tbl(size_support)
  
  species_stage_support <- fish_size %>%
    group_by(Species, Sci_Name, feeding_guild) %>%
    summarise(
      occupied_surveys = n_distinct(survey_id[stage_Count > 0]),
      occupied_dates = n_distinct(Date[stage_Count > 0]),
      occupied_reef_types = n_distinct(type[stage_Count > 0]),
      juvenile_individuals = sum(stage_Count[life_stage == "juvenile"]),
      adult_individuals = sum(stage_Count[life_stage == "adult"]),
      juvenile_surveys = n_distinct(survey_id[life_stage == "juvenile" & stage_Count > 0]),
      adult_surveys = n_distinct(survey_id[life_stage == "adult" & stage_Count > 0]),
      .groups = "drop"
    ) %>%
    mutate(
      descriptive_support =
        occupied_surveys >= 10 &
        occupied_dates >= 5 &
        occupied_reef_types == 2 &
        juvenile_individuals >= 10 &
        adult_individuals >= 10 &
        juvenile_surveys >= 3 &
        adult_surveys >= 3
    ) %>%
    arrange(desc(descriptive_support), desc(occupied_surveys))
  
  print_tbl(species_stage_support)
  
  
  #### 12. Survey timeline ####
  
  heading("12. Survey timeline")
  
  timeline <- ggplot(survey_dates, aes(Date, site, shape = type)) +
    geom_point(size = 2, alpha = 0.75) +
    geom_vline(
      data = deployment_lookup %>% filter(!is.na(deployment_date)),
      aes(xintercept = deployment_date),
      linetype = 2,
      inherit.aes = FALSE
    ) +
    facet_wrap(~pair, scales = "free_y", ncol = 1) +
    scale_shape_manual(values = c("Artificial" = 16, "Natural" = 1)) +
    labs(x = "Survey date", y = NULL, shape = "Reef type") +
    theme_clean
  
  print(timeline)
  
  if (save_outputs) {
    ggsave(
      file.path(exploration_dir, paste0("survey_timeline_", analysis_date, ".png")),
      timeline,
      width = 10,
      height = 7,
      dpi = 300,
      bg = "white"
    )
    
    readr::write_csv(extreme_surveys, file.path(exploration_dir, paste0("extreme_surveys_", analysis_date, ".csv")))
    message("Saved exploration report, timeline, and extreme-survey list to: ", exploration_dir)
  }
  
  invisible(NULL)
}


#### Run exploration ####

run_exploration(fish_long, fish_counts, fish_size, deployment_lookup)