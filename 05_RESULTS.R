### 05_RESULTS.R ###

#### 1. Helpers ####

trajectory_label <- function(estimate, p_value) {
  case_when(
    p_value < 0.05 & estimate > 0 ~ "increase",
    p_value < 0.05 & estimate < 0 ~ "decrease",
    TRUE ~ "no significant change"
  )
}


#### 2. Abundance results ####

extract_abundance_results <- function(model, data, guild_adjust = "holm") {
  overall_emm <- emmeans::emtrends(model, ~type, var = "date_s", weights = "equal")
  
  overall_slopes <- summary(overall_emm, infer = c(TRUE, TRUE), adjust = "none") %>%
    as_tibble() %>%
    rename(slope = date_s.trend, p_value = p.value) %>%
    mutate(trajectory = trajectory_label(slope, p_value))
  
  overall_contrast <- contrast(
    overall_emm, method = list("Artificial - Natural" = c(-1, 1))
  ) %>%
    summary(infer = c(TRUE, TRUE), adjust = "none") %>%
    as_tibble() %>%
    rename(slope_difference = estimate, p_value = p.value)
  
  pair_emm <- emmeans::emtrends(model, ~type | pair, var = "date_s", weights = "equal")
  
  pair_slopes <- summary(update(pair_emm, by = NULL), infer = c(TRUE, TRUE), adjust = "none") %>%
    as_tibble() %>%
    rename(slope = date_s.trend, p_value = p.value) %>%
    mutate(trajectory = trajectory_label(slope, p_value))
  
  pair_contrasts <- contrast(
    pair_emm, method = list("Artificial - Natural" = c(-1, 1)), by = "pair"
  ) %>%
    update(by = NULL) %>%
    summary(infer = c(TRUE, TRUE), adjust = "none") %>%
    as_tibble() %>%
    rename(slope_difference = estimate, p_value = p.value)
  
  guild_emm <- emmeans::emtrends(model, ~type | pair * feeding_guild, var = "date_s")
  
  guild_slopes <- summary(update(guild_emm, by = NULL), infer = c(TRUE, TRUE), adjust = guild_adjust) %>%
    as_tibble() %>%
    rename(slope = date_s.trend, p_adjusted = p.value) %>%
    mutate(trajectory = trajectory_label(slope, p_adjusted))
  
  guild_contrasts <- contrast(
    guild_emm, method = list("Artificial - Natural" = c(-1, 1)),
    by = c("pair", "feeding_guild")
  ) %>%
    update(by = NULL) %>%
    summary(infer = c(TRUE, TRUE), adjust = guild_adjust) %>%
    as_tibble() %>%
    rename(slope_difference = estimate, p_adjusted = p.value)
  
  list(
    overall_slopes = overall_slopes,
    overall_contrast = overall_contrast,
    pair_slopes = pair_slopes,
    pair_contrasts = pair_contrasts,
    guild_slopes = guild_slopes,
    guild_contrasts = guild_contrasts
  )
}

abundance_results <- extract_abundance_results(abundance_final, abundance_data)


#### 3. Total abundance results ####

total_pair_emm <- emmeans::emtrends(total_abundance_final, ~type | pair, var = "date_s")

total_abundance_results <- list(
  pair_slopes = summary(update(total_pair_emm, by = NULL), infer = c(TRUE, TRUE), adjust = "none") %>%
    as_tibble() %>%
    rename(slope = date_s.trend, p_value = p.value) %>%
    mutate(trajectory = trajectory_label(slope, p_value)),
  
  pair_contrasts = contrast(
    total_pair_emm, method = list("Artificial - Natural" = c(-1, 1)), by = "pair"
  ) %>%
    update(by = NULL) %>%
    summary(infer = c(TRUE, TRUE), adjust = "none") %>%
    as_tibble() %>%
    rename(slope_difference = estimate, p_value = p.value)
)


#### 4. Shannon diversity results ####

diversity_slope_emm <- emmeans::emtrends(shannon_final, ~type, var = "date_s", weights = "equal")

diversity_slopes <- summary(diversity_slope_emm, infer = c(TRUE, TRUE), adjust = "none") %>%
  as_tibble() %>%
  rename(slope = date_s.trend, p_value = p.value) %>%
  mutate(trajectory = trajectory_label(slope, p_value))

diversity_slope_contrast <- contrast(
  diversity_slope_emm, method = list("Artificial - Natural" = c(-1, 1))
) %>%
  summary(infer = c(TRUE, TRUE), adjust = "none") %>%
  as_tibble() %>%
  rename(slope_difference = estimate, p_value = p.value)

diversity_period_emm <- emmeans(shannon_final, ~type * period, weights = "equal")
diversity_period_contrast <- contrast(
  diversity_period_emm,
  interaction = c("revpairwise", "revpairwise")
) %>%
  summary(infer = c(TRUE, TRUE), adjust = "none") %>%
  as_tibble()

diversity_time_values <- c(start = min(diversity_data$date_s), end = max(diversity_data$date_s))

diversity_type_emm <- emmeans(
  shannon_final, ~type | date_s,
  at = list(date_s = unname(diversity_time_values), period = "Post"),
  weights = "equal"
)

diversity_type_estimates <- summary(update(diversity_type_emm, by = NULL), infer = c(TRUE, TRUE), adjust = "none") %>%
  as_tibble() %>%
  mutate(timepoint = if_else(near(date_s, diversity_time_values["start"]), "Start", "End")) %>%
  rename(estimated_shannon = emmean, p_value = p.value)

diversity_type_contrasts <- contrast(
  diversity_type_emm, method = list("Artificial - Natural" = c(-1, 1)), by = "date_s"
) %>%
  summary(infer = c(TRUE, TRUE), adjust = "none") %>%
  as_tibble() %>%
  mutate(timepoint = if_else(near(date_s, diversity_time_values["start"]), "Start", "End")) %>%
  rename(diversity_difference = estimate, p_value = p.value)

diversity_results <- list(
  slopes = diversity_slopes,
  slope_contrast = diversity_slope_contrast,
  period_contrast = diversity_period_contrast,
  type_estimates = diversity_type_estimates,
  type_contrasts = diversity_type_contrasts
)


#### 5. Community life-stage results ####

juvenile_tests <- emmeans::joint_tests(juvenile_final)

stage_emm <- emmeans(juvenile_final, ~life_stage | type * pair, type = "response")

stage_predictions <- summary(stage_emm, infer = c(TRUE, TRUE), adjust = "none") %>%
  as_tibble()

stage_proportions <- stage_predictions %>%
  select(type, pair, life_stage, response) %>%
  pivot_wider(names_from = life_stage, values_from = response) %>%
  mutate(total = juvenile + adult, juvenile_prop = juvenile / total, adult_prop = adult / total)

reef_emm <- emmeans(juvenile_final, ~type | life_stage * pair, type = "response")

stage_contrasts <- contrast(
  reef_emm, method = list("Artificial / Natural" = c(-1, 1)),
  by = c("life_stage", "pair")
) %>%
  summary(infer = c(TRUE, TRUE), adjust = "none") %>%
  as_tibble()

life_stage_contrasts <- contrast(
  stage_emm, method = list("Adult / Juvenile" = c(-1, 1)),
  by = c("type", "pair")
) %>%
  summary(infer = c(TRUE, TRUE), adjust = "none") %>%
  as_tibble()

stage_grid_pair <- emmeans(juvenile_final, ~life_stage * type | pair)

proportion_contrasts_pair <- contrast(
  stage_grid_pair, interaction = c("revpairwise", "revpairwise"), by = "pair"
) %>%
  summary(infer = c(TRUE, TRUE), type = "response", adjust = "none") %>%
  as_tibble()

stage_grid_overall <- emmeans(juvenile_final, ~life_stage * type)

proportion_contrast_overall <- contrast(
  stage_grid_overall, interaction = c("revpairwise", "revpairwise")
) %>%
  summary(infer = c(TRUE, TRUE), type = "response", adjust = "none") %>%
  as_tibble()

juvenile_results <- list(
  tests = juvenile_tests,
  predictions = stage_predictions,
  proportions = stage_proportions,
  reef_contrasts = stage_contrasts,
  stage_contrasts = life_stage_contrasts,
  pair_composition = proportion_contrasts_pair,
  overall_composition = proportion_contrast_overall
)


#### 6. Species-specific results ####

extract_species_results <- function(models, adjust = "holm") {
  omnibus <- imap_dfr(models, \(model, species) {
    car::Anova(model, type = 3) %>%
      as.data.frame() %>%
      rownames_to_column("term") %>%
      transmute(
        Species = species, term, chisq = Chisq, df = Df,
        p_value = `Pr(>Chisq)`,
        significance = if_else(p_value < 0.05, "Significant", "Non-significant")
      )
  })
  
  emm <- map(models, ~ emmeans(.x, ~life_stage * type * pair))
  
  predicted_counts <- imap_dfr(emm, \(x, species) {
    summary(x, type = "response", infer = c(TRUE, TRUE)) %>%
      as.data.frame() %>%
      mutate(Species = species, .before = 1)
  })
  
  extract_contrast <- function(x, species, by = NULL, interaction = NULL) {
    if (is.null(interaction)) {
      out <- contrast(x, method = "revpairwise", by = by, adjust = adjust)
    } else {
      out <- contrast(x, interaction = interaction, by = by, adjust = adjust)
    }
    
    summary(out, type = "response", infer = c(TRUE, TRUE)) %>%
      as.data.frame() %>%
      mutate(
        Species = species,
        significance = if_else(p.value < 0.05, "Significant", "Non-significant"),
        .before = 1
      )
  }
  
  reef_contrasts <- imap_dfr(emm, ~ extract_contrast(.x, .y, by = c("life_stage", "pair")))
  stage_contrasts <- imap_dfr(emm, ~ extract_contrast(.x, .y, by = c("type", "pair")))
  composition_contrasts <- imap_dfr(
    emm, ~ extract_contrast(.x, .y, by = "pair", interaction = c("revpairwise", "revpairwise"))
  )
  
  response_col <- intersect(c("response", "rate", "prob", "emmean"), names(predicted_counts))[1]
  if (is.na(response_col)) stop("Could not identify the predicted-response column.")
  
  juvenile_proportions <- predicted_counts %>%
    select(Species, life_stage, type, pair, predicted_count = all_of(response_col)) %>%
    pivot_wider(names_from = life_stage, values_from = predicted_count) %>%
    mutate(juvenile_proportion = juvenile / (juvenile + adult))
  
  list(
    omnibus = omnibus,
    predicted_counts = predicted_counts,
    reef_contrasts = reef_contrasts,
    stage_contrasts = stage_contrasts,
    juvenile_proportions = juvenile_proportions,
    composition_contrasts = composition_contrasts
  )
}

species_results <- extract_species_results(species_models)


#### 7. Save final result tables ####

result_tables <- c(
  abundance_results,
  list(
    total_pair_slopes = total_abundance_results$pair_slopes,
    total_pair_contrasts = total_abundance_results$pair_contrasts,
    diversity_slopes = diversity_results$slopes,
    diversity_slope_contrast = diversity_results$slope_contrast,
    diversity_period_contrast = diversity_results$period_contrast,
    diversity_type_estimates = diversity_results$type_estimates,
    diversity_type_contrasts = diversity_results$type_contrasts,
    juvenile_tests = as_tibble(juvenile_results$tests),
    juvenile_predictions = juvenile_results$predictions,
    juvenile_proportions = juvenile_results$proportions,
    juvenile_reef_contrasts = juvenile_results$reef_contrasts,
    juvenile_stage_contrasts = juvenile_results$stage_contrasts,
    juvenile_pair_composition = juvenile_results$pair_composition,
    juvenile_overall_composition = juvenile_results$overall_composition,
    species_omnibus = species_results$omnibus,
    species_predicted_counts = species_results$predicted_counts,
    species_reef_contrasts = species_results$reef_contrasts,
    species_stage_contrasts = species_results$stage_contrasts,
    species_juvenile_proportions = species_results$juvenile_proportions,
    species_composition_contrasts = species_results$composition_contrasts,
    species_model_summary = species_model_summary
  )
)

iwalk(result_tables, ~ readr::write_csv(.x, file.path(tables_dir, paste0(.y, "_", analysis_date, ".csv"))))
message("Saved final result tables to: ", tables_dir)
