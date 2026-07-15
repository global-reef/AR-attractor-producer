#### 3. Shannon diversity trajectories ####

community_wide <- fish_long %>%
  select(
    survey_id, site, pair, type, period,
    Date, date_s, Species, Count
  ) %>%
  pivot_wider(
    names_from = Species,
    values_from = Count,
    values_fill = 0
  )

species_columns <- setdiff(
  names(community_wide),
  c(
    "survey_id", "site", "pair", "type",
    "period", "Date", "date_s"
  )
)

shannon_data <- community_wide %>%
  mutate(
    shannon = vegan::diversity(
      as.matrix(pick(all_of(species_columns))),
      index = "shannon"
    )
  ) %>%
  select(
    survey_id, site, pair, type,
    period, Date, date_s, shannon
  )

shannon_model <- lmer(
  shannon ~ type * (date_s + period) +
    (1 + date_s || site),
  data = shannon_data
)

cat("\n### Shannon diversity model ###\n")
print(summary(shannon_model))
print(performance::check_collinearity(shannon_model))
print(performance::check_singularity(shannon_model))
plot(performance::check_model(shannon_model))


#### 4. Community juvenile-adult patterns ####

size_dates <- fish_size %>%
  distinct(survey_id, Date) %>%
  mutate(
    date_num = as.numeric(Date - min(Date)),
    date_s = as.numeric(scale(date_num))
  )

community_stage_data <- fish_size %>%
  group_by(
    survey_id, site, pair, type,
    Date, life_stage
  ) %>%
  summarise(
    Count = sum(stage_Count),
    .groups = "drop"
  ) %>%
  complete(
    nesting(survey_id, site, pair, type, Date),
    life_stage,
    fill = list(Count = 0)
  ) %>%
  left_join(size_dates, by = c("survey_id", "Date"))

community_stage_model <- glmmTMB(
  Count ~ life_stage * type * pair +
    date_s +
    (1 | survey_id),
  family = nbinom2,
  data = community_stage_data
)

cat("\n### Community juvenile-adult model ###\n")
print(summary(community_stage_model))
community_stage_residuals <- check_glmm(
  community_stage_model,
  community_stage_data,
  time = community_stage_data$date_num,
  group = community_stage_data$survey_id
)


#### 5. Target-taxa juvenile-adult patterns ####

missing_target_taxa <- setdiff(
  target_taxa,
  unique(as.character(fish_size$Species))
)

if (length(missing_target_taxa) > 0) {
  warning(
    "Target taxa absent from fish_size: ",
    paste(missing_target_taxa, collapse = ", "),
    call. = FALSE
  )
}

target_stage_data <- fish_size %>%
  filter(Species %in% target_taxa) %>%
  group_by(
    survey_id, site, pair, type,
    Date, Species, life_stage
  ) %>%
  summarise(
    Count = sum(stage_Count),
    .groups = "drop"
  ) %>%
  complete(
    nesting(survey_id, site, pair, type, Date),
    Species,
    life_stage,
    fill = list(Count = 0)
  ) %>%
  left_join(size_dates, by = c("survey_id", "Date")) %>%
  mutate(
    Species = droplevels(factor(Species)),
    survey_taxon_id = interaction(
      survey_id, Species,
      drop = TRUE
    )
  )

target_stage_model <- glmmTMB(
  Count ~ Species * life_stage * type +
    pair + date_s +
    (1 | survey_taxon_id),
  family = nbinom2,
  data = target_stage_data
)

cat("\n### Target-taxa juvenile-adult model ###\n")
print(summary(target_stage_model))
target_stage_residuals <- check_glmm(
  target_stage_model,
  target_stage_data,
  time = target_stage_data$date_num,
  group = target_stage_data$survey_id
)


#### 6. Essential model estimates ####

abundance_slopes <- emtrends(
  abundance_final,
  ~ feeding_guild * type | pair,
  var = "date_s"
)

shannon_slopes <- emtrends(
  shannon_model,
  ~ type,
  var = "date_s"
)

community_stage_contrasts <- emmeans(
  community_stage_model,
  ~ life_stage * type | pair,
  type = "response"
)

target_stage_contrasts <- emmeans(
  target_stage_model,
  ~ life_stage * type | Species,
  type = "response"
)

cat("\n### Abundance temporal slopes ###\n")
print(abundance_slopes)

cat("\n### Shannon temporal slopes ###\n")
print(shannon_slopes)

cat("\n### Community juvenile-adult estimates ###\n")
print(community_stage_contrasts)

cat("\n### Target-taxa juvenile-adult estimates ###\n")
print(target_stage_contrasts)


#### 7. Save final model objects only ####

saveRDS(
  abundance_final,
  file.path(models_dir, paste0("abundance_final_", analysis_date, ".rds"))
)

saveRDS(
  shannon_model,
  file.path(models_dir, paste0("shannon_final_", analysis_date, ".rds"))
)

saveRDS(
  community_stage_model,
  file.path(models_dir, paste0("community_stage_final_", analysis_date, ".rds"))
)

saveRDS(
  target_stage_model,
  file.path(models_dir, paste0("target_stage_final_", analysis_date, ".rds"))
)

message("Four final model objects saved.")










#### trying juvenile modelling 3.4 #### 

## 4. Species-standardised life-stage abundance ####

juvenile_results <- local({
  
  #### 4.1 Analysis data ####
  
  survey_frame <- fish_size %>%
    distinct(survey_id, site, pair, type, Date)
  
  species_frame <- fish_size %>%
    distinct(Species)
  
  stage_counts <- fish_size %>%
    group_by(
      survey_id, site, pair, type, Date,
      Species, life_stage) %>%
    summarise(
      stage_count = sum(stage_Count),
      .groups = "drop")
  
  juvenile_data <- tidyr::crossing(
    survey_frame,
    species_frame,
    life_stage = c("juvenile", "adult")) %>%
    left_join(
      stage_counts,
      by = c(
        "survey_id", "site", "pair", "type", "Date",
        "Species", "life_stage")) %>%
    mutate(
      stage_count = as.integer(replace_na(stage_count, 0)),
      life_stage = factor(
        life_stage,
        levels = c("juvenile", "adult")),
      type = factor(
        type,
        levels = c("Natural", "Artificial")),
      pair = factor(pair),
      Species = factor(Species),
      survey_id = factor(survey_id))
  
  cat("\n### Life-stage count summary ###\n")
  juvenile_data %>%
    summarise(
      n_rows = n(),
      n_surveys = n_distinct(survey_id),
      n_species = n_distinct(Species),
      n_zero = sum(stage_count == 0),
      zero_proportion = mean(stage_count == 0),
      total_juveniles = sum(
        stage_count[life_stage == "juvenile"]),
      total_adults = sum(
        stage_count[life_stage == "adult"])) %>%
    print()
  
  
  #### 4.2 Distribution comparison ####
  
  poisson_model <- glmmTMB(
    stage_count ~ life_stage * type * pair +
      (1 | Species) + (1 | survey_id),
    family = poisson,
    data = juvenile_data)
  
  nbinom1_model <- update(
    poisson_model,
    family = nbinom1)
  
  nbinom2_model <- update(
    poisson_model,
    family = nbinom2)
  
  distribution_comparison <- AIC(
    poisson_model,
    nbinom1_model,
    nbinom2_model)
  
  cat("\n### Distribution comparison ###\n")
  print(distribution_comparison)
  
  
  #### 4.3 Primary conditional model ####
  
  # Fixed a priori from the ecological question:
  # Do juvenile and adult counts differ between ARs and NRs,
  # and does that pattern vary among site pairs?
  #
  # Species is included as a random intercept to control for
  # baseline abundance differences among taxa.
  # Survey is included for repeated species-stage observations.
  
  conditional_model <- nbinom2_model
  
  
  #### 4.4 Dispersion and zero-inflation structure ####
  
  candidate_models <- list(
    constant = conditional_model,
    dispersion_type = update(
      conditional_model,
      dispformula = ~type),
    dispersion_pair = update(
      conditional_model,
      dispformula = ~pair),
    dispersion_type_pair = update(
      conditional_model,
      dispformula = ~type * pair),
    zero_inflated = update(
      conditional_model,
      ziformula = ~1),
    zero_inflated_type = update(
      conditional_model,
      ziformula = ~type))
  
  candidate_comparison <- bind_rows(
    lapply(names(candidate_models), function(x) {
      tibble(
        model = x,
        df = attr(
          logLik(candidate_models[[x]]),
          "df"),
        AIC = AIC(candidate_models[[x]]))
    })) %>%
    arrange(AIC)
  
  cat("\n### Dispersion and zero-inflation comparison ###\n")
  print(candidate_comparison)
  
  # Retain the best-supported nuisance variance structure,
  # provided the model converges and residual diagnostics pass.
  
  juvenile_final <- candidate_models$dispersion_type
  
  
  #### 4.5 Final-model diagnostics ####
  
  final_residuals <- DHARMa::simulateResiduals(
    juvenile_final,
    n = 1000,
    plot = FALSE)
  
  diagnostic_tests <- list(
    uniformity = DHARMa::testUniformity(
      final_residuals),
    dispersion = DHARMa::testDispersion(
      final_residuals),
    outliers = DHARMa::testOutliers(
      final_residuals,
      type = "bootstrap"),
    zero_inflation = DHARMa::testZeroInflation(
      final_residuals))
  
  cat("\n### Final life-stage model diagnostics ###\n")
  print(diagnostic_tests$uniformity)
  print(diagnostic_tests$dispersion)
  print(diagnostic_tests$outliers)
  print(diagnostic_tests$zero_inflation)
  
  
  #### 4.6 Omnibus tests ####
  
  stage_tests <- emmeans::joint_tests(
    juvenile_final)
  
  cat("\n### Omnibus fixed-effect tests ###\n")
  print(stage_tests)
  
  
  #### 4.7 Predicted counts and derived proportions ####
  
  stage_emm <- emmeans::emmeans(
    juvenile_final,
    ~life_stage | type * pair,
    type = "response")
  
  stage_predictions <- summary(
    stage_emm,
    infer = c(TRUE, TRUE),
    adjust = "none") %>%
    as_tibble()
  
  stage_proportions <- stage_predictions %>%
    select(
      type, pair, life_stage, response) %>%
    pivot_wider(
      names_from = life_stage,
      values_from = response) %>%
    mutate(
      total = juvenile + adult,
      juvenile_prop = juvenile / total,
      adult_prop = adult / total)
  
  cat("\n### Predicted stage counts ###\n")
  print(stage_predictions)
  
  cat("\n### Derived life-stage proportions ###\n")
  print(stage_proportions)
  
  
  #### 4.8 Artificial versus natural contrasts ####
  
  reef_emm <- emmeans::emmeans(
    juvenile_final,
    ~type | life_stage * pair,
    type = "response")
  
  stage_contrasts <- emmeans::contrast(
    reef_emm,
    method = list(
      "Artificial / Natural" = c(-1, 1)),
    by = c("life_stage", "pair")) %>%
    summary(
      infer = c(TRUE, TRUE),
      adjust = "none") %>%
    as_tibble()
  
  cat("\n### Artificial versus natural counts ###\n")
  print(stage_contrasts)
  
  
  #### 4.9 Adult versus juvenile contrasts ####
  
  life_stage_contrasts <- emmeans::contrast(
    stage_emm,
    method = list(
      "Adult / Juvenile" = c(-1, 1)),
    by = c("type", "pair")) %>%
    summary(
      infer = c(TRUE, TRUE),
      adjust = "none") %>%
    as_tibble()
  
  cat("\n### Adult versus juvenile counts ###\n")
  print(life_stage_contrasts)
  
  
  #### 4.10 Pair-specific life-stage composition ####
  
  stage_grid_pair <- emmeans::emmeans(
    juvenile_final,
    ~life_stage * type | pair)
  
  proportion_contrasts_pair <- emmeans::contrast(
    stage_grid_pair,
    interaction = c(
      "revpairwise",
      "revpairwise"),
    by = "pair") %>%
    summary(
      infer = c(TRUE, TRUE),
      type = "response",
      adjust = "none") %>%
    as_tibble()
  
  cat("\n### Reef differences in life-stage balance by pair ###\n")
  print(proportion_contrasts_pair)
  
  
  #### 4.11 Overall life-stage composition ####
  
  stage_grid_overall <- emmeans::emmeans(
    juvenile_final,
    ~life_stage * type)
  
  proportion_contrast_overall <- emmeans::contrast(
    stage_grid_overall,
    interaction = c(
      "revpairwise",
      "revpairwise")) %>%
    summary(
      infer = c(TRUE, TRUE),
      type = "response",
      adjust = "none") %>%
    as_tibble()
  
  cat("\n### Overall reef difference in life-stage balance ###\n")
  print(proportion_contrast_overall)
  
  
  #### 4.12 Return outputs ####
  
  list(
    data = juvenile_data,
    model = juvenile_final,
    distribution_comparison =
      distribution_comparison,
    candidate_comparison =
      candidate_comparison,
    diagnostics = diagnostic_tests,
    tests = stage_tests,
    predictions = stage_predictions,
    proportions = stage_proportions,
    reef_contrasts = stage_contrasts,
    stage_contrasts =
      life_stage_contrasts,
    pair_composition =
      proportion_contrasts_pair,
    overall_composition =
      proportion_contrast_overall)
})


#### 4.13 Save outputs ####

juvenile_data <- juvenile_results$data
juvenile_final <- juvenile_results$model

juvenile_distribution_comparison <-
  juvenile_results$distribution_comparison

juvenile_candidate_comparison <-
  juvenile_results$candidate_comparison

juvenile_diagnostics <-
  juvenile_results$diagnostics

juvenile_tests <-
  juvenile_results$tests

stage_predictions <-
  juvenile_results$predictions

stage_proportions <-
  juvenile_results$proportions

stage_contrasts <-
  juvenile_results$reef_contrasts

life_stage_contrasts <-
  juvenile_results$stage_contrasts

proportion_contrasts_pair <-
  juvenile_results$pair_composition

proportion_contrast_overall <-
  juvenile_results$overall_composition

rm(juvenile_results)

summary(juvenile_final)
