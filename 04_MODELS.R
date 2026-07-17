## 04_MODELS.R ####

library(dplyr)
library(tidyr)
library(glmmTMB)
library(lme4)
library(lmerTest)
library(mgcv)
library(DHARMa)
library(performance)
library(emmeans)
library(vegan)


## 1. Settings and helpers ####

target_taxa <- c(
  "Damsels - Regal Demoiselle",
  "Rabbit - Java",
  "Snapper - Russells",
  "Snapper - Mangrove"
)

# A GAMM must improve temporal holdout log-RMSE by at least 2% to replace
# the simpler GLMM, unless GLMM diagnostics clearly fail.
gamm_improvement_threshold <- 2

log_rmse <- function(observed, predicted) {
  sqrt(mean((log1p(observed) - log1p(predicted))^2))
}

check_glmm <- function(model, data, time = NULL, group = NULL) {
  cat("\n")
  print(performance::check_collinearity(model))
  print(performance::check_singularity(model))

  residuals <- DHARMa::simulateResiduals(model, plot = FALSE)
  plot(residuals)
  print(DHARMa::testDispersion(residuals))
  print(DHARMa::testZeroInflation(residuals))
  print(DHARMa::testOutliers(residuals))

  if (!is.null(time) && !is.null(group)) {
    grouped_residuals <- DHARMa::recalculateResiduals(
      residuals,
      group = group
    )

    grouped_time <- tapply(
      time,
      group,
      function(x) x[1]
    )

    print(
      DHARMa::testTemporalAutocorrelation(
        grouped_residuals,
        time = as.numeric(grouped_time)
      )
    )
  }

  invisible(residuals)
}



## 2. Feeding-guild abundance trajectories ####

abundance_results <- local({
  
  #### 2.1 Analysis data ####
  
  abundance_data <- fish_counts %>%
    mutate(
      row_id = row_number(),
      trajectory = interaction(
        pair, type, feeding_guild,
        drop = TRUE
      ),
      time_series = interaction(
        site, feeding_guild,
        drop = TRUE
      )
    )
  
  
  #### 2.2 Negative-binomial family ####
  
  nb1 <- glmmTMB(
    Count ~ feeding_guild * type * pair * date_s +
      (1 | survey_id),
    family = nbinom1,
    data = abundance_data
  )
  
  nb2 <- glmmTMB(
    Count ~ feeding_guild * type * pair * date_s +
      (1 | survey_id),
    family = nbinom2,
    data = abundance_data
  )
  
  cat("\n### Negative-binomial family comparison ###\n")
  print(AIC(nb1, nb2))
  
  res_nb1 <- DHARMa::simulateResiduals(
    nb1,
    plot = FALSE
  )
  
  res_nb2 <- DHARMa::simulateResiduals(
    nb2,
    plot = FALSE
  )
  
  cat("\n### nbinom1 diagnostics ###\n")
  print(DHARMa::testDispersion(res_nb1))
  print(DHARMa::testZeroInflation(res_nb1))
  print(DHARMa::testOutliers(
    res_nb1,
    type = "bootstrap"
  ))
  
  cat("\n### nbinom2 diagnostics ###\n")
  print(DHARMa::testDispersion(res_nb2))
  print(DHARMa::testZeroInflation(res_nb2))
  print(DHARMa::testOutliers(
    res_nb2,
    type = "bootstrap"
  ))
  
  
  #### 2.3 Zero-inflation structure ####
  
  zi_constant <- update(
    nb2,
    ziformula = ~1
  )
  
  zi_guild <- update(
    nb2,
    ziformula = ~feeding_guild
  )
  
  cat("\n### Zero-inflation comparison ###\n")
  print(AIC(
    nb2,
    zi_constant,
    zi_guild
  ))
  
  res_zi_constant <- DHARMa::simulateResiduals(
    zi_constant,
    plot = FALSE
  )
  
  res_zi_guild <- DHARMa::simulateResiduals(
    zi_guild,
    plot = FALSE
  )
  
  cat("\n### Constant zero-inflation diagnostics ###\n")
  print(DHARMa::testDispersion(res_zi_constant))
  print(DHARMa::testZeroInflation(res_zi_constant))
  print(DHARMa::testOutliers(
    res_zi_constant,
    type = "bootstrap"
  ))
  
  cat("\n### Guild-specific zero-inflation diagnostics ###\n")
  print(DHARMa::testDispersion(res_zi_guild))
  print(DHARMa::testZeroInflation(res_zi_guild))
  print(DHARMa::testOutliers(
    res_zi_guild,
    type = "bootstrap"
  ))
  
  
  #### 2.4 Fixed-effect structure ####
  
  original <- glmmTMB(
    Count ~ type * pair * date_s +
      feeding_guild +
      (1 | survey_id),
    ziformula = ~1,
    family = nbinom2,
    data = abundance_data
  )
  
  reduced <- glmmTMB(
    Count ~
      (feeding_guild + type + pair + date_s)^3 +
      (1 | survey_id),
    ziformula = ~1,
    family = nbinom2,
    data = abundance_data
  )
  
  full <- zi_constant
  
  cat("\n### Original versus full structure ###\n")
  print(AIC(original, full))
  print(anova(original, full))
  
  cat("\n### Reduced versus full interaction structure ###\n")
  print(AIC(reduced, full))
  print(anova(reduced, full))
  
  
  #### 2.5 GLMM versus GAMM temporal holdout ####
  
  holdout_surveys <- abundance_data %>%
    distinct(site, survey_id, Date) %>%
    group_by(site) %>%
    arrange(Date, .by_group = TRUE) %>%
    mutate(
      test = row_number() > floor(0.8 * n())
    ) %>%
    ungroup() %>%
    filter(test) %>%
    pull(survey_id)
  
  train <- abundance_data %>%
    filter(!survey_id %in% holdout_surveys)
  
  test <- abundance_data %>%
    filter(survey_id %in% holdout_surveys)
  
  glmm_train <- glmmTMB(
    Count ~ feeding_guild * type * pair * date_s +
      (1 | survey_id),
    ziformula = ~1,
    family = nbinom2,
    data = train
  )
  
  gamm_train <- bam(
    Count ~ trajectory +
      s(date_s, by = trajectory, k = 4) +
      s(survey_id, bs = "re"),
    family = nb(),
    method = "fREML",
    discrete = TRUE,
    data = train
  )
  
  glmm_prediction <- predict(
    glmm_train,
    newdata = test,
    type = "response",
    re.form = NA,
    allow.new.levels = TRUE
  )
  
  test_gamm <- test %>%
    mutate(
      survey_id = factor(
        levels(train$survey_id)[1],
        levels = levels(train$survey_id)
      )
    )
  
  gamm_prediction <- predict(
    gamm_train,
    newdata = test_gamm,
    type = "response",
    exclude = "s(survey_id)"
  )
  
  model_comparison <- tibble(
    model = c("GLMM", "GAMM"),
    test_log_rmse = c(
      log_rmse(test$Count, glmm_prediction),
      log_rmse(test$Count, gamm_prediction)
    )
  ) %>%
    mutate(
      difference_from_glmm_pct =
        100 *
        (
          test_log_rmse -
            test_log_rmse[model == "GLMM"]
        ) /
        test_log_rmse[model == "GLMM"]
    )
  
  cat("\n### GLMM versus GAMM temporal holdout ###\n")
  print(model_comparison)
  
  
  #### 2.6 Final-model diagnostics ####
  
  abundance_final <- full
  
  final_residuals <- DHARMa::simulateResiduals(
    abundance_final,
    plot = FALSE
  )
  
  plot(final_residuals)
  
  cat("\n### Final abundance diagnostics ###\n")
  print(DHARMa::testDispersion(final_residuals))
  print(DHARMa::testZeroInflation(final_residuals))
  print(DHARMa::testOutliers(
    final_residuals,
    type = "bootstrap"
  ))
  
  
  #### 2.7 Temporal autocorrelation ####
  
  abundance_temporal_ac <- abundance_data %>%
    group_by(time_series) %>%
    group_modify(~{
      
      test <- DHARMa::testTemporalAutocorrelation(
        final_residuals$scaledResiduals[
          .x$row_id
        ],
        time = .x$date_num
      )
      
      tibble(
        p_value = test$p.value
      )
    }) %>%
    ungroup() %>%
    mutate(
      p_adjusted = p.adjust(
        p_value,
        method = "BH"
      )
    )
  
  cat("\n### Temporal autocorrelation by site and feeding guild ###\n")
  print(abundance_temporal_ac, n = Inf)
  
  abundance_data %>%
    mutate(
      residual = final_residuals$scaledResiduals
    ) %>%
    filter(
      time_series %in%
        abundance_temporal_ac$time_series[
          abundance_temporal_ac$p_adjusted < 0.05
        ]
    ) %>%
    ggplot(aes(Date, residual)) +
    geom_hline(
      yintercept = 0.5,
      linetype = 2
    ) +
    geom_point() +
    geom_line() +
    facet_wrap(
      ~time_series,
      scales = "free_x"
    ) +
    labs(
      x = "Survey date",
      y = "DHARMa residual"
    ) +
    theme_clean %>%
    print()
  
  
  list(
    data = abundance_data,
    model = abundance_final,
    autocorrelation = abundance_temporal_ac
  )
})


abundance_data <- abundance_results$data
abundance_final <- abundance_results$model
abundance_temporal_ac <- abundance_results$autocorrelation

rm(abundance_results)


## 2.8 Total fish abundance ####

total_abundance_results <- local({
  
  #### 2.8.1 Survey-level totals ####
  
  total_data <- fish_counts %>%
    group_by(
      survey_id, site, pair, type, period,
      Date, date_num, date_s
    ) %>%
    summarise(
      Count = sum(Count),
      .groups = "drop"
    ) %>%
    mutate(
      row_id = row_number()
    )
  
  
  #### 2.8.2 Negative-binomial family ####
  
  nb1 <- glmmTMB(
    Count ~ type * pair * date_s,
    family = nbinom1,
    data = total_data
  )
  
  nb2 <- glmmTMB(
    Count ~ type * pair * date_s,
    family = nbinom2,
    data = total_data
  )
  
  cat("\n### Total abundance family comparison ###\n")
  print(AIC(nb1, nb2))
  
  res_nb1 <- DHARMa::simulateResiduals(nb1, plot = FALSE)
  res_nb2 <- DHARMa::simulateResiduals(nb2, plot = FALSE)
  
  cat("\n### Total abundance nbinom1 diagnostics ###\n")
  print(DHARMa::testDispersion(res_nb1))
  print(DHARMa::testZeroInflation(res_nb1))
  print(DHARMa::testOutliers(
    res_nb1,
    type = "bootstrap"
  ))
  
  cat("\n### Total abundance nbinom2 diagnostics ###\n")
  print(DHARMa::testDispersion(res_nb2))
  print(DHARMa::testZeroInflation(res_nb2))
  print(DHARMa::testOutliers(
    res_nb2,
    type = "bootstrap"
  ))
  
  family_models <- list(
    nbinom1 = nb1,
    nbinom2 = nb2
  )
  
  family_aic <- AIC(nb1, nb2)
  
  base_model <- family_models[[
    rownames(family_aic)[which.min(family_aic$AIC)]
  ]]
  
  
  #### 2.8.3 Fixed-effect structure ####
  
  full <- nb1
  
  reduced <- glmmTMB(
    Count ~ (type + pair + date_s)^2,
    family = nbinom1,
    data = total_data
  )
  
  cat("\n### Total abundance fixed-effect comparison ###\n")
  print(AIC(reduced, full))
  print(anova(reduced, full))
  
  
  #### 2.8.4 Final diagnostics ####
  
  total_final <- full
  
  final_residuals <- DHARMa::simulateResiduals(
    total_final,
    plot = FALSE
  )
  
  plot(final_residuals)
  
  cat("\n### Final total abundance diagnostics ###\n")
  print(DHARMa::testDispersion(final_residuals))
  print(DHARMa::testZeroInflation(final_residuals))
  print(DHARMa::testOutliers(
    final_residuals,
    type = "bootstrap"
  ))
  
  #### 2.8.5 Temporal autocorrelation ####
  
  total_ac <- total_data %>%
    group_by(site) %>%
    arrange(Date, .by_group = TRUE) %>%
    group_modify(~{
      
      test <- DHARMa::testTemporalAutocorrelation(
        final_residuals$scaledResiduals[.x$row_id],
        time = .x$date_num
      )
      
      tibble(
        p_value = test$p.value
      )
    }) %>%
    ungroup() %>%
    mutate(
      p_adjusted = p.adjust(
        p_value,
        method = "BH"
      )
    )
  
  cat("\n### Total abundance temporal autocorrelation ###\n")
  print(total_ac, n = Inf)
  
  
  list(
    data = total_data,
    model = total_final,
    autocorrelation = total_ac
  )
})


total_abundance_data <- total_abundance_results$data
total_abundance_final <- total_abundance_results$model
total_abundance_ac <- total_abundance_results$autocorrelation

rm(total_abundance_results)

### final abundance models #### 
# with f groups 
# abundance_final <- glmmTMB(
#   Count ~ feeding_guild * type * pair * date_s +
#     (1 | survey_id),
#   ziformula = ~1,
#   family = nbinom2,
#   data = abundance_data
# )

# total counts
# total_abundance_final <- glmmTMB(
#   Count ~ type * pair * date_s,
#   family = nbinom1,
#   data = total_abundance_data
# )




## 3. Shannon diversity trajectories ####

diversity_results <- local({
  
  #### 3.1 Survey-level Shannon diversity ####
  
  community_wide <- fish_long %>%
    select(
      survey_id, site, pair, type, period,
      Date, date_num, date_s,
      Species, Count
    ) %>%
    pivot_wider(
      names_from = Species,
      values_from = Count,
      values_fill = 0
    )
  
  species_cols <- setdiff(
    names(community_wide),
    c(
      "survey_id", "site", "pair", "type",
      "period", "Date", "date_num", "date_s"
    )
  )
  
  diversity_data <- community_wide %>%
    mutate(
      shannon = vegan::diversity(
        as.matrix(pick(all_of(species_cols))),
        index = "shannon"
      )
    ) %>%
    select(
      survey_id, site, pair, type, period,
      Date, date_num, date_s, shannon
    )
  
  
  #### 3.2 Fixed-effect structure ####
  
  no_step <- lm(
    shannon ~ type * date_s + pair,
    data = diversity_data
  )
  
  common_step <- lm(
    shannon ~ type * date_s + period + pair,
    data = diversity_data
  )
  
  type_step <- lm(
    shannon ~ type * (date_s + period) + pair,
    data = diversity_data
  )
  
  cat("\n### Shannon fixed-effect comparison ###\n")
  print(AIC(no_step, common_step, type_step))
  
  cat("\n### Nested model tests ###\n")
  print(anova(no_step, common_step, type_step))
  
  
  #### 3.3 Variance structure ####
  
  variance_models <- list(
    constant = glmmTMB(
      shannon ~ type * (date_s + period) + pair,
      family = gaussian,
      data = diversity_data
    ),
    
    type = glmmTMB(
      shannon ~ type * (date_s + period) + pair,
      dispformula = ~ type,
      family = gaussian,
      data = diversity_data
    ),
    
    period = glmmTMB(
      shannon ~ type * (date_s + period) + pair,
      dispformula = ~ period,
      family = gaussian,
      data = diversity_data
    ),
    
    type_period = glmmTMB(
      shannon ~ type * (date_s + period) + pair,
      dispformula = ~ type * period,
      family = gaussian,
      data = diversity_data
    )
  )
  
  variance_comparison <- bind_rows(
    lapply(names(variance_models), function(x) {
      tibble(
        model = x,
        df = attr(logLik(variance_models[[x]]), "df"),
        AIC = AIC(variance_models[[x]])
      )
    })
  )
  
  cat("\n### Shannon variance comparison ###\n")
  print(variance_comparison)
  
  shannon_final <- variance_models[[
    variance_comparison$model[
      which.min(variance_comparison$AIC)
    ]
  ]]
  
  
  #### 3.4 Final-model diagnostics ####
  
  shannon_residuals <- DHARMa::simulateResiduals(
    shannon_final,
    plot = FALSE
  )
  
  plot(shannon_residuals)
  
  cat("\n### Final Shannon diagnostics ###\n")
  print(DHARMa::testUniformity(shannon_residuals))
  print(DHARMa::testDispersion(shannon_residuals))
  print(
    DHARMa::testOutliers(
      shannon_residuals,
      type = "bootstrap"
    )
  )
  
  
  #### 3.5 Temporal autocorrelation ####
  
  shannon_ac <- diversity_data %>%
    mutate(
      residual = shannon_residuals$scaledResiduals
    ) %>%
    group_by(site) %>%
    arrange(Date, .by_group = TRUE) %>%
    group_modify(~{
      
      test <- DHARMa::testTemporalAutocorrelation(
        .x$residual,
        time = .x$date_num
      )
      
      tibble(
        p_value = test$p.value
      )
    }) %>%
    ungroup() %>%
    mutate(
      p_adjusted = p.adjust(
        p_value,
        method = "BH"
      )
    )
  
  cat("\n### Shannon temporal autocorrelation ###\n")
  print(shannon_ac, n = Inf)
  
  
  list(
    data = diversity_data,
    model = shannon_final,
    autocorrelation = shannon_ac
  )
})


diversity_data <- diversity_results$data
shannon_final <- diversity_results$model
shannon_ac <- diversity_results$autocorrelation

rm(diversity_results)

###  final diversity model #### 
# shannon_final <- glmmTMB(
#   shannon ~ type * (date_s + period) + pair,
#   dispformula = ~ type,
#   family = gaussian,
#   data = diversity_data
# )

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



