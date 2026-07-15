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



