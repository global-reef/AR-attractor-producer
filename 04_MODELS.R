### 04_MODELS.R ###

#### 1. Helpers ####

log_rmse <- function(observed, predicted) sqrt(mean((log1p(observed) - log1p(predicted))^2))

run_dharma <- function(model, n = 1000, seed = 42, plot = TRUE) {
  set.seed(seed)
  residuals <- DHARMa::simulateResiduals(model, n = n, plot = FALSE)
  if (plot) plot(residuals)
  
  list(
    residuals = residuals,
    uniformity = DHARMa::testUniformity(residuals),
    dispersion = DHARMa::testDispersion(residuals),
    zero_inflation = DHARMa::testZeroInflation(residuals),
    outliers = DHARMa::testOutliers(residuals, type = "bootstrap")
  )
}

temporal_ac <- function(data, residuals, group_var, time_var = "date_num", adjust = "holm") {
  data %>%
    mutate(.residual = residuals$scaledResiduals) %>%
    group_by(across(all_of(group_var))) %>%
    arrange(.data[[time_var]], .by_group = TRUE) %>%
    group_modify(~ {
      test <- DHARMa::testTemporalAutocorrelation(.x$.residual, time = .x[[time_var]])
      tibble(p_value = test$p.value)
    }) %>%
    ungroup() %>%
    mutate(p_adjusted = p.adjust(p_value, method = adjust))
}


#### 2. Feeding-guild abundance model ####

abundance_data <- fish_counts %>%
  mutate(
    row_id = row_number(),
    trajectory = interaction(pair, type, feeding_guild, drop = TRUE),
    time_series = interaction(site, feeding_guild, drop = TRUE)
  )

abundance_nb1 <- glmmTMB(
  Count ~ feeding_guild * type * pair * date_s + (1 | survey_id),
  ziformula = ~feeding_guild, family = nbinom1, data = abundance_data
)

abundance_nb2 <- update(abundance_nb1, family = nbinom2)
abundance_family_comparison <- AIC(abundance_nb1, abundance_nb2)

abundance_final <- abundance_nb2
abundance_diagnostics <- run_dharma(abundance_final)
abundance_temporal_ac <- temporal_ac(
  abundance_data, abundance_diagnostics$residuals,
  group_var = "time_series", adjust = "holm"
)

cat("\n### Abundance family comparison ###\n")
print(abundance_family_comparison)
cat("\n### Final abundance diagnostics ###\n")
print(abundance_diagnostics[c("uniformity", "dispersion", "zero_inflation", "outliers")])
cat("\n### Abundance temporal autocorrelation ###\n")
print(abundance_temporal_ac, n = Inf)


#### 2.1 GLMM-GAMM temporal holdout sensitivity ####

holdout_surveys <- abundance_data %>%
  distinct(site, survey_id, Date) %>%
  group_by(site) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(test = row_number() > floor(0.8 * n())) %>%
  ungroup() %>%
  filter(test) %>%
  pull(survey_id)

abundance_train <- filter(abundance_data, !survey_id %in% holdout_surveys)
abundance_test <- filter(abundance_data, survey_id %in% holdout_surveys)

abundance_glmm_train <- glmmTMB(
  Count ~ feeding_guild * type * pair * date_s + (1 | survey_id),
  ziformula = ~feeding_guild, family = nbinom2, data = abundance_train
)

abundance_gamm_train <- mgcv::bam(
  Count ~ trajectory + s(date_s, by = trajectory, k = 4) + s(survey_id, bs = "re"),
  family = mgcv::nb(), method = "fREML", discrete = TRUE, data = abundance_train
)

glmm_prediction <- predict(
  abundance_glmm_train, newdata = abundance_test, type = "response",
  re.form = NA, allow.new.levels = TRUE
)

gamm_test_data <- abundance_test %>%
  mutate(survey_id = factor(levels(abundance_train$survey_id)[1], levels = levels(abundance_train$survey_id)))

gamm_prediction <- predict(
  abundance_gamm_train, newdata = gamm_test_data,
  type = "response", exclude = "s(survey_id)"
)

abundance_holdout_comparison <- tibble(
  model = c("GLMM", "GAMM"),
  test_log_rmse = c(
    log_rmse(abundance_test$Count, glmm_prediction),
    log_rmse(abundance_test$Count, gamm_prediction)
  )
) %>%
  mutate(difference_from_glmm_pct = 100 * (test_log_rmse - first(test_log_rmse)) / first(test_log_rmse))

cat("\n### GLMM-GAMM temporal holdout ###\n")
print(abundance_holdout_comparison)


#### 3. Total abundance model ####

total_abundance_data <- fish_counts %>%
  group_by(survey_id, site, pair, type, period, Date, date_num, date_s) %>%
  summarise(Count = sum(Count), .groups = "drop")

total_nb1 <- glmmTMB(Count ~ type * pair * date_s, family = nbinom1, data = total_abundance_data)
total_nb2 <- update(total_nb1, family = nbinom2)
total_family_comparison <- AIC(total_nb1, total_nb2)

total_abundance_final <- total_nb1
total_abundance_diagnostics <- run_dharma(total_abundance_final)
total_abundance_ac <- temporal_ac(
  total_abundance_data, total_abundance_diagnostics$residuals,
  group_var = "site", adjust = "holm"
)

cat("\n### Total abundance family comparison ###\n")
print(total_family_comparison)
cat("\n### Final total abundance diagnostics ###\n")
print(total_abundance_diagnostics[c("uniformity", "dispersion", "zero_inflation", "outliers")])


#### 4. Shannon diversity model ####

community_wide <- fish_long %>%
  select(survey_id, site, pair, type, period, Date, date_num, date_s, Species, Count) %>%
  pivot_wider(names_from = Species, values_from = Count, values_fill = 0)

diversity_metadata <- c("survey_id", "site", "pair", "type", "period", "Date", "date_num", "date_s")
diversity_species <- setdiff(names(community_wide), diversity_metadata)

diversity_data <- community_wide %>%
  mutate(shannon = vegan::diversity(as.matrix(pick(all_of(diversity_species))), index = "shannon")) %>%
  select(all_of(diversity_metadata), shannon)

shannon_constant <- glmmTMB(
  shannon ~ type * (date_s + period) + pair,
  family = gaussian, data = diversity_data
)

shannon_type_dispersion <- update(shannon_constant, dispformula = ~type)
shannon_variance_comparison <- AIC(shannon_constant, shannon_type_dispersion)

shannon_final <- shannon_type_dispersion
shannon_diagnostics <- run_dharma(shannon_final)
shannon_ac <- temporal_ac(
  diversity_data, shannon_diagnostics$residuals,
  group_var = "site", adjust = "holm"
)

cat("\n### Shannon variance comparison ###\n")
print(shannon_variance_comparison)
cat("\n### Final Shannon diagnostics ###\n")
print(shannon_diagnostics[c("uniformity", "dispersion", "outliers")])


#### 5. Community life-stage model ####

survey_frame <- fish_size %>% distinct(survey_id, site, pair, type, Date)
species_frame <- fish_size %>% distinct(Species)

stage_counts <- fish_size %>%
  group_by(survey_id, site, pair, type, Date, Species, life_stage) %>%
  summarise(stage_count = sum(stage_Count), .groups = "drop")

juvenile_data <- crossing(
  survey_frame,
  species_frame,
  life_stage = factor(life_stage_levels, levels = life_stage_levels)
) %>%
  left_join(
    stage_counts,
    by = c("survey_id", "site", "pair", "type", "Date", "Species", "life_stage")
  ) %>%
  mutate(
    stage_count = as.integer(replace_na(stage_count, 0)),
    life_stage = factor(life_stage, levels = life_stage_levels),
    type = factor(type, levels = reef_type_levels),
    pair = factor(pair, levels = pair_levels),
    Species = factor(Species),
    survey_id = factor(survey_id)
  )

juvenile_poisson <- glmmTMB(
  stage_count ~ life_stage * type * pair + (1 | Species) + (1 | survey_id),
  family = poisson, data = juvenile_data
)

juvenile_nb1 <- update(juvenile_poisson, family = nbinom1)
juvenile_nb2 <- update(juvenile_poisson, family = nbinom2)
juvenile_distribution_comparison <- AIC(juvenile_poisson, juvenile_nb1, juvenile_nb2)

juvenile_constant <- juvenile_nb2
juvenile_type_dispersion <- update(juvenile_constant, dispformula = ~type)
juvenile_candidate_comparison <- AIC(juvenile_constant, juvenile_type_dispersion)

juvenile_final <- juvenile_type_dispersion
juvenile_diagnostics <- run_dharma(juvenile_final)

cat("\n### Life-stage distribution comparison ###\n")
print(juvenile_distribution_comparison)
cat("\n### Life-stage dispersion comparison ###\n")
print(juvenile_candidate_comparison)
cat("\n### Final life-stage diagnostics ###\n")
print(juvenile_diagnostics[c("uniformity", "dispersion", "zero_inflation", "outliers")])


#### 6. Species-specific life-stage models ####

focal_taxa <- c(
  "Damsels - Regal Demoiselle",
  "Damsels - Alexanders",
  "Snapper - Russells",
  "Snapper - Brownstripe",
  "Snapper - Mangrove"
)

focal_observed <- fish_size %>%
  filter(Species %in% focal_taxa) %>%
  group_by(survey_id, site, pair, type, Date, Species, Sci_Name, life_stage) %>%
  summarise(stage_count = sum(stage_Count), .groups = "drop")

focal_species_key <- fish_size %>%
  filter(Species %in% focal_taxa) %>%
  distinct(Species, Sci_Name)

focal_stage <- survey_frame %>%
  crossing(focal_species_key, life_stage = factor(life_stage_levels, levels = life_stage_levels)) %>%
  left_join(
    focal_observed,
    by = c("survey_id", "site", "pair", "type", "Date", "Species", "Sci_Name", "life_stage")
  ) %>%
  mutate(
    stage_count = as.integer(replace_na(stage_count, 0)),
    Species = factor(Species, levels = focal_taxa),
    pair = factor(pair, levels = pair_levels),
    type = factor(type, levels = reef_type_levels),
    life_stage = factor(life_stage, levels = c("adult", "juvenile")),
    survey_id = factor(survey_id)
  )

species_data <- focal_stage %>%
  droplevels() %>%
  group_split(Species) %>%
  set_names(map_chr(., ~ as.character(unique(.x$Species))))

fit_species_model <- function(data) {
  glmmTMB(
    stage_count ~ life_stage * type * pair + (1 | survey_id),
    family = nbinom2, data = data,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
  )
}

species_models <- map(species_data, fit_species_model)

species_model_summary <- imap_dfr(species_models, \(model, species) {
  set.seed(42)
  residuals <- DHARMa::simulateResiduals(model, n = 1000, plot = FALSE)
  survey_variance <- VarCorr(model)$cond$survey_id
  
  tibble(
    Species = species,
    convergence_code = model$fit$convergence,
    convergence_message = model$fit$message,
    positive_definite_hessian = model$sdr$pdHess,
    singular = performance::check_singularity(model),
    survey_sd = attr(survey_variance, "stddev"),
    dispersion = sigma(model),
    uniformity_p = DHARMa::testUniformity(residuals, plot = FALSE)$p.value,
    dispersion_p = DHARMa::testDispersion(residuals, plot = FALSE)$p.value,
    zero_inflation_p = DHARMa::testZeroInflation(residuals, plot = FALSE)$p.value,
    outlier_p = DHARMa::testOutliers(residuals, plot = FALSE)$p.value
  )
})

cat("\n### Species-model checks ###\n")
print(species_model_summary, n = Inf)


#### 7. Save final model summaries ####

save_model_summary(abundance_final, "abundance_final")
save_model_summary(total_abundance_final, "total_abundance_final")
save_model_summary(shannon_final, "shannon_final")
save_model_summary(juvenile_final, "juvenile_final")

iwalk(species_models, ~ save_model_summary(.x, paste0("species_", make.names(.y))))
message("Model fitting complete.")
