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