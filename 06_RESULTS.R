#### 06_RESULTS.R ####


## 1.1 Extract abundance results ####

extract_abundance_answers <- function(
    model,
    data,
    guild_adjust = "holm",
    n_sim = 5000,
    seed = 42
) {
  
  #### 1. Overall artificial versus natural trajectories ####
  
  overall_emm <- emmeans::emtrends(
    model,
    ~ type,
    var = "date_s",
    weights = "equal"
  )
  
  overall_slopes <- summary(
    overall_emm,
    infer = c(TRUE, TRUE),
    adjust = "none"
  ) %>%
    as_tibble() %>%
    rename(
      slope = date_s.trend,
      p_value = p.value
    ) %>%
    mutate(
      trajectory = case_when(
        p_value < 0.05 & slope > 0 ~ "increase",
        p_value < 0.05 & slope < 0 ~ "decrease",
        TRUE ~ "no significant change"
      )
    )
  
  overall_contrast <- emmeans::contrast(
    overall_emm,
    method = list(
      "Artificial - Natural" = c(1, -1)
    )
  ) %>%
    summary(
      infer = c(TRUE, TRUE),
      adjust = "none"
    ) %>%
    as_tibble() %>%
    rename(
      slope_difference = estimate,
      p_value = p.value
    )
  
  
  #### 2. Artificial and natural trajectories within each pair ####
  
  pair_emm <- emmeans::emtrends(
    model,
    ~ type | pair,
    var = "date_s",
    weights = "equal"
  )
  
  pair_slopes <- summary(
    update(pair_emm, by = NULL),
    infer = c(TRUE, TRUE),
    adjust = "none"
  ) %>%
    as_tibble() %>%
    rename(
      slope = date_s.trend,
      p_value = p.value
    ) %>%
    mutate(
      trajectory = case_when(
        p_value < 0.05 & slope > 0 ~ "increase",
        p_value < 0.05 & slope < 0 ~ "decrease",
        TRUE ~ "no significant change"
      )
    )
  
  pair_contrasts <- emmeans::contrast(
    pair_emm,
    method = list(
      "Artificial - Natural" = c(1, -1)
    ),
    by = "pair"
  ) %>%
    update(by = NULL) %>%
    summary(
      infer = c(TRUE, TRUE),
      adjust = "none"
    ) %>%
    as_tibble() %>%
    rename(
      slope_difference = estimate,
      p_value = p.value
    )
  
  #### 3. Combined artificial + natural abundance within each pair ####
  
  paired_total_change <- bind_rows(
    lapply(
      levels(data$pair),
      function(pair_name) {
        
        pair_data <- data %>%
          filter(pair == pair_name)
        
        time_values <- c(
          min(pair_data$date_s),
          max(pair_data$date_s)
        )
        
        total_emm <- emmeans::emmeans(
          model,
          ~ date_s,
          at = list(
            pair = pair_name,
            date_s = time_values
          ),
          weights = "equal",
          component = "response"
        ) %>%
          emmeans::regrid(
            transform = "response"
          )
        
        emmeans::contrast(
          total_emm,
          method = list(
            "End - Start" = c(-8, 8)
          )
        ) %>%
          summary(
            infer = c(TRUE, TRUE),
            adjust = "none"
          ) %>%
          as_tibble() %>%
          mutate(
            pair = pair_name
          )
      }
    )
  ) %>%
    mutate(
      total_trajectory = case_when(
        p.value < 0.05 & estimate > 0 ~ "increase",
        p.value < 0.05 & estimate < 0 ~ "decrease",
        TRUE ~ "no significant change"
      )
    )
  #### 4. Feeding-guild decomposition ####
  
  guild_emm <- emmeans::emtrends(
    model,
    ~ type | pair * feeding_guild,
    var = "date_s"
  )
  
  guild_slopes <- summary(
    update(guild_emm, by = NULL),
    infer = c(TRUE, TRUE),
    adjust = guild_adjust
  ) %>%
    as_tibble() %>%
    rename(
      slope = date_s.trend,
      p_adjusted = p.value
    ) %>%
    mutate(
      trajectory = case_when(
        p_adjusted < 0.05 & slope > 0 ~ "increase",
        p_adjusted < 0.05 & slope < 0 ~ "decrease",
        TRUE ~ "no significant change"
      )
    )
  
  guild_contrasts <- emmeans::contrast(
    guild_emm,
    method = list(
      "Artificial - Natural" = c(1, -1)
    ),
    by = c("pair", "feeding_guild")
  ) %>%
    update(by = NULL) %>%
    summary(
      infer = c(TRUE, TRUE),
      adjust = guild_adjust
    ) %>%
    as_tibble() %>%
    rename(
      slope_difference = estimate,
      p_adjusted = p.value
    )
  
  
  #### Print results ####
  
  cat(
    "\n### 1. Overall artificial and natural trajectories ###\n"
  )
  print(overall_slopes, n = Inf)
  
  cat(
    "\n### Overall artificial-natural difference ###\n"
  )
  print(overall_contrast, n = Inf)
  
  cat(
    "\n### 2. Reef trajectories within each pair ###\n"
  )
  print(pair_slopes, n = Inf)
  
  cat(
    "\n### Artificial-natural difference within each pair ###\n"
  )
  print(pair_contrasts, n = Inf)
  
  cat(
    "\n### 3. Combined paired-abundance change ###\n"
  )
  print(paired_total_change, n = Inf)
  
  cat(
    "\n### 4. Feeding-guild trajectories ###\n"
  )
  print(guild_slopes, n = Inf)
  
  cat(
    "\n### Guild-specific artificial-natural differences ###\n"
  )
  print(guild_contrasts, n = Inf)
  
  
  list(
    overall_slopes = overall_slopes,
    overall_contrast = overall_contrast,
    pair_slopes = pair_slopes,
    pair_contrasts = pair_contrasts,
    paired_total_change = paired_total_change,
    guild_slopes = guild_slopes,
    guild_contrasts = guild_contrasts
  )
}


### 1.2. Run abundance extraction ####

abundance_answers <- extract_abundance_answers(
  model = abundance_final,
  data = abundance_data
)


### 1.3. Access individual results ####

abundance_answers$overall_slopes
abundance_answers$overall_contrast

abundance_answers$pair_slopes
abundance_answers$pair_contrasts

abundance_answers$paired_total_change

print(abundance_answers$guild_slopes, n=Inf)
abundance_answers$guild_contrasts


## 2.1 Extract Shannon diversity results ####

extract_diversity_answers <- function(
    model,
    data
) {
  
  #### 1. Artificial and natural diversity trajectories ####
  
  slope_emm <- emmeans::emtrends(
    model,
    ~ type,
    var = "date_s",
    weights = "equal"
  )
  
  slopes <- summary(
    slope_emm,
    infer = c(TRUE, TRUE),
    adjust = "none"
  ) %>%
    as_tibble() %>%
    rename(
      slope = date_s.trend,
      p_value = p.value
    ) %>%
    mutate(
      trajectory = case_when(
        p_value < 0.05 & slope > 0 ~ "increase",
        p_value < 0.05 & slope < 0 ~ "decrease",
        TRUE ~ "no significant change"
      )
    )
  
  
  #### 2. Artificial-natural difference in trajectories ####
  
  slope_contrast <- emmeans::contrast(
    slope_emm,
    method = list(
      "Artificial - Natural" = c(1, -1)
    )
  ) %>%
    summary(
      infer = c(TRUE, TRUE),
      adjust = "none"
    ) %>%
    as_tibble() %>%
    rename(
      slope_difference = estimate,
      p_value = p.value
    ) %>%
    mutate(
      trajectory_relationship = case_when(
        p_value >= 0.05 ~ "parallel trajectories",
        slope_difference > 0 ~
          "artificial trajectory more positive",
        slope_difference < 0 ~
          "natural trajectory more positive"
      )
    )
  
  
  #### 3. Artificial-natural difference at start and end ####
  
  time_values <- c(
    start = min(data$date_s),
    end = max(data$date_s)
  )
  
  type_emm <- emmeans::emmeans(
    model,
    ~ type | date_s,
    at = list(
      date_s = unname(time_values),
      period = "Post"
    ),
    weights = "equal"
  )
  
  type_estimates <- summary(
    update(type_emm, by = NULL),
    infer = c(TRUE, TRUE),
    adjust = "none"
  ) %>%
    as_tibble() %>%
    mutate(
      timepoint = case_when(
        near(date_s, time_values["start"]) ~ "Start",
        near(date_s, time_values["end"]) ~ "End",
        TRUE ~ as.character(date_s)
      )
    ) %>%
    rename(
      estimated_shannon = emmean,
      p_value = p.value
    )
  
  type_contrasts <- emmeans::contrast(
    type_emm,
    method = list(
      "Artificial - Natural" = c(1, -1)
    ),
    by = "date_s"
  ) %>%
    summary(
      infer = c(TRUE, TRUE),
      adjust = "none"
    ) %>%
    as_tibble() %>%
    mutate(
      timepoint = case_when(
        near(date_s, time_values["start"]) ~ "Start",
        near(date_s, time_values["end"]) ~ "End",
        TRUE ~ as.character(date_s)
      )
    ) %>%
    rename(
      diversity_difference = estimate,
      p_value = p.value
    )
  
  
  #### Print results ####
  
  cat(
    "\n### Artificial and natural Shannon diversity trajectories ###\n"
  )
  print(slopes, n = Inf)
  
  cat(
    "\n### Artificial-natural difference in trajectories ###\n"
  )
  print(slope_contrast, n = Inf)
  
  cat(
    "\n### Estimated Shannon diversity at the start and end ###\n"
  )
  print(type_estimates, n = Inf)
  
  cat(
    "\n### Artificial-natural difference at the start and end ###\n"
  )
  print(type_contrasts, n = Inf)
  
  
  list(
    slopes = slopes,
    slope_contrast = slope_contrast,
    type_estimates = type_estimates,
    type_contrasts = type_contrasts
  )
}
## 2.2 Run diversity extraction ####

diversity_answers <- extract_diversity_answers(
  model = shannon_final,
  data = diversity_data
)

## 2.3 Access individual results ####

diversity_answers$slopes
diversity_answers$slope_contrast
diversity_answers$type_estimates
diversity_answers$type_contrasts