library(dplyr)
library(tidyr)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(performance)

focal_taxa <- c(
  "Damsels - Regal Demoiselle",
  "Damsels - Alexanders",
  "Snapper - Russells",
  "Snapper - Brownstripe",
  "Snapper - Mangrove"
)

# aggregate 
stage_observed <- fish_size %>%
  filter(Species %in% focal_taxa) %>%
  group_by(
    survey_id, site, pair, type, Date,
    Species, Sci_Name, life_stage
  ) %>%
  summarise(
    stage_count = sum(stage_Count),
    .groups = "drop"
  )
stage_observed %>%
  summarise(
    n_noninteger = sum(stage_count %% 1 != 0),
    n_negative = sum(stage_count < 0),
    total = sum(stage_count)
  )

# restore structure 
survey_key <- fish_size %>%
  distinct(survey_id, site, pair, type, Date)

species_key <- fish_size %>%
  filter(Species %in% focal_taxa) %>%
  distinct(Species, Sci_Name)

focal_stage <- survey_key %>%
  crossing(
    species_key,
    life_stage = factor(
      c("juvenile", "adult"),
      levels = c("juvenile", "adult")
    )
  ) %>%
  left_join(
    stage_observed,
    by = c(
      "survey_id", "site", "pair", "type", "Date",
      "Species", "Sci_Name", "life_stage"
    )
  ) %>%
  mutate(
    stage_count = replace_na(stage_count, 0),
    stage_count = as.integer(stage_count),
    Species = factor(Species, levels = focal_taxa),
    pair = droplevels(pair),
    type = factor(type, levels = c("Natural", "Artificial")),
    life_stage = factor(
      life_stage,
      levels = c("adult", "juvenile")
    )
  )

focal_summary <- focal_stage %>%
  group_by(Species) %>%
  summarise(
    total_count = sum(stage_count),
    juvenile_count = sum(stage_count[life_stage == "juvenile"]),
    adult_count = sum(stage_count[life_stage == "adult"]),
    juvenile_surveys = sum(
      life_stage == "juvenile" & stage_count > 0
    ),
    adult_surveys = sum(
      life_stage == "adult" & stage_count > 0
    ),
    surveys_present = n_distinct(survey_id[stage_count > 0]),
    zero_pct = mean(stage_count == 0) * 100,
    .groups = "drop"
  )

focal_summary

focal_cells <- focal_stage %>%
  group_by(Species, pair, type, life_stage) %>%
  summarise(
    total_count = sum(stage_count),
    positive_surveys = sum(stage_count > 0),
    mean_count = mean(stage_count),
    maximum_count = max(stage_count),
    .groups = "drop"
  )

focal_cells



library(dplyr)
library(purrr)
library(tidyr)
library(glmmTMB)
library(DHARMa)
library(broom.mixed)

species_data <- focal_stage %>%
  mutate(
    Species = droplevels(Species),
    life_stage = factor(life_stage, levels = c("adult", "juvenile")),
    type = factor(type, levels = c("Natural", "Artificial")),
    pair = factor(pair)
  ) %>%
  group_split(Species) %>%
  set_names(map_chr(., ~ as.character(unique(.x$Species))))

fit_species_model <- function(dat) {
  glmmTMB(
    stage_count ~ life_stage * type * pair + (1 | survey_id),
    family = nbinom2,
    data = dat
  )
}

species_models <- map(species_data, fit_species_model)

fit_species_model <- function(dat) {
  glmmTMB(
    stage_count ~ life_stage * type * pair + (1 | survey_id),
    family = nbinom2,
    data = dat,
    control = glmmTMBControl(
      optimizer = optim,
      optArgs = list(method = "BFGS")
    )
  )
}

species_models_bfgs <- map(species_data, fit_species_model)
model_checks_bfgs <- imap_dfr(
  species_models_bfgs,
  \(mod, spp) tibble(
    Species = spp,
    convergence_code = mod$fit$convergence,
    convergence_message = mod$fit$message,
    positive_definite_hessian = mod$sdr$pdHess,
    singular = performance::check_singularity(mod)
  )
)

model_checks_bfgs


model_checks <- imap_dfr(
  species_models,
  \(mod, spp) tibble(
    Species = spp,
    convergence_code = mod$fit$convergence,
    convergence_message = mod$fit$message,
    positive_definite_hessian = mod$sdr$pdHess,
    singular = performance::check_singularity(mod)
  )
)

model_checks

model_parameters <- imap_dfr(
  species_models,
  \(mod, spp) {
    vc <- VarCorr(mod)$cond$survey_id
    
    tibble(
      Species = spp,
      survey_sd = attr(vc, "stddev"),
      dispersion = sigma(mod)
    )
  }
)

model_parameters

set.seed(123)

species_dharma <- map(
  species_models,
  ~ simulateResiduals(
    fittedModel = .x,
    n = 1000,
    plot = FALSE
  )
)
dharma_checks <- imap_dfr(
  species_dharma,
  \(res, spp) {
    uniformity <- testUniformity(res, plot = FALSE)
    dispersion <- testDispersion(res, plot = FALSE)
    zero_inflation <- testZeroInflation(res, plot = FALSE)
    outliers <- testOutliers(res, plot = FALSE)
    
    tibble(
      Species = spp,
      uniformity_p = uniformity$p.value,
      dispersion_p = dispersion$p.value,
      zero_inflation_p = zero_inflation$p.value,
      outlier_p = outliers$p.value
    )
  }
)

dharma_checks

species_model_summary <- model_checks %>%
  left_join(model_parameters, by = "Species") %>%
  left_join(dharma_checks, by = "Species")

species_model_summary

### for results 
extract_species_results <- function(models, adjust = "holm") {
  
  ## Omnibus Type III tests
  omnibus <- purrr::imap_dfr(
    models,
    \(mod, spp) {
      car::Anova(mod, type = 3) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("term") %>%
        transmute(
          Species = spp,
          term,
          chisq = Chisq,
          df = Df,
          p_value = `Pr(>Chisq)`,
          significance = if_else(
            p_value < 0.05,
            "Significant",
            "Non-significant"
          )
        )
    }
  )
  
  ## Estimated counts for every stage × type × pair combination
  emm <- purrr::map(
    models,
    ~ emmeans::emmeans(
      .x,
      ~ life_stage * type * pair
    )
  )
  
  predicted_counts <- purrr::imap_dfr(
    emm,
    \(x, spp) {
      summary(
        x,
        type = "response",
        infer = c(TRUE, TRUE)
      ) %>%
        as.data.frame() %>%
        mutate(Species = spp, .before = 1)
    }
  )
  
  ## Artificial versus natural within each life stage and pair
  reef_contrasts <- purrr::imap_dfr(
    emm,
    \(x, spp) {
      emmeans::contrast(
        x,
        method = "revpairwise",
        by = c("life_stage", "pair"),
        adjust = adjust
      ) %>%
        summary(
          type = "response",
          infer = c(TRUE, TRUE)
        ) %>%
        as.data.frame() %>%
        mutate(
          Species = spp,
          significance = if_else(
            p.value < 0.05,
            "Significant",
            "Non-significant"
          ),
          .before = 1
        )
    }
  )
  
  ## Juvenile versus adult within each reef type and pair
  stage_contrasts <- purrr::imap_dfr(
    emm,
    \(x, spp) {
      emmeans::contrast(
        x,
        method = "revpairwise",
        by = c("type", "pair"),
        adjust = adjust
      ) %>%
        summary(
          type = "response",
          infer = c(TRUE, TRUE)
        ) %>%
        as.data.frame() %>%
        mutate(
          Species = spp,
          significance = if_else(
            p.value < 0.05,
            "Significant",
            "Non-significant"
          ),
          .before = 1
        )
    }
  )
  
  ## Difference in juvenile-to-adult composition between reef types
  composition_contrasts <- purrr::imap_dfr(
    emm,
    \(x, spp) {
      emmeans::contrast(
        x,
        interaction = c("revpairwise", "revpairwise"),
        by = "pair",
        adjust = adjust
      ) %>%
        summary(
          type = "response",
          infer = c(TRUE, TRUE)
        ) %>%
        as.data.frame() %>%
        mutate(
          Species = spp,
          significance = if_else(
            p.value < 0.05,
            "Significant",
            "Non-significant"
          ),
          .before = 1
        )
    }
  )
  
  ## Derived juvenile proportions from predicted stage counts
  response_col <- intersect(
    c("response", "rate", "prob", "emmean"),
    names(predicted_counts)
  )[1]
  
  if (is.na(response_col)) {
    stop("Could not identify the predicted-response column.")
  }
  
  juvenile_proportions <- predicted_counts %>%
    select(
      Species, life_stage, type, pair,
      predicted_count = all_of(response_col)
    ) %>%
    pivot_wider(
      names_from = life_stage,
      values_from = predicted_count
    ) %>%
    mutate(
      juvenile_proportion = juvenile / (juvenile + adult)
    )
  
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

# Access each results table with 
species_results$omnibus
species_results$predicted_counts
species_results$reef_contrasts
species_results$stage_contrasts
species_results$juvenile_proportions
species_results$composition_contrasts


### plotting ###
library(dplyr)
library(ggplot2)
library(purrr)
library(cowplot)


composition_sig <- species_results$composition_contrasts %>%
  transmute(
    Sci_Name = recode(as.character(Species), !!!species_labels),
    pair,
    significant = p.value < 0.05
  )

plot_species_group <- function(
    species_names,
    raw_data = raw_stage_props,
    model_data = juvenile_plot_data,
    sig_data = composition_sig
) {
  
  raw_plot <- raw_data %>%
    filter(Sci_Name %in% species_names) %>%
    mutate(
      Sci_Name = factor(Sci_Name, levels = species_names),
      pair = factor(pair, levels = pair_levels),
      type = factor(type, levels = c("Natural", "Artificial"))
    )
  
  model_plot <- model_data %>%
    filter(Sci_Name %in% species_names) %>%
    left_join(sig_data, by = c("Sci_Name", "pair")) %>%
    mutate(
      Sci_Name = factor(Sci_Name, levels = species_names),
      pair = factor(pair, levels = pair_levels),
      type = factor(type, levels = c("Natural", "Artificial"))
    )
  
  sig_plot <- model_plot %>%
    distinct(Sci_Name, pair, significant) %>%
    filter(significant)
  
  p <- ggplot() +
    geom_hline(
      yintercept = 0.5,
      linetype = 2,
      linewidth = 0.35,
      colour = "grey70"
    ) +
    
    ## Raw survey-level proportions
    geom_jitter(
      data = raw_plot,
      aes(x = type, y = raw_juvenile_prop),
      width = 0.08,
      height = 0,
      size = 1,
      alpha = 0.25,
      colour = "grey35"
    ) +
    
    ## Modelled NR-to-AR difference
    geom_line(
      data = model_plot,
      aes(
        x = type,
        y = juvenile_proportion,
        group = interaction(Sci_Name, pair)
      ),
      linewidth = 0.8,
      colour = "grey20"
    ) +
    geom_point(
      data = model_plot,
      aes(
        x = type,
        y = juvenile_proportion,
        fill = type
      ),
      shape = 21,
      size = 2,
      stroke = 0.7,
      colour = "black"
    )
  
  ## Add significance markers only when present
  if (nrow(sig_plot) > 0) {
    p <- p +
      geom_text(
        data = sig_plot,
        aes(x = 1.5, y = 1.04, label = "*"),
        inherit.aes = FALSE,
        size = 5
      )
  }
  
  p +
    facet_grid(
      rows = vars(Sci_Name),
      cols = vars(pair),
      drop = FALSE
    ) +
    scale_fill_manual(
      values = reef_cols,
      guide = "none"
    ) +
    scale_x_discrete(
      labels = c(
        "Natural" = "NR",
        "Artificial" = "AR"
      )
    ) +
    scale_y_continuous(
      limits = c(0, 1.08),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = scales::percent,
      oob = scales::squish
    ) +
    labs(
      x = NULL,
      y = "Juvenile proportion"
    ) +
    theme_classic(base_family = "Arial") +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(
        face = "bold",
        size = 10
      ),
      strip.text.y = element_text(
        angle = 0,
        face = "italic",
        size = 10
      ),
      axis.text.x = element_text(size = 9),
      axis.title.y = element_text(
        margin = margin(r = 10)
      ),
      panel.spacing.x = unit(0.8, "lines"),
      panel.spacing.y = unit(0.7, "lines"),
      plot.margin = margin(5, 5, 5, 10)
    )
}

fig_snappers_stage <- plot_species_group(
  c(
    "Lutjanus russellii",
    "Lutjanus vitta",
    "Lutjanus argentimaculatus"
  )
)

fig_damsels_stage <- plot_species_group(
  c(
    "Neopomacentrus cyanomos",
    "Pomacentrus alexanderae"
  )
)

fig_snappers_stage
fig_damsels_stage


species_order <- c(
  "Neopomacentrus cyanomos",
  "Pomacentrus alexanderae",
  "Lutjanus russellii",
  "Lutjanus argentimaculatus",
  "Lutjanus vitta"
)

fig_species_stage <- plot_species_group(species_order)