### 06_PLOTS.R ###

#### 1. Prediction helpers ####

date_origin <- unique(abundance_data$Date - abundance_data$date_num)[1]
date_scaling <- coef(lm(date_s ~ date_num, data = abundance_data %>% distinct(date_num, date_s)))

to_date_s <- function(Date) {
  date_num <- as.numeric(Date - date_origin)
  date_scaling[1] + date_scaling[2] * date_num
}

predict_count_model <- function(model, newdata) {
  prediction <- predict(
    model, newdata = newdata, type = "link", se.fit = TRUE,
    re.form = NA, allow.new.levels = TRUE
  )
  
  newdata %>%
    mutate(
      predicted = exp(prediction$fit),
      lower = exp(prediction$fit - 1.96 * prediction$se.fit),
      upper = exp(prediction$fit + 1.96 * prediction$se.fit)
    )
}


#### 2. Figure 3: abundance trajectories ####

total_predictions <- total_abundance_data %>%
  group_by(pair, type) %>%
  reframe(Date = seq(min(Date), max(Date), length.out = 100)) %>%
  mutate(date_s = to_date_s(Date)) %>%
  predict_count_model(total_abundance_final, newdata = .)

guild_grid <- abundance_data %>%
  group_by(pair, type, feeding_guild) %>%
  reframe(Date = seq(min(Date), max(Date), length.out = 100)) %>%
  mutate(
    date_s = to_date_s(Date),
    survey_id = factor(levels(abundance_data$survey_id)[1], levels = levels(abundance_data$survey_id))
  )

guild_predictions <- predict_count_model(abundance_final, guild_grid)

fig3_total <- ggplot(total_predictions, aes(Date, predicted, colour = type, fill = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~pair, nrow = 1) +
  scale_colour_manual(values = reef_cols) +
  scale_fill_manual(values = reef_cols) +
  labs(x = NULL, y = "Predicted total abundance", colour = "Reef type", fill = "Reef type") +
  theme_clean +
  theme(legend.position = "none", strip.text = element_text(size = 11))

fig3_guild <- ggplot(guild_predictions, aes(Date, predicted, colour = type, fill = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1) +
  facet_grid(pair ~ feeding_guild, scales = "free_y", switch = "y") +
  scale_colour_manual(values = reef_cols) +
  scale_fill_manual(values = reef_cols) +
  labs(x = "Date", y = "Predicted feeding-guild abundance", colour = "Reef type", fill = "Reef type") +
  theme_clean +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    legend.position = "bottom"
  )

fig3_abundance <- fig3_total / fig3_guild +
  plot_layout(heights = c(0.65, 2), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

ggsave(
  file.path(figures_dir, "Fig3_Abundance.png"),
  fig3_abundance, width = 13, height = 5.5, dpi = 600, bg = "white"
)


#### 3. Figure 4: Shannon diversity ####

diversity_origin <- unique(diversity_data$Date - diversity_data$date_num)[1]
diversity_scaling <- coef(lm(date_s ~ date_num, data = diversity_data %>% distinct(date_num, date_s)))

diversity_grid <- diversity_data %>%
  group_by(pair, type) %>%
  reframe(Date = seq(min(Date), max(Date), length.out = 100)) %>%
  mutate(
    date_num = as.numeric(Date - diversity_origin),
    date_s = diversity_scaling[1] + diversity_scaling[2] * date_num,
    period = factor(
      if_else(pair == "Sattakut" | Date >= as.Date("2023-09-01"), "Post", "Pre"),
      levels = period_levels
    )
  )

diversity_prediction <- predict(shannon_final, newdata = diversity_grid, type = "response", se.fit = TRUE)

diversity_predictions <- diversity_grid %>%
  mutate(
    predicted = diversity_prediction$fit,
    lower = predicted - 1.96 * diversity_prediction$se.fit,
    upper = predicted + 1.96 * diversity_prediction$se.fit
  )

fig4a_diversity <- ggplot() +
  geom_vline(
    data = deployment_lookup %>% filter(!is.na(deployment_date)),
    aes(xintercept = deployment_date), linetype = 2, colour = "grey40"
  ) +
  geom_point(
    data = diversity_data, aes(Date, shannon, colour = type),
    alpha = 0.45, size = 1.8, position = position_jitter(width = 8, height = 0)
  ) +
  geom_ribbon(
    data = diversity_predictions, aes(Date, ymin = lower, ymax = upper, fill = type),
    alpha = 0.15, colour = NA
  ) +
  geom_line(
    data = diversity_predictions, aes(Date, predicted, colour = type),
    linewidth = 1
  ) +
  facet_wrap(~pair, nrow = 1) +
  scale_colour_manual(values = reef_cols) +
  scale_fill_manual(values = reef_cols) +
  labs(x = "Date", y = "Shannon diversity", colour = "Reef type", fill = "Reef type") +
  theme_clean +
  theme(legend.position = "bottom", strip.text = element_text(size = 11))

fig4b_diversity <- ggplot(diversity_data, aes(type, shannon, colour = type)) +
  geom_boxplot(width = 0.62, linewidth = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.14, height = 0, alpha = 0.5, size = 1.6) +
  scale_colour_manual(values = reef_cols) +
  labs(x = "Reef type", y = "Shannon diversity") +
  theme_clean +
  theme(legend.position = "none")

fig4_diversity <- (
  fig4a_diversity + guides(fill = "none") |
    fig4b_diversity + guides(colour = "none")
) +
  plot_layout(widths = c(3.2, 1)) +
  plot_annotation(tag_levels = "a")

ggsave(
  file.path(figures_dir, "Fig4_Shannon_diversity.png"),
  fig4_diversity, width = 13, height = 5.5, dpi = 600, bg = "white"
)


#### 4. Figure 5: community juvenile-adult balance ####

stage_balance_plot <- life_stage_contrasts %>%
  transmute(
    type, pair,
    adult_juvenile_ratio = ratio,
    adult_juvenile_lower = asymp.LCL,
    adult_juvenile_upper = asymp.UCL
  ) %>%
  mutate(
    juvenile_adult_ratio = 1 / adult_juvenile_ratio,
    juvenile_adult_lower = 1 / adult_juvenile_upper,
    juvenile_adult_upper = 1 / adult_juvenile_lower,
    pair = factor(pair, levels = pair_levels),
    type = factor(type, levels = reef_type_levels),
    pair_y = recode(as.character(pair), "Aow Mao" = 3, "No Name" = 2, "Sattakut" = 1),
    y = pair_y + if_else(type == "Natural", 0.11, -0.11),
    log2_ratio = log2(juvenile_adult_ratio),
    log2_lower = log2(juvenile_adult_lower),
    log2_upper = log2(juvenile_adult_upper)
  )

stage_balance_lines <- stage_balance_plot %>%
  select(pair, type, log2_ratio, y) %>%
  pivot_wider(names_from = type, values_from = c(log2_ratio, y))

fig5_juvenile_balance <- ggplot() +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45, colour = "grey65") +
  geom_segment(
    data = stage_balance_lines,
    aes(x = log2_ratio_Natural, xend = log2_ratio_Artificial, y = y_Natural, yend = y_Artificial),
    linewidth = 0.6, colour = "grey70"
  ) +
  geom_errorbarh(
    data = stage_balance_plot,
    aes(y = y, xmin = log2_lower, xmax = log2_upper),
    height = 0.05, linewidth = 0.7, colour = "grey25"
  ) +
  geom_point(
    data = stage_balance_plot, aes(log2_ratio, y, fill = type),
    shape = 21, size = 5, stroke = 0.8, colour = "black"
  ) +
  scale_fill_manual(
    values = scales::alpha(reef_cols, 0.75),
    labels = c("Natural" = "Natural reef", "Artificial" = "Artificial reef"),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = c(-1, 0, 1, 2),
    labels = c("2× adults", "Equal", "2× juveniles", "4× juveniles"),
    expand = expansion(mult = c(0.12, 0.12))
  ) +
  scale_y_continuous(
    breaks = c(1, 2, 3), labels = c("Sattakut", "No Name", "Aow Mao"),
    limits = c(0.55, 3.45)
  ) +
  labs(x = "Modelled juvenile-to-adult abundance ratio", y = NULL) +
  theme_classic(base_family = "Arial") +
  theme(
    axis.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(size = 10, margin = margin(t = 10)),
    legend.position = "bottom",
    panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey90"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(6, 8, 5, 6)
  )

ggsave(
  file.path(figures_dir, "Fig5_Juvenile_balance.png"),
  fig5_juvenile_balance, width = 8, height = 4.8, dpi = 600, bg = "white"
)


#### 5. Figure 6: species-specific juvenile proportions ####

species_order <- c(
  "Neopomacentrus cyanomos",
  "Pomacentrus alexanderae",
  "Lutjanus russellii",
  "Lutjanus argentimaculatus",
  "Lutjanus vitta"
)

species_strip_labels <- c(
  "Neopomacentrus cyanomos" = "Neopomacentrus\ncyanomos",
  "Pomacentrus alexanderae" = "Pomacentrus\nalexanderae",
  "Lutjanus russellii" = "Lutjanus\nrussellii",
  "Lutjanus argentimaculatus" = "Lutjanus\nargentimaculatus",
  "Lutjanus vitta" = "Lutjanus\nvitta"
)

species_labels <- focal_stage %>% distinct(Species, Sci_Name) %>% deframe()

raw_stage_props <- focal_stage %>%
  select(survey_id, Species, Sci_Name, pair, type, life_stage, stage_count) %>%
  pivot_wider(names_from = life_stage, values_from = stage_count, values_fill = 0) %>%
  mutate(total_count = juvenile + adult, raw_juvenile_prop = juvenile / total_count) %>%
  filter(total_count > 0)

juvenile_plot_data <- species_results$juvenile_proportions %>%
  mutate(Sci_Name = recode(as.character(Species), !!!species_labels)) %>%
  left_join(
    species_results$stage_contrasts %>%
      filter(contrast == "juvenile / adult") %>%
      transmute(
        Species, type, pair,
        juvenile_prop_lower = asymp.LCL / (1 + asymp.LCL),
        juvenile_prop_upper = asymp.UCL / (1 + asymp.UCL)
      ),
    by = c("Species", "type", "pair")
  )

composition_sig <- species_results$composition_contrasts %>%
  transmute(
    Sci_Name = recode(as.character(Species), !!!species_labels),
    pair, p_value = p.value, significant = p.value < 0.05
  )

raw_plot <- raw_stage_props %>%
  filter(Sci_Name %in% species_order) %>%
  mutate(
    Sci_Name = factor(Sci_Name, levels = species_order),
    pair = factor(pair, levels = pair_levels),
    type = factor(type, levels = reef_type_levels)
  )

model_plot <- juvenile_plot_data %>%
  filter(Sci_Name %in% species_order) %>%
  left_join(composition_sig, by = c("Sci_Name", "pair")) %>%
  mutate(
    Sci_Name = factor(Sci_Name, levels = species_order),
    pair = factor(pair, levels = pair_levels),
    type = factor(type, levels = reef_type_levels)
  )

sig_plot <- model_plot %>%
  filter(significant) %>%
  distinct(Sci_Name, pair, p_value) %>%
  mutate(x = 1.5, y = 0.8, label = paste0("p = ", format.pval(p_value, digits = 2)))

fig6_species_stage <- ggplot() +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.35, colour = "grey75") +
  geom_jitter(
    data = raw_plot, aes(type, raw_juvenile_prop),
    width = 0.07, height = 0, size = 0.8, alpha = 0.15, colour = "grey30"
  ) +
  geom_line(
    data = model_plot,
    aes(type, juvenile_proportion, group = interaction(Sci_Name, pair)),
    linewidth = 0.5, colour = "grey30"
  ) +
  geom_errorbar(
    data = model_plot,
    aes(type, ymin = juvenile_prop_lower, ymax = juvenile_prop_upper),
    width = 0.08, linewidth = 0.55, colour = "grey30"
  ) +
  geom_point(
    data = model_plot, aes(type, juvenile_proportion, fill = type),
    shape = 21, size = 2.6, stroke = 0.7, colour = "grey30"
  ) +
  geom_text(
    data = sig_plot, aes(x, y, label = label),
    inherit.aes = FALSE, size = 3, fontface = "bold"
  ) +
  facet_grid(
    rows = vars(Sci_Name), cols = vars(pair),
    labeller = labeller(Sci_Name = species_strip_labels), drop = FALSE
  ) +
  scale_fill_manual(values = scales::alpha(reef_cols, 0.7), guide = "none") +
  scale_x_discrete(labels = c("Natural" = "NR", "Artificial" = "AR")) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.08))
  ) +
  coord_cartesian(ylim = c(0, 1.07), clip = "off") +
  labs(x = NULL, y = "Juvenile proportion") +
  theme_classic(base_family = "Arial") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(face = "bold", size = 10),
    strip.text.y = element_text(
      angle = 0, face = "italic", size = 9.5, hjust = 0,
      lineheight = 0.9, margin = margin(l = 5)
    ),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 8.5),
    axis.title.y = element_text(size = 10, margin = margin(r = 10)),
    panel.spacing.x = unit(0.8, "lines"),
    panel.spacing.y = unit(0.65, "lines"),
    plot.margin = margin(6, 8, 5, 8)
  )

ggsave(
  file.path(figures_dir, "Fig6_Species_juvenile_proportions.png"),
  fig6_species_stage, width = 10, height = 10, dpi = 600, bg = "white"
)

message("Saved final figures to: ", figures_dir)
