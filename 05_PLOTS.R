## 05_PLOTS.R ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)


## 1. Prediction helpers ####

abundance_date_origin <- unique(abundance_data$Date - abundance_data$date_num)[1]

abundance_date_scaling <- coef(
  lm(date_s ~ date_num,
     data = abundance_data %>% distinct(date_num, date_s))
)

to_abundance_date_s <- function(Date) {
  date_num <- as.numeric(Date - abundance_date_origin)
  abundance_date_scaling[1] + abundance_date_scaling[2] * date_num
}

predict_count_model <- function(model, newdata) {
  prediction <- predict(model, newdata = newdata, type = "link",
                        se.fit = TRUE, re.form = NA, allow.new.levels = TRUE)
  
  newdata %>%
    mutate(
      predicted = exp(prediction$fit),
      lower = exp(prediction$fit - 1.96 * prediction$se.fit),
      upper = exp(prediction$fit + 1.96 * prediction$se.fit)
    )
}


## 2. Figure 3: abundance trajectories ####

abundance_predictions <- local({
  
  #### Total abundance ####
  
  total_grid <- total_abundance_data %>%
    group_by(pair, type) %>%
    reframe(Date = seq(min(Date), max(Date), length.out = 100)) %>%
    mutate(date_s = to_abundance_date_s(Date))
  
  total_predictions <- predict_count_model(total_abundance_final, total_grid)
  
  
  #### Feeding-guild abundance ####
  
  guild_grid <- abundance_data %>%
    group_by(pair, type, feeding_guild) %>%
    reframe(Date = seq(min(Date), max(Date), length.out = 100)) %>%
    mutate(
      date_s = to_abundance_date_s(Date),
      survey_id = factor(levels(abundance_data$survey_id)[1],
                         levels = levels(abundance_data$survey_id))
    )
  
  guild_predictions <- predict_count_model(abundance_final, guild_grid)
  
  list(total = total_predictions, guild = guild_predictions)
})


#### 2.1 Total abundance panel ####

fig3_total <- ggplot(abundance_predictions$total,
                     aes(Date, predicted, colour = type, fill = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~pair, nrow = 1) +
  scale_colour_manual(values = reef_cols) +
  scale_fill_manual(values = reef_cols) +
  labs(
    x = NULL,
    y = "Predicted total abundance",
    colour = "Reef type",
    fill = "Reef type"
  ) +
  theme_clean +
  theme(legend.position = "none",
        strip.text = element_text(size = 11))


#### 2.2 Feeding-guild panel ####

fig3_guild <- ggplot(abundance_predictions$guild,
                     aes(Date, predicted, colour = type, fill = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1) +
  facet_grid(pair ~ feeding_guild, scales = "free_y", switch = "y") +
  scale_colour_manual(values = reef_cols) +
  scale_fill_manual(values = reef_cols) +
  labs(
    x = "Date",
    y = "Predicted feeding-guild abundance",
    colour = "Reef type",
    fill = "Reef type"
  ) +
  theme_clean +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    legend.position = "bottom"
  )


#### 2.3 Combined Figure 3 ####

fig3_abundance <- fig3_total / fig3_guild +
  plot_layout(heights = c(0.65, 2), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

print(fig3_abundance)

ggsave(file.path(figures_dir, "Fig3_Abundance.png"),
       plot = fig3_abundance, width = 13, height = 5.5,
       dpi = 600, bg = "white")

rm(abundance_predictions)


## 3. Figure 4: Shannon diversity ####

diversity_date_origin <- unique(diversity_data$Date - diversity_data$date_num)[1]

diversity_date_scaling <- coef(
  lm(date_s ~ date_num,
     data = diversity_data %>% distinct(date_num, date_s))
)


### 3.1 Model predictions ####

diversity_predictions <- local({
  
  prediction_grid <- diversity_data %>%
    group_by(pair, type) %>%
    reframe(Date = seq(min(Date), max(Date), length.out = 100)) %>%
    mutate(
      date_num = as.numeric(Date - diversity_date_origin),
      date_s = diversity_date_scaling[1] + diversity_date_scaling[2] * date_num,
      period = case_when(
        pair == "Sattakut" ~ "Post",
        Date < as.Date("2023-09-01") ~ "Pre",
        TRUE ~ "Post"
      ),
      period = factor(period, levels = c("Pre", "Post"))
    )
  
  prediction <- predict(shannon_final, newdata = prediction_grid,
                        type = "response", se.fit = TRUE)
  
  prediction_grid %>%
    mutate(
      predicted = prediction$fit,
      lower = prediction$fit - 1.96 * prediction$se.fit,
      upper = prediction$fit + 1.96 * prediction$se.fit
    )
})


### 3.2 Panel A: temporal trajectories ####

fig4a_diversity <- ggplot() +
  geom_vline(
    data = deployment_lookup %>% filter(!is.na(deployment_date)),
    aes(xintercept = deployment_date),
    linetype = 2,
    colour = "grey40"
  ) +
  geom_point(
    data = diversity_data,
    aes(Date, shannon, colour = type),
    alpha = 0.45,
    size = 1.8,
    position = position_jitter(width = 8, height = 0)
  ) +
  geom_ribbon(
    data = diversity_predictions,
    aes(Date, ymin = lower, ymax = upper, fill = type),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(
    data = diversity_predictions,
    aes(Date, predicted, colour = type),
    linewidth = 1
  ) +
  facet_wrap(~pair, nrow = 1) +
  scale_colour_manual(values = reef_cols) +
  scale_fill_manual(values = reef_cols) +
  labs(
    x = "Date",
    y = "Shannon diversity",
    colour = "Reef type",
    fill = "Reef type"
  ) +
  theme_clean +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 11))


### 3.3 Panel B: overall reef-type distributions ####

fig4b_diversity <- ggplot(diversity_data,
                          aes(type, shannon, colour = type)) +
  geom_boxplot(width = 0.62, linewidth = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.14, height = 0, alpha = 0.50, size = 1.6) +
  scale_colour_manual(values = reef_cols) +
  labs(
    x = "Reef type",
    y = "Shannon diversity",
    colour = "Reef type"
  ) +
  theme_clean +
  theme(legend.position = "none")


### 3.4 Combine panels ####
fig4_diversity <- (fig4a_diversity + guides(fill = "none") | fig4b_diversity) +
  plot_layout(widths = c(3.2, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

fig4_diversity <- (
  fig4a_diversity + guides(fill = "none") +
    theme(legend.position = "bottom") |
    fig4b_diversity + guides(colour = "none")
) +
  plot_layout(widths = c(3.2, 1)) +
  plot_annotation(tag_levels = "a")

print(fig4_diversity)

ggsave(file.path(figures_dir, "Fig4_Shannon_diversity.png"),
       plot = fig4_diversity, width = 13, height = 5.5,
       dpi = 600, bg = "white")

rm(diversity_predictions)