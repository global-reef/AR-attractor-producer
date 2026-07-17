library(dplyr)
library(ggplot2)

pair_levels <- c("Aow Mao", "No Name", "Sattakut")

stage_balance_plot <- life_stage_contrasts %>%
  transmute(
    type,
    pair,
    adult_juvenile_ratio = ratio,
    adult_juvenile_lower = asymp.LCL,
    adult_juvenile_upper = asymp.UCL,
    p_value = p.value
  ) %>%
  mutate(
    juvenile_adult_ratio = 1 / adult_juvenile_ratio,
    juvenile_adult_lower = 1 / adult_juvenile_upper,
    juvenile_adult_upper = 1 / adult_juvenile_lower
  ) %>%
  left_join(
    stage_proportions %>%
      select(type, pair, total),
    by = c("type", "pair")
  ) %>%
  mutate(
    pair = factor(pair, levels = pair_levels),
    type = factor(type, levels = c("Natural", "Artificial")),
    pair_y = as.numeric(pair),
    y = pair_y + if_else(type == "Natural", 0.11, -0.11),
    log2_ratio = log2(juvenile_adult_ratio),
    log2_lower = log2(juvenile_adult_lower),
    log2_upper = log2(juvenile_adult_upper)
  ) %>%
  mutate(
    pair_y = case_when(
      pair == "Aow Mao" ~ 3,
      pair == "No Name" ~ 2,
      pair == "Sattakut" ~ 1
    ),
    y = pair_y + if_else(type == "Natural", 0.11, -0.11)
  )

stage_balance_lines <- stage_balance_plot %>%
  select(pair, type, log2_ratio, y) %>%
  tidyr::pivot_wider(
    names_from = type,
    values_from = c(log2_ratio, y)
  )




fig_juvenile_balance <- ggplot() +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.45,
    colour = "grey65"
  ) +
  
  geom_segment(
    data = stage_balance_lines,
    aes(
      x = log2_ratio_Natural,
      xend = log2_ratio_Artificial,
      y = y_Natural,
      yend = y_Artificial
    ),
    linewidth = 0.6,
    colour = "grey70"
  ) +
  
  geom_errorbarh(
    data = stage_balance_plot,
    aes(
      y = y,
      xmin = log2_lower,
      xmax = log2_upper
    ),
    height = 0.05,
    linewidth = 0.7,
    colour = "grey25"
  ) +
  
  geom_point(
    data = stage_balance_plot,
    aes(
      x = log2_ratio,
      y = y,
      fill = type
    ),
    shape = 21,
    size = 5,
    stroke = 0.8,
    colour = "black"
  ) +
  
  scale_fill_manual(
    values = scales::alpha(reef_cols, 0.75),
    labels = c(
      "Natural" = "Natural reef",
      "Artificial" = "Artificial reef"
    ),
    name = NULL
  ) +
  
  scale_x_continuous(
    breaks = c(-1, 0, 1, 2),
    labels = c(
      "2× adults",
      "Equal",
      "2× juveniles",
      "4× juveniles"
    ),
    expand = expansion(mult = c(0.12, 0.12))
  ) +
  scale_y_continuous(
    breaks = c(1, 2, 3),
    labels = c("Sattakut", "No Name", "Aow Mao"),
    limits = c(0.55, 3.45)
  ) + 
  labs(
    x = "Modelled juvenile-to-adult abundance ratio",
    y = NULL
  ) +
  
  theme_classic(base_family = "Arial") +
  
  theme(
    axis.text.y = element_text(
      face = "bold",
      size = 10
    ),
    
    axis.text.x = element_text(
      size = 9
    ),
    
    axis.title.x = element_text(
      size = 10,
      margin = margin(t = 10)
    ),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    
    panel.grid.major.y = element_line(
      linewidth = 0.3,
      colour = "grey90"
    ),
    
    panel.grid.minor = element_blank(),
    
    plot.margin = margin(
      t = 6,
      r = 8,
      b = 5,
      l = 6
    )
  )

fig_juvenile_balance
