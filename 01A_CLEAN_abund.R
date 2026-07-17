### 01A_CLEAN_abund.R ###

#### 1. Lookups ####

abundance_species <- c(
  "Parrotfish", "Rabbitfish", "Butterflyfish", "Angelfish", "Cleaner_Wrasse",
  "Batfish", "Thicklip", "Red_Breast", "Slingjaw", "Sweetlips",
  "Squirrel.Soldier", "Triggerfish", "Porcupine.Puffer", "Ray",
  "Brown_Stripe_Snapper", "Russels_Snapper", "lrg_Snapper", "Eel",
  "Trevally", "Emperorfish", "sml_Grouper", "lrg_Grouper", "Barracuda"
)

abundance_guilds <- tibble::tribble(
  ~Species,                ~feeding_guild,
  "Parrotfish",            "Herbivore",
  "Rabbitfish",            "Herbivore",
  "Butterflyfish",         "Invertivore",
  "Angelfish",             "Invertivore",
  "Cleaner_Wrasse",        "Invertivore",
  "Batfish",               "Invertivore",
  "Thicklip",              "Invertivore",
  "Red_Breast",            "Invertivore",
  "Slingjaw",              "Invertivore",
  "Sweetlips",             "Invertivore",
  "Squirrel.Soldier",      "Invertivore",
  "Triggerfish",           "Invertivore",
  "Porcupine.Puffer",      "Invertivore",
  "Ray",                   "Mesopredator",
  "Brown_Stripe_Snapper",  "Mesopredator",
  "Russels_Snapper",       "Mesopredator",
  "lrg_Snapper",           "HTLP",
  "Eel",                   "Mesopredator",
  "Trevally",              "HTLP",
  "Emperorfish",           "Mesopredator",
  "sml_Grouper",           "Mesopredator",
  "lrg_Grouper",           "HTLP",
  "Barracuda",             "HTLP"
)


#### 2. Cleaning function ####

clean_abundance_data <- function(file_path, end_date = as.Date("2025-11-30")) {
  raw_fish <- read.csv(file_path, stringsAsFactors = FALSE, strip.white = TRUE, check.names = FALSE)
  raw_fish[raw_fish == ""] <- NA
  raw_fish <- raw_fish[rowSums(!is.na(raw_fish)) > 0, colSums(!is.na(raw_fish)) > 0, drop = FALSE]
  
  species_cols <- intersect(abundance_species, names(raw_fish))
  raw_fish[species_cols] <- lapply(raw_fish[species_cols], \(x) suppressWarnings(as.numeric(x)))
  
  fish_long <- raw_fish %>%
    mutate(
      Date = as.Date(Date, "%m/%d/%Y"),
      Site = recode(Site, "No Name" = "No Name Pinnacle")
    ) %>%
    filter(!is.na(Date), Date <= end_date) %>%
    pivot_longer(all_of(species_cols), names_to = "Species", values_to = "Count") %>%
    mutate(Count = ceiling(replace_na(Count, 0))) %>%
    left_join(abundance_guilds, by = "Species") %>%
    left_join(
      site_lookup %>% rename(Site = site, Type_lookup = type, pair = pair),
      by = "Site"
    ) %>%
    filter(!is.na(pair), !Species %in% c("Eel", "Ray")) %>%
    mutate(
      Type = coalesce(Type_lookup, as.character(Type)),
      period = if_else(pair == "Sattakut" | Date >= as.Date("2023-09-01"), "Post", "Pre"),
      t_since = pmax(
        0,
        lubridate::time_length(
          lubridate::interval(as.Date("2023-09-01"), Date),
          unit = "month"
        )
      ),
      date_num = as.numeric(Date - min(Date, na.rm = TRUE)),
      type = factor(Type, levels = reef_type_levels),
      period = factor(period, levels = period_levels),
      pair = factor(pair, levels = pair_levels),
      site = factor(Site),
      Species = factor(Species),
      feeding_guild = factor(feeding_guild, levels = feeding_guild_levels)
    ) %>%
    mutate(date_s = as.numeric(scale(date_num))) %>%
    group_by(site, pair, type, period, t_since, date_num, Date, date_s, Species, feeding_guild) %>%
    summarise(
      Count = ceiling(mean(Count, na.rm = TRUE)),
      n_observer_rows = n(),
      n_observers = n_distinct(Researcher, na.rm = TRUE),
      researchers = paste(sort(unique(na.omit(Researcher))), collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(
      researchers = na_if(researchers, ""),
      survey_id = factor(paste(site, Date, sep = "_"))
    ) %>%
    select(
      site, pair, type, period, t_since, date_num, Date, date_s, survey_id,
      researchers, n_observer_rows, n_observers, Species, feeding_guild, Count
    )
  
  stopifnot(
    !anyNA(fish_long$Count),
    all(fish_long$Count >= 0),
    all(fish_long$Count == floor(fish_long$Count)),
    !anyNA(fish_long$feeding_guild)
  )
  
  saveRDS(
    fish_long,
    file.path(data_processed_dir, paste0("abundance_species_cleaned_", analysis_date, ".rds"))
  )
  
  message(
    "Saved cleaned abundance data: ",
    n_distinct(fish_long$survey_id), " surveys, ",
    sum(fish_long$Count), " fish."
  )
  
  fish_long
}


#### 3. Run cleaning and construct analysis data ####

fish_long <- clean_abundance_data(abundance_raw_path)

fish_counts <- fish_long %>%
  group_by(site, pair, type, period, t_since, date_num, Date, date_s, survey_id, feeding_guild) %>%
  summarise(Count = sum(Count), .groups = "drop")

fish_long %>%
  summarise(
    n_surveys = n_distinct(survey_id),
    total_fish = sum(Count),
    observation_hours = n_surveys * 8 / 60
  ) %>%
  print()