## 01A_CLEAN_abund.R ###

library(dplyr)
library(tidyr)
library(forcats)
library(lubridate)

clean_abundance_data <- function(file_path) {
  
  ## Load raw data ##
  raw_fish <- read.csv(file_path, stringsAsFactors = TRUE, strip.white = TRUE)
  
  ## Remove empty rows and columns ##
  raw_fish[raw_fish == ""] <- NA
  raw_fish <- raw_fish[, colSums(!is.na(raw_fish)) > 0]
  raw_fish <- raw_fish[rowSums(!is.na(raw_fish)) > 0, ]
  
  ## Identify species columns ##
  species_cols <- c("Parrotfish", "Rabbitfish", "Butterflyfish", "Angelfish",
                    "Cleaner_Wrasse", "Batfish", "Thicklip", "Red_Breast",
                    "Slingjaw", "Sweetlips", "Squirrel.Soldier", "Triggerfish",
                    "Porcupine.Puffer", "Ray", "Brown_Stripe_Snapper",
                    "Russels_Snapper", "lrg_Snapper", "Eel", "Trevally",
                    "Emperorfish", "sml_Grouper", "lrg_Grouper", "Barracuda")
  
  species_cols <- species_cols[species_cols %in% colnames(raw_fish)]
  
  ## Fix data types ##
  raw_fish[species_cols] <- lapply(raw_fish[species_cols], as.numeric)
  raw_fish$Date <- as.Date(as.character(raw_fish$Date), format = "%m/%d/%Y")
  
  ## Pivot to long format ##
  fish_long <- raw_fish %>%
    pivot_longer(
      cols = all_of(species_cols),
      names_to = "Species",
      values_to = "Count"
    )
  
  ## Assign feeding guilds ##
  feeding_guilds <- tibble::tribble(
    ~Species, ~feeding_guild,
    "Parrotfish", "Herbivore",
    "Rabbitfish", "Herbivore",
    "Butterflyfish", "Invertivore",
    "Angelfish", "Invertivore",
    "Cleaner_Wrasse", "Invertivore",
    "Batfish", "Invertivore",
    "Thicklip", "Invertivore",
    "Red_Breast", "Invertivore",
    "Slingjaw", "Invertivore",
    "Sweetlips", "Invertivore",
    "Squirrel.Soldier", "Invertivore",
    "Triggerfish", "Invertivore",
    "Porcupine.Puffer", "Invertivore",
    "Ray", "Mesopredator",
    "Brown_Stripe_Snapper", "Mesopredator",
    "Russels_Snapper", "Mesopredator",
    "lrg_Snapper", "HTLP",
    "Eel", "Mesopredator",
    "Trevally", "HTLP",
    "Emperorfish", "Mesopredator",
    "sml_Grouper", "Mesopredator",
    "lrg_Grouper", "HTLP",
    "Barracuda", "HTLP"
  )
  
  fish_long <- fish_long %>%
    left_join(feeding_guilds, by = "Species") %>%
    mutate(
      Count = ceiling(replace_na(Count, 0)),
      feeding_guild = factor(
        feeding_guild,
        levels = c("Herbivore", "Invertivore", "Mesopredator", "HTLP")
      )
    )
  
  ## Fix site naming ##
  fish_long$Site <- fct_recode(fish_long$Site, "No Name Pinnacle" = "No Name")
  
  ## Assign site pairs ##
  fish_long <- fish_long %>%
    mutate(
      pair = case_when(
        Site %in% c("Aow Mao", "Aow Mao Wreck") ~ "Aow Mao",
        Site %in% c("No Name Pinnacle", "No Name Wreck") ~ "No Name",
        Site %in% c("Hin Pee Wee", "Sattakut") ~ "Sattakut"
      )
    ) %>%
    filter(!is.na(pair))
  
  ## Check structure ##
  message("Site-pair structure:")
  print(fish_long %>% count(pair, Site))
  
  ## Deployment period flag ##
  deployment_date <- as.Date("2023-09-01")
  
  fish_long <- fish_long %>%
    mutate(
      deployment_period = case_when(
        pair == "Sattakut" ~ "Post",
        Date < deployment_date ~ "Pre",
        TRUE ~ "Post"
      ),
      Count = replace_na(Count, 0)
    )
  
  ## Compute months since deployment ##
  fish_long <- fish_long %>%
    mutate(
      months_since_deployment = if_else(
        Date < deployment_date,
        0,
        interval(deployment_date, Date) / months(1)
      )
    )
  
  ## Apply standardised naming ##
  fish_long <- fish_long %>%
    mutate(
      period = factor(deployment_period, levels = c("Pre", "Post")),
      type = factor(Type, levels = c("Artificial", "Natural")),
      t_since = pmax(0, months_since_deployment),
      pair = factor(pair),
      site = factor(Site),
      Species = factor(Species),
      feeding_guild = factor(
        feeding_guild,
        levels = c("Herbivore", "Invertivore", "Mesopredator", "HTLP")
      ),
      Date = ymd(Date),
      date_num = as.numeric(Date - min(Date))
    )
  
  fish_long$date_s <- as.numeric(scale(fish_long$date_num))
  
  ## Average observer rows for the same survey event ##
  fish_long <- fish_long %>%
    group_by(site, pair, type, period, t_since, date_num, Date, date_s,
             Species, feeding_guild) %>%
    summarise(
      Count = ceiling(mean(Count, na.rm = TRUE)),
      n_observer_rows = n(),
      n_observers = n_distinct(Researcher, na.rm = TRUE),
      researchers = paste(sort(unique(na.omit(as.character(Researcher)))),
                          collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(
      researchers = na_if(researchers, ""),
      survey_id = paste(site, Date, sep = "_"),
      survey_id = factor(survey_id)
    )
  
  ## Remove zero-inflated species after exploration ##
  fish_long <- fish_long %>%
    filter(Species != "Eel", Species != "Ray")
  
  ## Final selection and ordering ##
  fish_long <- fish_long %>%
    select(site, pair, type, period, t_since, date_num, Date, date_s,
           survey_id, researchers, n_observer_rows, n_observers,
           Species, feeding_guild, Count)
  
  ## Save cleaned file ##
  saveRDS(
    fish_long,
    file.path(data_processed_dir, paste0("abundance_species_cleaned_", analysis_date, ".rds"))
  )
  
  message("Saved cleaned abundance data.")
  
  return(fish_long)
}

## Run cleaning ##
fish_long <- clean_abundance_data(abundance_raw_path)

# restrict to Nov 2025 end 
fish_long <- fish_long %>%
  filter(Date <= as.Date("2025-11-30"))

str(fish_long)

# dataset for analysis 
fish_counts <- fish_long %>%
  group_by(site, pair, type, period, t_since, date_num, Date, date_s,
           survey_id, feeding_guild) %>%
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop")

# total fish before cleaning 
# abundance_raw %>%
#   summarise(total_fish = sum(total_N, na.rm = TRUE))
# and after 
fish_long %>%
  summarise(total_fish = sum(Count, na.rm = TRUE)) 

# observation hours 
fish_long %>%
  distinct(survey_id) %>%
  summarise(
    n_surveys = n(),
    observation_hours = n_surveys * 8 / 60
  )
