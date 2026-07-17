## 01B_CLEAN_size.R #### 

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

clean_size_data <- function(file_path) {
  ## Load raw data ##
  raw_fish <- read.csv(file_path, stringsAsFactors = TRUE, strip.white = TRUE)
  
  ## Standardise site names and attach site metadata ##
  fish_size <- raw_fish %>%
    mutate(
      Site_raw = Site,
      Site = recode(
        Site,
        "Aow Mao Wall" = "Aow Mao",
        "No Name" = "No Name Pinnacle"
      ),
      Site = factor(Site)
    )
  
  site_lookup_size <- tibble::tribble(
    ~Site,              ~Type,        ~Pair,
    "Aow Mao",          "Natural",    "Aow Mao",
    "Aow Mao Wreck",    "Artificial", "Aow Mao",
    "Hin Pee Wee",      "Natural",    "Sattakut",
    "Sattakut",         "Artificial", "Sattakut",
    "No Name Pinnacle", "Natural",    "No Name",
    "No Name Wreck",    "Artificial", "No Name"
  )
  
  fish_size <- fish_size %>%
    left_join(site_lookup_size, by = "Site") %>%
    mutate(
      Type = factor(Type, levels = c("Artificial", "Natural")),
      Pair = factor(Pair),
      Species = factor(Species),
      Site = factor(Site)
    )
  
  ## Dates ##
  fish_size <- fish_size %>%
    mutate(
      Date = as.Date(as.character(Date_mm.dd.yy), format = "%m/%d/%Y"),
      Date = if_else(
        Date < as.Date("2024-01-01") & !is.na(Date),
        as.Date(paste0("2025-", format(Date, "%m-%d"))),
        Date
      ),
      Month_Year = format(Date, "%Y-%m")
    ) %>%
    select(-Date_mm.dd.yy, -Site_raw) %>%
    mutate(Month_Year = factor(Month_Year))
  
  ## Fix size bins ##
  fish_size <- fix_size_bins(fish_size)
  
  ## Structural size-bin NAs ##
  spp_max20 <- c(
    "Damsels - Regal Demoiselle",
    "Damsels - Alexanders"
  )
  
  spp_max50 <- c(
    "Parrotfish - Surf",
    "Rabbit - Java",
    "Butterfly - Weibels",
    "Butterfly - Longfin bannerfish",
    "Butterfly - Eight banded",
    "Angel - Blue-ringed"
  )
  
  fish_size <- fish_size %>%
    mutate(
      bin_20_50 = ifelse(Species %in% spp_max20, NA_integer_, bin_20_50),
      bin_50_100 = ifelse(Species %in% c(spp_max20, spp_max50),
                          NA_integer_, bin_50_100),
      bin_100 = ifelse(Species %in% c(spp_max20, spp_max50),
                       NA_integer_, bin_100)
    )
  ## Clean visibility - surveys with less than 3m vis excluded ##
  fish_size <- fish_size %>%
    mutate(
      Visibility_m = as.character(Visibility_m),
      Visibility_m = stringr::str_replace(Visibility_m, "m", ""),
      Visibility_m = as.numeric(Visibility_m)
    ) %>%
    filter(is.na(Visibility_m) | Visibility_m >= 3)
  
  ## Assign feeding guilds ##
  feeding_guilds <- tibble::tribble(
    ~Species,                          ~feeding_guild,
    "Angel - Blue-ringed",             "Invertivore",
    "Barracuda - Chevron",             "Mesopredator",
    "Barracuda - Great",               "HTLP",
    "Barracuda - Pickhandle",          "HTLP",
    "Barracuda - Yellowtail",          "Mesopredator",
    "Batfish - ALL",                   "Invertivore",
    "Butterfly - Eight banded",        "Invertivore",
    "Butterfly - Longfin bannerfish",  "Invertivore",
    "Butterfly - Weibels",             "Invertivore",
    "Cleaner - Blue-streaked",         "Invertivore",
    "Damsels - Alexanders",            "Herbivore",
    "Damsels - Regal Demoiselle",      "Herbivore",
    "Fusiliers - Yellowback",          "Herbivore",
    "Grouper - Blacktip",              "Mesopredator",
    "Grouper - Brown marbled",         "HTLP",
    "Grouper - Coral groupers (all)",  "HTLP",
    "Parrotfish - Surf",               "Herbivore",
    "Rabbit - Java",                   "Herbivore",
    "Rabbit - Virgate",                "Herbivore",
    "Snapper - Brownstripe",           "Mesopredator",
    "Snapper - Mangrove",              "HTLP",
    "Snapper - Russells",              "Mesopredator",
    "Sweetlips - Harlequin",           "Mesopredator",
    "Sweetlips - Harry hotlips",       "Invertivore",
    "Sweetlips - Painted",             "Mesopredator",
    "Tigger - Titan",                  "Invertivore",
    "Trigger - Titan",                 "Invertivore",
    "Trevally - Gold spotted",         "HTLP",
    "Trevally - Golden",               "HTLP",
    "Trevally - Orange spotted",       "HTLP",
    "Wrasse - Moon",                   "Invertivore",
    "Wrasse - Redbreasted",            "Invertivore"
  )
  
  fish_size <- fish_size %>%
    left_join(feeding_guilds, by = "Species") %>%
    mutate(
      feeding_guild = factor(
        feeding_guild,
        levels = c("Herbivore", "Invertivore", "Mesopredator", "HTLP")
      )
    )
  
  ## Fix scientific names ##
  fish_size <- fish_size %>%
    mutate(
      Sci_Name = case_when(
        Sci_Name == "Lutjanus griseus" ~ "Lutjanus argentimaculatus",
        Sci_Name == "Carangoides bajad" ~ "Flavocaranx bajad",
        Sci_Name == "Carangoides fulvoguttatus" ~ "Turrum fulvoguttatum",
        Sci_Name == "Neopomacentrus cyanos" ~ "Neopomacentrus cyanomos",
        Sci_Name == "Chaetodon wiebeli" ~ "Chaetodon weibeli",
        TRUE ~ Sci_Name
      )
    )
  
  ## Apply standardised naming ##
  fish_size <- fish_size %>%
    mutate(
      type = factor(Type, levels = c("Natural", "Artificial")),
      pair = factor(Pair),
      site = factor(Site),
      sci_name = factor(Sci_Name)
    )

  
  ## Pivot size bins to long format ##
  size_bins <- tibble::tribble(
    ~Size_Class, ~lower, ~upper,
    "0-1",       0,       1,
    "1-2",       1,       2,
    "2-5",       2,       5,
    "5-10",      5,       10,
    "10-15",     10,      15,
    "15-20",     15,      20,
    "20-50",     20,      50,
    "50-100",    50,      100,
    "100+",      100,     Inf
  )
  
  fish_long <- fish_size %>%
    pivot_longer(
      cols = c(bin_0_1, bin_1_2, bin_2_5, bin_5_10,
               bin_10_15, bin_15_20, bin_20_50,
               bin_50_100, bin_100),
      names_to = "Size_Bin",
      values_to = "Count"
    ) %>%
    mutate(
      Size_Class = recode(
        Size_Bin,
        "bin_0_1" = "0-1",
        "bin_1_2" = "1-2",
        "bin_2_5" = "2-5",
        "bin_5_10" = "5-10",
        "bin_10_15" = "10-15",
        "bin_15_20" = "15-20",
        "bin_20_50" = "20-50",
        "bin_50_100" = "50-100",
        "bin_100" = "100+"
      ),
      Size_Class = factor(Size_Class, levels = size_bins$Size_Class),
      Count = replace_na(Count, 0)
    ) %>%
    left_join(size_bins, by = "Size_Class")
  
  ## Create survey ID ##
  fish_long <- fish_long %>%
    mutate(
      survey_id = paste(site, Date, sep = "_"),
      survey_id = factor(survey_id)
    )
  
  ## Life stage: probabilistic ##
  fish_long_life_prob <- fish_long %>%
    left_join(maturity_lookup, by = "Sci_Name") %>%
    mutate(
      p_juv = case_when(
        is.na(Lmat_cm) ~ NA_real_,
        upper <= Lmat_cm ~ 1,
        lower >= Lmat_cm ~ 0,
        TRUE ~ (Lmat_cm - lower) / (upper - lower)
      ),
      p_adult = if_else(is.na(p_juv), NA_real_, 1 - p_juv),
      n_juv = round(Count * p_juv),
      n_adult = Count - n_juv
    ) %>%
    pivot_longer(
      cols = c(n_juv, n_adult),
      names_to = "life_stage",
      values_to = "stage_Count"
    ) %>%
    filter(stage_Count > 0) %>%
    mutate(
      life_stage = recode(
        life_stage,
        "n_juv" = "juvenile",
        "n_adult" = "adult"
      ),
      life_stage = factor(life_stage, levels = c("juvenile", "adult"))
    )
  
  ## Final selection and ordering ##
  fish_long_life_prob <- fish_long_life_prob %>%
    select(site, pair, type, Date, survey_id,
           Species, Sci_Name, sci_name, feeding_guild,
           Size_Class, lower, upper, Count, Lmat_cm,
           p_juv, p_adult, life_stage, stage_Count,
           Count.Type, Inclusion_m,
           Depth_m, Visibility_m, Current, Boats)
  
  ## Save cleaned file ##
  saveRDS(
    fish_long_life_prob,
    file.path(data_processed_dir, paste0("size_life_stage_probabilistic_", analysis_date, ".rds"))
  )
  
  message("Saved probabilistic juvenile/adult size data.")
  
  return(fish_long_life_prob)
}

## Run cleaning ##
fish_size <- clean_size_data(size_raw_path)

str(fish_size)
# and after 
fish_size %>%
  mutate(Species = factor(Species)) %>% 
  summarise(total_fish = sum(Count, na.rm = TRUE)) 

names(fish_size)
unique(fish_size$Species)
unique(fish_size$sci_name)
