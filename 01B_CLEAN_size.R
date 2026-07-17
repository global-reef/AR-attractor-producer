### 01B_CLEAN_size.R ###

#### 1. Lookups ####

size_guilds <- tibble::tribble(
  ~Species,                         ~feeding_guild,
  "Angel - Blue-ringed",            "Invertivore",
  "Barracuda - Chevron",            "Mesopredator",
  "Barracuda - Great",              "HTLP",
  "Barracuda - Pickhandle",         "HTLP",
  "Barracuda - Yellowtail",         "Mesopredator",
  "Batfish - ALL",                  "Invertivore",
  "Butterfly - Eight banded",       "Invertivore",
  "Butterfly - Longfin bannerfish", "Invertivore",
  "Butterfly - Weibels",            "Invertivore",
  "Cleaner - Blue-streaked",        "Invertivore",
  "Damsels - Alexanders",           "Herbivore",
  "Damsels - Regal Demoiselle",     "Herbivore",
  "Fusiliers - Yellowback",         "Herbivore",
  "Grouper - Blacktip",             "Mesopredator",
  "Grouper - Brown marbled",        "HTLP",
  "Grouper - Coral groupers (all)", "HTLP",
  "Parrotfish - Surf",              "Herbivore",
  "Rabbit - Java",                  "Herbivore",
  "Rabbit - Virgate",               "Herbivore",
  "Snapper - Brownstripe",          "Mesopredator",
  "Snapper - Mangrove",             "HTLP",
  "Snapper - Russells",             "Mesopredator",
  "Sweetlips - Harlequin",          "Mesopredator",
  "Sweetlips - Harry hotlips",      "Invertivore",
  "Sweetlips - Painted",            "Mesopredator",
  "Tigger - Titan",                 "Invertivore",
  "Trigger - Titan",                "Invertivore",
  "Trevally - Gold spotted",        "HTLP",
  "Trevally - Golden",              "HTLP",
  "Trevally - Orange spotted",      "HTLP",
  "Wrasse - Moon",                  "Invertivore",
  "Wrasse - Redbreasted",           "Invertivore"
)

size_bins <- tibble::tribble(
  ~Size_Bin,   ~Size_Class, ~lower, ~upper,
  "bin_0_1",   "0-1",            0,      1,
  "bin_1_2",   "1-2",            1,      2,
  "bin_2_5",   "2-5",            2,      5,
  "bin_5_10",  "5-10",           5,     10,
  "bin_10_15", "10-15",         10,     15,
  "bin_15_20", "15-20",         15,     20,
  "bin_20_50", "20-50",         20,     50,
  "bin_50_100","50-100",        50,    100,
  "bin_100",   "100+",          100,    Inf
)

spp_max20 <- c("Damsels - Regal Demoiselle", "Damsels - Alexanders")

spp_max50 <- c(
  "Parrotfish - Surf", "Rabbit - Java", "Butterfly - Weibels",
  "Butterfly - Longfin bannerfish", "Butterfly - Eight banded",
  "Angel - Blue-ringed"
)


#### 2. Cleaning function ####

clean_size_data <- function(file_path) {
  raw_fish <- read.csv(file_path, stringsAsFactors = FALSE, strip.white = TRUE, check.names = TRUE)
  
  fish_size <- raw_fish %>%
    mutate(
      Site = recode(Site, "Aow Mao Wall" = "Aow Mao", "No Name" = "No Name Pinnacle"),
      Date = as.Date(as.character(Date_mm.dd.yy), format = "%m/%d/%Y"),
      Date = if_else(
        !is.na(Date) & Date < as.Date("2024-01-01"),
        as.Date(paste0("2025-", format(Date, "%m-%d"))),
        Date
      ),
      visibility_raw = str_trim(as.character(Visibility_m)),
      visibility_below_3 = str_detect(visibility_raw, "^<\\s*3"),
      Visibility_m = readr::parse_number(visibility_raw)
    ) %>%
    filter(!visibility_below_3, is.na(Visibility_m) | Visibility_m >= 3) %>%
    select(-Date_mm.dd.yy, -visibility_raw, -visibility_below_3) %>%
    left_join(
      site_lookup %>% rename(Site = site, type = type, pair = pair),
      by = "Site"
    ) %>%
    filter(!is.na(pair)) %>%
    mutate(
      Month_Year = factor(format(Date, "%Y-%m")),
      site = factor(Site),
      pair = factor(pair, levels = pair_levels),
      type = factor(type, levels = reef_type_levels),
      Species = factor(Species)
    )
  
  fish_size <- fix_size_bins(fish_size) %>%
    mutate(
      bin_20_50 = if_else(Species %in% spp_max20, NA_integer_, bin_20_50),
      bin_50_100 = if_else(Species %in% c(spp_max20, spp_max50), NA_integer_, bin_50_100),
      bin_100 = if_else(Species %in% c(spp_max20, spp_max50), NA_integer_, bin_100)
    ) %>%
    left_join(size_guilds, by = "Species") %>%
    mutate(
      feeding_guild = factor(feeding_guild, levels = feeding_guild_levels),
      Sci_Name = recode(
        Sci_Name,
        "Lutjanus griseus" = "Lutjanus argentimaculatus",
        "Carangoides bajad" = "Flavocaranx bajad",
        "Carangoides fulvoguttatus" = "Turrum fulvoguttatum",
        "Neopomacentrus cyanos" = "Neopomacentrus cyanomos",
        "Chaetodon wiebeli" = "Chaetodon weibeli"
      ),
      sci_name = factor(Sci_Name)
    )
  
  fish_long <- fish_size %>%
    pivot_longer(all_of(size_bins$Size_Bin), names_to = "Size_Bin", values_to = "Count") %>%
    left_join(size_bins, by = "Size_Bin") %>%
    mutate(
      Size_Class = factor(Size_Class, levels = size_bins$Size_Class),
      Count = replace_na(Count, 0),
      survey_id = factor(paste(site, Date, sep = "_"))
    )
  
  fish_life_stage <- fish_long %>%
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
    pivot_longer(c(n_juv, n_adult), names_to = "life_stage", values_to = "stage_Count") %>%
    filter(stage_Count > 0) %>%
    mutate(
      life_stage = recode(life_stage, "n_juv" = "juvenile", "n_adult" = "adult"),
      life_stage = factor(life_stage, levels = life_stage_levels)
    ) %>%
    select(
      site, pair, type, Date, survey_id, Species, Sci_Name, sci_name, feeding_guild,
      Size_Class, lower, upper, Count, Lmat_cm, p_juv, p_adult, life_stage,
      stage_Count, Count.Type, Inclusion_m, Depth_m, Visibility_m, Current, Boats
    )
  
  stopifnot(
    all(fish_life_stage$stage_Count >= 0),
    all(fish_life_stage$stage_Count == floor(fish_life_stage$stage_Count)),
    !anyNA(fish_life_stage$feeding_guild)
  )
  
  saveRDS(
    fish_life_stage,
    file.path(data_processed_dir, paste0("size_life_stage_probabilistic_", analysis_date, ".rds"))
  )
  
  message(
    "Saved probabilistic life-stage data: ",
    n_distinct(fish_life_stage$survey_id), " surveys, ",
    n_distinct(fish_life_stage$Species), " species."
  )
  
  fish_life_stage
}


#### 3. Run cleaning ####

fish_size <- clean_size_data(size_raw_path)

fish_size %>%
  summarise(
    n_surveys = n_distinct(survey_id),
    n_species = n_distinct(Species),
    total_fish = sum(stage_Count)
  ) %>%
  print()

