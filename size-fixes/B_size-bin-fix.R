################################################################################
# B_size-bin-fix.R
#
# Project: AR-attractor-producer
# Purpose: Reconstruct early 10–20 cm size-bin data into 10–15 and 15–20 cm
#          bins using species-specific proportions from later surveys.
################################################################################

fix_size_bins <- function(fish_size) {
  
  message("Applying size-bin fix...")
  
  #### 1. Rename old raw size-bin columns if needed ####
  
  size_name_lookup <- c(
    "X0.1"      = "bin_0_1",
    "X.0.1"    = "bin_0_1",
    "X.1.2"    = "bin_1_2",
    "X.2.5"    = "bin_2_5",
    "X.5.10"   = "bin_5_10",
    "X.10.15"  = "bin_10_15",
    "X.15.20"  = "bin_15_20",
    "X.20.50"  = "bin_20_50",
    "X.50.100" = "bin_50_100",
    "X.100"    = "bin_100",
    "X.10.20"  = "bin_10_20"
  )
  
  present_old_names <- intersect(names(size_name_lookup), names(fish_size))
  
  names(fish_size)[match(present_old_names, names(fish_size))] <-
    unname(size_name_lookup[present_old_names])
  
  
  #### 2. Convert size-bin columns to integers ####
  
  clean_int <- function(x) {
    x_chr <- as.character(x)
    x_chr[x_chr %in% c("", " ", "NA", "N/A", "na", "n/a")] <- NA
    suppressWarnings(as.integer(x_chr))
  }
  
  size_cols <- c(
    "bin_0_1",
    "bin_1_2",
    "bin_2_5",
    "bin_5_10",
    "bin_10_15",
    "bin_15_20",
    "bin_20_50",
    "bin_50_100",
    "bin_100",
    "bin_10_20"
  )
  
  present_bins <- intersect(size_cols, names(fish_size))
  
  fish_size <- fish_size %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(present_bins), clean_int)
    )
  
  
  #### 3. Stop here if there is no old 10–20 cm bin ####
  
  if (!"bin_10_20" %in% names(fish_size)) {
    message("No old 10–20 cm bin found. Size-bin fix not needed.")
    return(fish_size)
  }
  
  if (!"bin_10_15" %in% names(fish_size)) {
    fish_size$bin_10_15 <- NA_integer_
  }
  
  if (!"bin_15_20" %in% names(fish_size)) {
    fish_size$bin_15_20 <- NA_integer_
  }
  
  fish_size <- fish_size %>%
    dplyr::mutate(bin_10_20_orig = bin_10_20)
  
  
  #### 4. Estimate species-specific 10–15 vs 15–20 proportions ####
  
  split_base <- fish_size %>%
    dplyr::filter(!is.na(bin_10_15), !is.na(bin_15_20)) %>%
    dplyr::mutate(total_10_20 = bin_10_15 + bin_15_20) %>%
    dplyr::filter(total_10_20 > 0)
  
  if (nrow(split_base) == 0) {
    warning("No later 10–15 / 15–20 data found. Old 10–20 cm bin was not split.")
    return(fish_size)
  }
  
  global_props <- split_base %>%
    dplyr::summarise(
      sum_10_15 = sum(bin_10_15, na.rm = TRUE),
      sum_15_20 = sum(bin_15_20, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      p10_global = sum_10_15 / (sum_10_15 + sum_15_20),
      p15_global = 1 - p10_global
    )
  
  p10_global <- global_props$p10_global[1]
  p15_global <- global_props$p15_global[1]
  
  species_props <- split_base %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(
      sum_10_15 = sum(bin_10_15, na.rm = TRUE),
      sum_15_20 = sum(bin_15_20, na.rm = TRUE),
      total_10_20 = sum_10_15 + sum_15_20,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      p10 = dplyr::if_else(total_10_20 > 0, sum_10_15 / total_10_20, NA_real_),
      p15 = 1 - p10
    )
  
  
  #### 5. Apply proportional split ####
  
  fish_size <- fish_size %>%
    dplyr::left_join(
      species_props %>% dplyr::select(Species, p10, p15),
      by = "Species"
    ) %>%
    dplyr::mutate(
      p10 = dplyr::coalesce(p10, p10_global),
      p15 = dplyr::coalesce(p15, p15_global),
      
      needs_10_20_split = !is.na(bin_10_20) &
        (is.na(bin_10_15) | is.na(bin_15_20)),
      
      bin_10_15 = dplyr::if_else(
        needs_10_20_split,
        as.integer(round(bin_10_20 * p10)),
        bin_10_15
      ),
      
      bin_15_20 = dplyr::if_else(
        needs_10_20_split,
        as.integer(bin_10_20 - bin_10_15),
        bin_15_20
      )
    ) %>%
    dplyr::select(-p10, -p15)
  
  
  #### 6. Remove old 10–20 cm modelling column ####
  # Keep bin_10_20_orig for provenance/QC.
  
  fish_size <- fish_size %>%
    dplyr::select(-bin_10_20)
  
  message("Size-bin fix complete.")
  
  return(fish_size)
}