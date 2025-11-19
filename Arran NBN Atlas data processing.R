library(tidyverse)
library(janitor)
library(stringr)

`%||%` <- function(x, y) if (is.null(x)) y else x

#----Converting NBN output to CSVs----
north_path <- "C:/Users/stuar/Documents/Masters/Arran/Data/NBN Atlas/North species list.txt"
south_path <- "C:/Users/stuar/Documents/Masters/Arran/Data/NBN Atlas/South Species_list_NBN Atlas.txt"
out_dir    <- "C:/Users/stuar/Documents/Masters/Arran/Data/NBN Atlas"
status_lookup <- file.path(out_dir, "status_lookup.csv")  # optional

read_nbn <- function(path, site_id) {
  x <- readr::read_csv(path, show_col_types = FALSE, locale = readr::locale(encoding = "UTF-8")) |>
    clean_names()
  
  if (!"lsid" %in% names(x)) stop("LSID column not found in: ", path)
  
  lsid <- str_split_fixed(x$lsid %||% "", pattern = "\\|", n = 5) |> as_tibble()
  names(lsid) <- c("lsid_sci_from_blob","nbn_key_from_blob","vernacular_from_blob","kingdom_from_blob","family_from_blob")
  x <- bind_cols(x, lsid)
  
  inv_raw <- x$invasive %||% ""
  inv_n   <- suppressWarnings(as.integer(str_extract(inv_raw, "\\d+")))
  inv_src <- str_match(inv_raw, '\\d+"?\\s*,\\s*"?([^"]*)"?$')[,2] |> na_if("")
  
  x <- x |>
    mutate(
      site_id = site_id,
      scientific_name = coalesce(species_name, lsid_sci_from_blob) |> str_squish(),
      vernacular_name = coalesce(vernacular_name, vernacular_from_blob) |> na_if(""),
      kingdom = coalesce(kingdom, kingdom_from_blob),
      family  = coalesce(family,  family_from_blob),
      nbn_key = nbn_key_from_blob,
      number_of_records = coalesce(suppressWarnings(as.integer(number_of_records)), inv_n),
      invasive_source = inv_src,
      invasive_flag = !is.na(invasive_source),
      genus = coalesce(genus, word(scientific_name, 1)),
      taxon_group = case_when(
        str_to_lower(class) %in% "aves" ~ "Birds",
        str_to_lower(class) %in% "mammalia" ~ "Mammals",
        str_to_lower(class) %in% "amphibia" ~ "Amphibians",
        str_to_lower(class) %in% "reptilia" ~ "Reptiles",
        str_to_lower(class) %in% c("actinopterygii","elasmobranchii") ~ "Fish",
        str_to_lower(class) %in% c("insecta","arachnida","malacostraca") ~ "Invertebrates",
        str_to_lower(kingdom) %in% "plantae" ~ "Plants",
        TRUE ~ "Other/Unknown"
      )
    ) |>
    filter(!is.na(scientific_name) & scientific_name != "") |>
    distinct()
  
  if (file.exists(status_lookup)) {
    st <- readr::read_csv(status_lookup, show_col_types = FALSE) |>
      clean_names() |>
      mutate(scientific_name = str_squish(scientific_name))
    x <- x |> left_join(st, by = "scientific_name")
  } else {
    x <- x |>
      mutate(
        eps_flag = NA, wca_schedule = NA, sbl_flag = NA,
        red_list_category = NA, annex_i_flag = NA, schedule_9 = NA
      )
  }
  
  rl_code <- str_extract(str_to_upper(x$red_list_category %||% ""), "CR|EN|VU|NT|LC|DD|NE")
  x <- x |>
    mutate(
      red_list_code = rl_code,
      endangered_flag = red_list_code %in% c("EN","CR"),
      threatened_flag = red_list_code %in% c("VU","EN","CR"),
      evidence_tier = "present"
    )
  
  x
}

nbn_north <- read_nbn(north_path, "North")
nbn_south <- read_nbn(south_path, "South")
nbn_all <- bind_rows(nbn_north, nbn_south)

readr::write_csv(nbn_north, file.path(out_dir, "nbn_species_clean_north.csv"))
readr::write_csv(nbn_south, file.path(out_dir, "nbn_species_clean_south.csv"))

summary_by_group <- nbn_all |>
  group_by(site_id, taxon_group) |>
  summarise(
    n_species = n_distinct(scientific_name),
    n_records_total = sum(number_of_records, na.rm = TRUE),
    n_inns = n_distinct(scientific_name[invasive_flag]),
    .groups = "drop"
  ) |>
  arrange(site_id, taxon_group)

important_species <- nbn_all |>
  filter(
    invasive_flag |
      !is.na(eps_flag) | !is.na(wca_schedule) | !is.na(sbl_flag) |
      !is.na(red_list_category) | !is.na(annex_i_flag) | !is.na(schedule_9) |
      threatened_flag
  ) |>
  select(
    site_id, scientific_name, vernacular_name, taxon_group,
    eps_flag, wca_schedule, sbl_flag, red_list_category, annex_i_flag, schedule_9,
    invasive_flag, invasive_source, endangered_flag, threatened_flag, number_of_records
  ) |>
  arrange(site_id, scientific_name)

endangered_species <- nbn_all |>
  filter(endangered_flag) |>
  select(site_id, scientific_name, vernacular_name, taxon_group, red_list_category, number_of_records) |>
  arrange(site_id, scientific_name)

readr::write_csv(summary_by_group, file.path(out_dir, "summary_by_taxon_group.csv"))
readr::write_csv(important_species, file.path(out_dir, "important_species.csv"))
readr::write_csv(endangered_species, file.path(out_dir, "endangered_species.csv"))

print(summary_by_group)
print(head(important_species, 15))
print(endangered_species)


##Plotting data##
#----Richness and abundance----
library(tidyverse)

out_dir <- "C:/Users/stuar/Documents/Masters/Arran/Data/NBN Atlas"

nbn <- bind_rows(
  readr::read_csv(file.path(out_dir, "nbn_species_clean_north.csv"), show_col_types = FALSE),
  readr::read_csv(file.path(out_dir, "nbn_species_clean_south.csv"), show_col_types = FALSE)
) |>
  mutate(number_of_records = as.integer(number_of_records)) |>
  mutate(number_of_records = replace_na(number_of_records, 0L))

rich <- nbn |>
  group_by(site_id) |>
  summarise(richness = n_distinct(scientific_name), .groups = "drop")

abund <- nbn |>
  group_by(site_id) |>
  summarise(abundance = sum(number_of_records, na.rm = TRUE), .groups = "drop")

p_rich <- ggplot(rich, aes(site_id, richness, fill = site_id)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  labs(x = NULL, y = "Species richness", title = "Species richness by site") +
  theme_minimal(base_size = 12)

p_abund <- ggplot(abund, aes(site_id, abundance, fill = site_id)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  labs(x = NULL, y = "Record abundance", title = "Record abundance by site") +
  theme_minimal(base_size = 12)

print(p_rich)
print(p_abund)

ggsave(file.path(out_dir, "plot_species_richness.png"), p_rich, width = 6.5, height = 4, dpi = 300)
ggsave(file.path(out_dir, "plot_record_abundance.png"), p_abund, width = 6.5, height = 4, dpi = 300)
