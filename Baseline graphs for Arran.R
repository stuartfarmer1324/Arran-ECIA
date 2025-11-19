library(tidyverse)
library(janitor)
library(stringr)
#----data formatting----
# data in (edit your path or use file.choose())
occurs_path <- "C:/Users/stuar/Documents/Masters/Arran/Data/associatedoccurences.csv"
# occurs_path <- file.choose()
stopifnot(file.exists(occurs_path))

# read + clean names
oc <- readr::read_csv(occurs_path, show_col_types = FALSE) |>
  clean_names()

# helper: coalesce from whichever name columns exist
coalesce_first <- function(df, candidates, new) {
  hit <- intersect(names(df), candidates)
  if (length(hit) == 0) {
    df[[new]] <- NA_character_
  } else {
    df[[new]] <- df[[hit[1]]]
    if (length(hit) > 1) for (i in hit[-1]) df[[new]] <- dplyr::coalesce(df[[new]], df[[i]])
  }
  df
}

# species-like test (binomial names; excludes sp./spp./rank words)
is_binomial <- function(x) {
  x <- str_trim(x)
  pat <- "^[A-Z][a-z]+\\s+[a-z][a-z\\-]+$"
  good <- str_detect(x, pat)
  bad  <- str_detect(x, regex("\\b(sp|spp|species|order|family|genus)\\b", ignore_case = TRUE))
  good & !bad
}

# standardise fields we need
if (!"taxonrank" %in% names(oc)) oc$taxonrank <- NA_character_
oc <- oc |>
  coalesce_first(c("scientificname","scientific_me","scientificme","scientific","taxon","verbatim_scientific_name"),
                 "scientific_name") |>
  coalesce_first(c("commonname","common_me","commonme","vernacular_name","common_name"),
                 "common_name") |>
  rename(order_name = any_of("order")) |>
  mutate(
    scientific_name    = scientific_name |> str_replace_all("_+", " ") |> str_squish(),
    common_name        = common_name     |> str_replace_all("_+", " ") |> str_squish(),
    occurrence_status  = str_to_sentence(str_squish(coalesce(occurrence_status, ""))),
    taxonrank          = str_to_sentence(str_squish(coalesce(taxonrank, ""))),
    kingdom            = str_to_sentence(str_squish(coalesce(kingdom, ""))),
    order_name         = str_squish(coalesce(order_name, "")),
    basis_of_record    = str_to_sentence(str_squish(coalesce(basis_of_record, ""))),
    hillside           = case_when(
      str_starts(tolower(event_id), "n") ~ "North",
      str_starts(tolower(event_id), "s") ~ "South",
      TRUE ~ NA_character_
    ),
    species_flag = case_when(
      !is.na(taxonrank) & tolower(taxonrank) == "species" ~ TRUE,
      !is.na(scientific_name) & is_binomial(scientific_name) ~ TRUE,
      TRUE ~ FALSE
    )
  )

# drop to present, species-like, hillside known
oc_sp <- oc |>
  filter(
    occurrence_status == "Present",
    species_flag,
    !is.na(hillside),
    !is.na(scientific_name),
    scientific_name != "",
    !str_detect(scientific_name, regex("^no records|unidentified", ignore_case = TRUE))
  )

# lightweight group mapping (rule-based; prioritised top->bottom)
mammal_orders <- c("artiodactyla","carnivora","rodentia","lagomorpha","eulipotyphla","chiroptera")  # bats handled earlier
aquatic_orders <- c("trichoptera","ephemeroptera","plecoptera","odonata","anisoptera","zygoptera")
terrestrial_orders <- c("lepidoptera","diptera","coleoptera","araneae","hymenoptera","hemiptera","thysanoptera","ixodida","julida","tipulidae","tipuloidea","homoptera")

# birds by event_id pattern (N1B, S2B, NB_I, SB_I, N3B_OoB, etc.)
bird_event_pat <- regex("^[NS](\\d+B(\\b|_)|B_)", ignore_case = TRUE)

oc_sp <- oc_sp |>
  mutate(
    oname_lc = tolower(order_name),
    cname_lc = tolower(coalesce(common_name, "")),
    sname_lc = tolower(coalesce(scientific_name, "")),
    bid = tolower(coalesce(basis_of_record, "")),
    evid = tolower(coalesce(event_id, "")),
    group = case_when(
      # plants
      kingdom == "Plantae" | str_detect(evid, "^[ns]\\d?p") ~ "Plants",
      # bats
      str_detect(oname_lc, "chiroptera") | str_detect(evid, "cr") | str_detect(bid, "audiomoth") ~ "Bats",
      # birds (use event IDs; taxonomic order often missing/messy)
      str_detect(evid, bird_event_pat) ~ "Birds",
      # aquatic inverts (orders + clear keywords + freshwater event IDs like *FI*)
      str_detect(evid, "fi") |
        str_detect(oname_lc, str_c(aquatic_orders, collapse = "|")) |
        str_detect(cname_lc, "caddis|mayfly|stonefly|dragonfly|damselfly|boatman|diving beetle|dytisc|corix") ~ "Aquatic invertebrates",
      # mammals (non-bat)
      str_detect(oname_lc, str_c(setdiff(mammal_orders, "chiroptera"), collapse = "|")) ~ "Mammals",
      # terrestrial inverts (orders + terrestrial sampling IDs)
      str_detect(oname_lc, str_c(terrestrial_orders, collapse = "|")) |
        str_detect(evid, "tin|tip|tit") ~ "Terrestrial invertebrates",
      TRUE ~ "Other"
    )
  )

# richness per group and hillside (distinct species per hillside)
rich_by_group <- oc_sp |>
  group_by(group, hillside) |>
  summarise(n_species = n_distinct(scientific_name), .groups = "drop")

# overall richness
rich_overall <- oc_sp |>
  group_by(hillside) |>
  summarise(n_species = n_distinct(scientific_name), .groups = "drop") |>
  mutate(group = "Overall")

plot_dat <- bind_rows(rich_by_group, rich_overall) |>
  filter(group %in% c("Birds","Bats","Mammals","Terrestrial invertebrates","Aquatic invertebrates","Plants","Overall"))

# order groups for plotting
group_levels <- c("Overall","Plants","Birds","Mammals","Bats","Terrestrial invertebrates","Aquatic invertebrates")
plot_dat <- plot_dat |>
  mutate(
    group = factor(group, levels = group_levels),
    hillside = factor(hillside, levels = c("North","South"))
  ) |>
  arrange(group, hillside)

#---- single figure: species richness by sampled group (+ overall), North vs South---- 
p <- ggplot(plot_dat, aes(x = n_species, y = group, fill = hillside)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = "white") +
  geom_text(aes(label = n_species),
            position = position_dodge(width = 0.75),
            hjust = -0.2, size = 3.8, colour = "grey20") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Species richness by sampled group and hillside",
    subtitle = "Distinct species per hillside (Present; species-level or binomial name)",
    x = "Species richness",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank()
  )

p 


#----Diverging richness plot----


# 1) Richness per group x hillside, plus overall
rich_by_group <- oc_sp |>
  dplyr::group_by(group, hillside) |>
  dplyr::summarise(n_species = dplyr::n_distinct(scientific_name), .groups = "drop")

rich_overall <- oc_sp |>
  dplyr::group_by(hillside) |>
  dplyr::summarise(n_species = dplyr::n_distinct(scientific_name), .groups = "drop") |>
  dplyr::mutate(group = "Overall")

# 2) Keep the groups we want 
groups_to_show <- c("Overall","Plants","Birds","Mammals","Bats","Terrestrial invertebrates","Aquatic invertebrates")
summ <- dplyr::bind_rows(rich_by_group, rich_overall) |>
  dplyr::filter(group %in% groups_to_show) |>
  tidyr::pivot_wider(names_from = hillside, values_from = n_species, values_fill = 0) |>
  dplyr::mutate(
    diff   = coalesce(South, 0) - coalesce(North, 0),           
    winner = dplyr::case_when(diff > 0 ~ "South", diff < 0 ~ "North", TRUE ~ "Tie"),
    group  = forcats::fct_reorder(group, abs(diff), .desc = TRUE)
  )

# 3) Symmetric x-limits for a clean diverging look
lim <- max(abs(summ$diff), na.rm = TRUE)
lims <- c(-lim * 1.15, lim * 1.15)
offset <- max(0.4, lim * 0.03)  

# 4) Plot
p_diverge <- ggplot(summ, aes(x = diff, y = group, fill = winner)) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey65") +
  geom_col(width = 0.7, colour = "white") +
  geom_text(aes(x = pmin(diff, 0) - offset, label = paste0("N=", coalesce(North, 0))),
            hjust = 1, size = 3.6, colour = "grey20") +
  geom_text(aes(x = pmax(diff, 0) + offset, label = paste0("S=", coalesce(South, 0))),
            hjust = 0, size = 3.6, colour = "grey20") +
  scale_fill_manual(values = c("North" = "#1f78b4", "South" = "#33a02c", "Tie" = "grey70")) +
  scale_x_continuous(limits = lims, breaks = scales::pretty_breaks()) +
  labs(
    title = "Species richness (diverging) by sampled group",
    subtitle = "Bar direction shows which hillside has higher richness; labels show North and South totals",
    x = "South − North (species)", y = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank()
  )

p_diverge  


#----abundance----
rank_abund <- oc_sp |>
  count(hillside, scientific_name, name = "records") |>
  group_by(hillside) |>
  arrange(desc(records), .by_group = TRUE) |>
  mutate(rank = row_number()) |>
  ungroup()

ggplot(rank_abund, aes(x = rank, y = records, colour = hillside)) +
  geom_line() + geom_point(size = 1.5) +
  labs(title = "Rank–abundance (by record frequency)",
       x = "Species rank (per hillside)", y = "Records (log scale)", colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")


#----simple turnover----
jaccard_by_group <- oc_sp |>
  distinct(group, hillside, scientific_name) |>
  group_by(group) |>
  summarise(
    north = n_distinct(scientific_name[hillside == "North"]),
    south = n_distinct(scientific_name[hillside == "South"]),
    shared = n_distinct(intersect(
      scientific_name[hillside == "North"],
      scientific_name[hillside == "South"]
    )),
    jaccard = ifelse(north + south - shared == 0, NA_real_, shared / (north + south - shared)),
    .groups = "drop"
  ) |>
  filter(group %in% group_levels)

ggplot(jaccard_by_group, aes(x = jaccard, y = factor(group, levels = rev(group_levels)))) +
  geom_vline(xintercept = 0.5, linetype = 3, colour = "grey70") +
  geom_point(size = 3, colour = "#2c7fb8") +
  scale_x_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Similarity between hillsides by group",
       subtitle = "Jaccard index = shared / (North + South − shared) [species richness]",
       x = "Similarity", y = NULL) +
  theme_minimal(base_size = 12)


# ---- Two plots (richness and rank–abundance) per group; paste below your existing code ----

make_group_plots <- function(df, grp_label) {
  df_grp <- df |> dplyr::filter(group == grp_label)
  
  if (nrow(df_grp) == 0) return(invisible(NULL))
  
  # Richness per hillside (distinct species)
  rich <- df_grp |>
    dplyr::distinct(hillside, scientific_name) |>
    dplyr::count(hillside, name = "n_species")
  
  p_rich <- ggplot2::ggplot(rich, ggplot2::aes(x = hillside, y = n_species, fill = hillside)) +
    ggplot2::geom_col(width = 0.7, colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = n_species), vjust = -0.3, size = 4) +
    ggplot2::labs(
      title = paste0(grp_label, " — species richness by hillside"),
      x = NULL, y = "Species richness"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none")
  
  print(p_rich)
  
  # Rank–abundance (by record frequency) per hillside
  rank_abund <- df_grp |>
    dplyr::count(hillside, scientific_name, name = "records") |>
    dplyr::group_by(hillside) |>
    dplyr::arrange(dplyr::desc(records), .by_group = TRUE) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::ungroup()
  
  if (nrow(rank_abund) == 0) return(invisible(NULL))
  
  p_abund <- ggplot2::ggplot(rank_abund, ggplot2::aes(x = rank, y = records, colour = hillside)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      title = paste0(grp_label, " — rank–abundance (by record frequency)"),
      x = "Species rank (per hillside)",
      y = "Records (log scale)",
      colour = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top")
  
  print(p_abund)
}

# Order requested: Birds first, then invertebrates etc.
groups_to_plot <- c("Birds", "Terrestrial invertebrates", "Aquatic invertebrates", "Mammals", "Bats", "Plants")
groups_to_plot <- groups_to_plot[groups_to_plot %in% unique(oc_sp$group)]

for (g in groups_to_plot) {
  make_group_plots(oc_sp, g)
}


#----Species specific differences between plots---



# ---- Individual species abundance per group (North vs South) — paste below your existing code ----
make_species_abundance <- function(df, grp_label, top_n = Inf) {
  df_grp <- df |>
    dplyr::filter(group == grp_label)
  
  if (nrow(df_grp) == 0) return(invisible(NULL))
  
  # Count records per species x hillside and ensure both hillsides appear per species (fill missing with 0)
  dat0 <- df_grp |>
    dplyr::count(hillside, scientific_name, name = "records")
  
  dat <- dat0 |>
    tidyr::complete(scientific_name, hillside = c("North","South"), fill = list(records = 0)) |>
    dplyr::mutate(hillside = factor(hillside, levels = c("North","South")))
  
  # Order species by total records (across both hillsides)
  totals <- dat |>
    dplyr::group_by(scientific_name) |>
    dplyr::summarise(total_records = sum(records), .groups = "drop")
  
  dat <- dat |>
    dplyr::left_join(totals, by = "scientific_name")
  
  # Optionally limit to top_n species for readability (keep all by default)
  if (is.finite(top_n)) {
    keep <- totals |>
      dplyr::arrange(dplyr::desc(total_records)) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::pull(scientific_name)
    dat <- dat |>
      dplyr::filter(scientific_name %in% keep)
  }
  
  dat <- dat |>
    dplyr::mutate(scientific_name = reorder(scientific_name, total_records))
  
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = scientific_name, y = records, fill = hillside)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7, colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = records),
                       position = ggplot2::position_dodge(width = 0.75),
                       hjust = -0.1, size = 3, colour = "grey25") +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::labs(
      title = paste0(grp_label, " — individual species records by hillside"),
      x = NULL, y = "Record count", fill = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top")
  
  print(p)
}

# Requested order: Birds first, then Mammals, Bats, Terrestrial inverts, Aquatic inverts, Plants
groups_to_plot <- c("Birds", "Mammals", "Bats", "Terrestrial invertebrates", "Aquatic invertebrates", "Plants")
groups_to_plot <- groups_to_plot[groups_to_plot %in% unique(oc_sp$group)]

# Make plots (set top_n to a number if too many species; by default shows all)
for (g in groups_to_plot) {
  make_species_abundance(oc_sp, g)          # all species
  # make_species_abundance(oc_sp, g, top_n = 25)  # uncomment to show top 25 per group
}




###Important species graphs###

# Ensure fixed, shared y-scales across North/South facets

# Richness (protected/important only)
ymax_rich <- max(imp_rich$richness, na.rm = TRUE)

p_imp_rich <- ggplot(imp_rich, aes(x = importance_class, y = richness, fill = importance_class)) +
  geom_col(width = 0.7, colour = "white") +
  facet_wrap(~ hillside, nrow = 1) +
  scale_fill_manual(values = imp_cols, guide = "none") +
  scale_y_continuous(
    limits = c(0, ymax_rich * 1.1),
    expand = expansion(mult = c(0, 0.05)),
    breaks = scales::pretty_breaks()
  ) +
  labs(title = "Important/protected species richness",
       subtitle = "Distinct species per hillside (fixed y-scale across facets)",
       x = NULL, y = "Species richness") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
p_imp_rich

# Abundance (protected/important only)
ymax_abund <- max(imp_abund$total_records, na.rm = TRUE)

p_imp_abund <- ggplot(imp_abund, aes(x = importance_class, y = total_records, fill = importance_class)) +
  geom_col(width = 0.7, colour = "white") +
  facet_wrap(~ hillside, nrow = 1) +  # fixed scales (no free_y)
  scale_fill_manual(values = imp_cols, guide = "none") +
  scale_y_continuous(
    limits = c(0, ymax_abund * 1.1),
    expand = expansion(mult = c(0, 0.05)),
    breaks = scales::pretty_breaks()
  ) +
  labs(title = "Important/protected species record abundance",
       subtitle = "Total records per hillside (fixed y-scale across facets)",
       x = NULL, y = "Records") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
p_imp_abund


# ---- species-level comparison (only important/protected), by hillside ----
imp_species <- imp %>%
  count(hillside, importance_class, scientific_name, name = "records") %>%
  group_by(scientific_name, importance_class) %>%
  mutate(total = sum(records)) %>%
  ungroup() %>%
  mutate(scientific_name = fct_reorder(scientific_name, total))

p_imp_species <- ggplot(imp_species, aes(x = scientific_name, y = records, fill = hillside)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, colour = "white") +
  coord_flip(clip = "off") +
  facet_wrap(~ importance_class, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("North" = "#1f78b4", "South" = "#33a02c")) +
  labs(title = "Important/protected species — North vs South",
       x = NULL, y = "Records", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"))
p_imp_species

