# Arran biodiversity analysis + important species pipeline
# Uses associatedoccurences.csv as the data source
# Imputes genus-only records to synthetic species without altering real species
# Produces all-taxa summaries plus important-species extracts

# ---- Libraries ----

required_pkgs <- c("tidyverse", "stringr", "fuzzyjoin", "stringdist", "forcats", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Install the missing packages before running: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

library(tidyverse)
library(stringr)
library(fuzzyjoin)
library(stringdist)
library(forcats)
library(scales)
# library(patchwork) # uncomment if you want multipanel layouts

# ---- Load occurrences ----

# Point to your occurrences CSV. If missing, fall back to a tiny demo dataset so
# the script still runs for testing/documentation.
occurs_path <- Sys.getenv(
  "OCCURS_PATH",
  unset = "C:/Users/stuar/Documents/Masters/Arran/Data/associatedoccurences.csv"
)

if (file.exists(occurs_path)) {
  message("Loading occurrences from: ", occurs_path)
  assoc <- readr::read_csv(occurs_path, show_col_types = FALSE)
} else {
  message(
    "Occurrences file not found at ", occurs_path, "; using built-in demo data. ",
    "Set OCCURS_PATH to your CSV to run with the real survey."
  )
  
  assoc <- tibble(
    eventID        = c("N1B", "S1B", "N2M", "S2M", "N3P", "S3P", "N4T", "S4A"),
    Commonme       = c("curlew", "curlew", "otter", "otter",
                       "scots pine", "scots pine", "bumblebee", "azure damselfly"),
    scientificme   = c("Numenius arquata", "Numenius arquata", "Lutra lutra", "Lutra lutra",
                       "Pinus sylvestris", "Pinus sylvestris", "Bombus sp.", "Coenagrion puella"),
    kingdom        = c("Animalia", "Animalia", "Animalia", "Animalia",
                       "Plantae", "Plantae", "Animalia", "Animalia"),
    Order          = c("Charadriiformes", "Charadriiformes", "Carnivora", "Carnivora",
                       "Pinales", "Pinales", "Hymenoptera", "Odonata"),
    taxonRank      = rep("Species", 8),
    occurrenceStatus = rep("Present", 8),
    basisOfRecord  = c("HumanObservation", "HumanObservation", "HumanObservation",
                       "HumanObservation", "HumanObservation", "HumanObservation",
                       "HumanObservation", "HumanObservation")
  )
}

# ---- Standardise fields ----

assoc <- assoc %>%
  rename(
    event_id        = eventID,
    common_name_raw = Commonme,
    scientific_raw  = scientificme
  ) %>%
  mutate(
    common_name_raw = str_squish(common_name_raw),
    scientific_raw  = str_squish(scientific_raw),
    common_name     = str_to_lower(common_name_raw),
    scientific_name = str_to_lower(scientific_raw),
    
    kingdom           = str_to_sentence(str_squish(kingdom)),
    taxonrank         = str_to_sentence(str_squish(taxonRank)),
    order_name        = str_squish(Order),
    occurrence_status = str_to_sentence(str_squish(occurrenceStatus)),
    basis_of_record   = str_to_sentence(str_squish(basisOfRecord)),
    
    # Hillside / region (same thing, N vs S)
    region = case_when(
      str_starts(event_id, "N") ~ "North",
      str_starts(event_id, "S") ~ "South",
      TRUE ~ NA_character_
    )
  )

# ---- Genus → synthetic species imputation ----

# Helper: identify if scientific name looks like a proper binomial
is_binomial <- function(x) {
  # Returns TRUE when a string looks like a species binomial and FALSE for
  # genus-only / "sp." style labels (keeps later logic clean and readable).
  x <- str_trim(x)
  pat <- "^[A-Z][a-z]+\\s+[a-z][a-z\\-]+$"
  good <- str_detect(x, pat)
  bad  <- str_detect(x, regex("\\b(sp|spp|species|cf|aff|indet|unknown)\\b", ignore_case = TRUE))
  good & !bad
}

# Work with Title-case genus for imputation
assoc <- assoc %>%
  mutate(
    sci_clean = case_when(
      # If we already have something that looks like a binomial, capitalise nicely
      !is.na(scientific_raw) & scientific_raw != "" &
        is_binomial(str_to_title(scientific_raw)) ~ str_to_title(scientific_raw),
      TRUE ~ str_to_title(scientific_raw)
    ),
    genus      = word(sci_clean, 1),
    epithet    = word(sci_clean, 2),
    epithet_lc = tolower(coalesce(epithet, "")),
    missing_epithet = is.na(epithet) |
      epithet_lc == "" |
      str_detect(epithet_lc, regex("^(sp|spp|sp\\.|spp\\.|cf\\.?|aff\\.?|indet\\.?|unknown)$",
                                   ignore_case = TRUE)),
    
    # taxon_unit:
    # - proper binomials → use Genus species (unchanged)
    # - genus-only / sp. / spp. → assign one synthetic species per genus
    taxon_unit = case_when(
      !is.na(genus) & !missing_epithet ~ str_c(genus, " ", epithet),
      !is.na(genus) & missing_epithet ~ str_c(genus, " species_", tolower(gsub("[^A-Za-z0-9]+", "", genus))),
      TRUE ~ NA_character_
    )
  )

# ---- Filter occurrences for analysis ----

assoc_use <- assoc %>%
  filter(
    occurrence_status == "Present",
    !is.na(region)
  )

# ---- Group assignment ----

mammal_orders      <- c("artiodactyla", "carnivora", "rodentia",
                        "lagomorpha", "eulipotyphla", "chiroptera")
aquatic_orders     <- c("trichoptera","ephemeroptera","plecoptera",
                        "odonata","anisoptera","zygoptera")
terrestrial_orders <- c("lepidoptera","diptera","coleoptera","araneae",
                        "hymenoptera","hemiptera","thysanoptera","ixodida",
                        "julida","tipulidae","tipuloidea","homoptera")

# Birds by event pattern (N1B, S2B, NB_I, etc.)
bird_event_pat <- regex("^[NS](\\d+B(\\b|_)|B_)", ignore_case = TRUE)

assoc_use <- assoc_use %>%
  mutate(
    oname_lc = tolower(order_name),
    cname_lc = tolower(coalesce(common_name, "")),
    sname_lc = tolower(coalesce(scientific_name, "")),
    bid      = tolower(coalesce(basis_of_record, "")),
    evid     = tolower(coalesce(event_id, "")),
    
    group = case_when(
      # Plants
      kingdom == "Plantae" | str_detect(evid, "^[ns]\\d?p") ~ "Plants",
      
      # Bats (by order, event code CR, or audiomoth)
      str_detect(oname_lc, "chiroptera") |
        str_detect(evid, "cr") |
        str_detect(bid, "audiomoth") ~ "Bats",
      
      # Birds (by event IDs)
      str_detect(evid, bird_event_pat) ~ "Birds",
      
      # Aquatic invertebrates
      str_detect(evid, "fi") |
        str_detect(oname_lc, str_c(aquatic_orders, collapse = "|")) |
        str_detect(cname_lc,
                   "caddis|mayfly|stonefly|dragonfly|damselfly|boatman|diving beetle|dytisc|corix") ~
        "Aquatic invertebrates",
      
      # Mammals (non-bats)
      str_detect(oname_lc, str_c(setdiff(mammal_orders, "chiroptera"), collapse = "|")) ~ "Mammals",
      
      # Terrestrial invertebrates
      str_detect(oname_lc, str_c(terrestrial_orders, collapse = "|")) |
        str_detect(evid, "tin|tip|tit") ~ "Terrestrial invertebrates",
      
      TRUE ~ "Other"
    ),
    
    hillside = region  # alias, to match your original plotting terminology
  )

# Only keep rows with usable taxon_unit
oc_use <- assoc_use %>%
  filter(
    !is.na(taxon_unit),
    taxon_unit != "",
    !str_detect(taxon_unit, regex("^No records|Unidentified", ignore_case = TRUE))
  )

# ---- All-taxa plots (not important-only) ----

group_levels <- c("Overall","Plants","Birds","Mammals","Bats",
                  "Terrestrial invertebrates","Aquatic invertebrates")

# 6.1 Richness by group (+ overall), North vs South
rich_by_group <- oc_use %>%
  group_by(group, hillside) %>%
  summarise(n_species = n_distinct(taxon_unit), .groups = "drop")

rich_overall <- oc_use %>%
  group_by(hillside) %>%
  summarise(n_species = n_distinct(taxon_unit), .groups = "drop") %>%
  mutate(group = "Overall")

plot_dat <- bind_rows(rich_by_group, rich_overall) %>%
  filter(group %in% group_levels) %>%
  mutate(
    group    = factor(group, levels = group_levels),
    hillside = factor(hillside, levels = c("North","South"))
  ) %>%
  arrange(group, hillside)

# Prepare a wide summary for the diverging plot
summ_diff <- plot_dat %>%
  select(group, hillside, n_species) %>%
  tidyr::complete(group, hillside, fill = list(n_species = 0)) %>%
  tidyr::pivot_wider(names_from = hillside, values_from = n_species,
                     values_fill = 0) %>%
  mutate(group = factor(group, levels = group_levels))

hillside_dodge <- position_dodge(width = 0.75)

p_rich_all <- ggplot(plot_dat, aes(x = n_species, y = group, fill = hillside)) +
  geom_col(position = hillside_dodge, width = 0.7, colour = "white") +
  geom_text(aes(label = n_species),
            position = hillside_dodge,
            hjust = -0.2, size = 3.8, colour = "grey20") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title    = "Richness by sampled group and hillside",
    subtitle = "Distinct species per hillside (genus-only records imputed to synthetic species)",
    x = "Richness (species)",
    y = "Group",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "top",
    panel.grid.major.y = element_blank(),
    
    # AXIS LINES + TICKS
    axis.line.x  = element_line(colour = "black", linewidth = 0.6),
    axis.line.y  = element_line(colour = "black", linewidth = 0.6),
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_line(colour = "black")
  )


p_rich_all


# 6.2 Diverging richness plot (South − North) by group
df_div <- summ_diff %>%
  mutate(
    diff = South - North,
    abs_diff = abs(diff),
    direction = ifelse(diff > 0, "South",
                       ifelse(diff < 0, "North", "Tie")),
    group = fct_reorder(group, diff),
    group = fct_relevel(group, "Overall", after = Inf)
  )

# ensure ±6 appear on axis
auto_breaks <- scales::pretty_breaks(5)(range(df_div$diff))
forced_breaks <- sort(unique(c(auto_breaks, -6, 6)))

p_diverge_final <- ggplot(df_div, aes(y = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_col(aes(x = diff, fill = direction),
           width = 0.65, colour = "white") +
  geom_text(
    data = df_div %>% filter(diff < 0),
    aes(x = diff - 0.25, label = paste0("N = ", North)),
    hjust = 1, size = 3.6, colour = "grey20"
  ) +
  geom_text(
    data = df_div %>% filter(diff > 0),
    aes(x = diff + 0.25, label = paste0("S = ", South)),
    hjust = 0, size = 3.6, colour = "grey20"
  ) +
  geom_text(
    data = df_div %>% filter(diff < 0),
    aes(x = 0.1, label = paste0("S = ", South)),
    hjust = 0, size = 3.4, colour = "grey40"
  ) +
  geom_text(
    data = df_div %>% filter(diff > 0),
    aes(x = -0.1, label = paste0("N = ", North)),
    hjust = 1, size = 3.4, colour = "grey40"
  ) +
  scale_fill_manual(values = c(
    "North" = "#74A9CF",
    "South" = "#A1D99B",
    "Tie"   = "grey80"
  )) +
  scale_x_continuous(
    labels = abs,
    breaks = forced_breaks,
    expand = expansion(mult = c(0.12, 0.12))
  ) +
  labs(
    title = "Richness difference by group (South minus North)",
    subtitle = "Positive bars = richer in South; negative bars = richer in North",
    x = "Difference in richness (species)",
    y = "Group",
    fill = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, colour = "grey25"),
    axis.line.x = element_line(colour = "grey40", linewidth = 0.6),
    axis.ticks.x = element_line(colour = "grey40"),
    axis.line.y = element_line(colour = "grey40", linewidth = 0.6),
    axis.ticks.y = element_line(colour = "grey40")
  )

p_diverge_final

# 6.3 Rank–abundance across all species, per hillside
rank_abund_all <- oc_use %>%
  count(hillside, taxon_unit, name = "records") %>%
  group_by(hillside) %>%
  arrange(desc(records), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

p_rank_all <- ggplot(rank_abund_all, aes(x = rank, y = records, colour = hillside)) +
  geom_line() +
  geom_point(size = 1.5) +
  scale_y_log10() +
  labs(
    title    = "Rank–abundance (by record frequency)",
    subtitle = "Species-level (genus-only records imputed)",
    x        = "Species rank (per hillside)",
    y        = "Records (log scale)",
    colour   = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

p_rank_all

# 6.4 Similarity (Jaccard) between hillsides by group
jaccard_by_group <- oc_use %>%
  distinct(group, hillside, taxon_unit) %>%
  group_by(group) %>%
  summarise(
    north  = n_distinct(taxon_unit[hillside == "North"]),
    south  = n_distinct(taxon_unit[hillside == "South"]),
    shared = n_distinct(intersect(
      taxon_unit[hillside == "North"],
      taxon_unit[hillside == "South"]
    )),
    jaccard = if_else(north + south - shared == 0,
                      NA_real_, shared / (north + south - shared)),
    .groups = "drop"
  ) %>%
  filter(group %in% group_levels)

p_jaccard_all <- ggplot(jaccard_by_group,
                        aes(x = jaccard,
                            y = factor(group, levels = rev(group_levels)))) +
  geom_vline(xintercept = 0.5, linetype = 3, colour = "grey70") +
  geom_point(size = 3, colour = "#2c7fb8") +
  scale_x_continuous(limits = c(0,1), labels = percent_format(accuracy = 1)) +
  labs(
    title    = "Similarity between hillsides by group",
    subtitle = "Jaccard index over species (genus-only records imputed)",
    x        = "Similarity",
    y        = NULL
  ) +
  theme_minimal(base_size = 12)

p_jaccard_all

# 6.5 Per-group abundance plots (helper function, all taxa)

grouped_species_abundance <- function(df, grp_label, top_n = Inf) {
  
  if (!"group" %in% names(df)) {
    stop("Error: dataset does not contain a 'group' column.")
  }
  if (!grp_label %in% unique(df$group)) {
    warning(paste0("Group '", grp_label, "' not found. Available groups: ",
                   paste(unique(df$group), collapse = ", ")))
    return(invisible(NULL))
  }
  if (!"taxon_unit" %in% names(df)) {
    stop("Error: dataset does not contain 'taxon_unit'.")
  }
  if (!"hillside" %in% names(df)) {
    stop("Error: dataset does not contain 'hillside'.")
  }
  
  df_grp <- df %>% filter(group == grp_label)
  if (nrow(df_grp) == 0) return(invisible(NULL))
  
  dat0 <- df_grp %>%
    count(hillside, taxon_unit, name = "records")
  
  dat <- dat0 %>%
    tidyr::complete(taxon_unit, hillside = c("North","South"),
                    fill = list(records = 0)) %>%
    mutate(hillside = factor(hillside, levels = c("North","South")))
  
  totals <- dat %>%
    group_by(taxon_unit) %>%
    summarise(total_records = sum(records), .groups = "drop")
  
  dat <- dat %>% left_join(totals, by = "taxon_unit")
  
  if (is.finite(top_n)) {
    keep <- totals %>%
      arrange(desc(total_records)) %>%
      slice_head(n = top_n) %>%
      pull(taxon_unit)
    dat <- dat %>% filter(taxon_unit %in% keep)
  }
  
  dat <- dat %>%
    mutate(taxon_unit = reorder(taxon_unit, total_records))
  
  if (!exists("hillside_dodge")) {
    hillside_dodge <- position_dodge(width = 0.75)
  }
  
  ggplot(dat, aes(x = taxon_unit, y = records, fill = hillside)) +
    geom_col(position = hillside_dodge, width = 0.7, colour = "white") +
    geom_text(aes(label = records),
              position = hillside_dodge,
              hjust = -0.1, size = 3, colour = "grey25") +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      title = paste0(grp_label, " — records by hillside (species)"),
      x     = NULL,
      y     = "Record count",
      fill  = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      axis.line.x = element_line(colour = "black", linewidth = 0.6),
      axis.line.y = element_line(colour = "black", linewidth = 0.6),
      axis.ticks.x = element_line(colour = "black"),
      axis.ticks.y = element_line(colour = "black")
    )
}


# Example call per group (comment/uncomment as needed)
# for (g in c("Birds","Mammals","Bats","Terrestrial invertebrates","Aquatic invertebrates","Plants")) {
#   if (g %in% unique(oc_use$group)) print(grouped_species_abundance(oc_use, g))
# }

# ---- Birds: BoCC5 lists ----

bocc5_part1 <- tribble(
  ~common_name, ~scientific_name, ~bocc_status, ~order, ~family,
  
  # RED LIST
  "house sparrow", "passer domesticus", "Red", "Passeriformes", "Passeridae",
  "tree sparrow", "passer montanus", "Red", "Passeriformes", "Passeridae",
  "skylark", "alauda arvensis", "Red", "Passeriformes", "Alaudidae",
  "yellow wagtail", "motacilla flava", "Red", "Passeriformes", "Motacillidae",
  "starling", "sturnus vulgaris", "Red", "Passeriformes", "Sturnidae",
  "mistle thrush", "turdus viscivorus", "Red", "Passeriformes", "Turdidae",
  "song thrush", "turdus philomelos", "Red", "Passeriformes", "Turdidae",
  "fieldfare", "turdus pilaris", "Red", "Passeriformes", "Turdidae",
  "redwing", "turdus iliacus", "Red", "Passeriformes", "Turdidae",
  "spotted flycatcher", "muscicapa striata", "Red", "Passeriformes", "Muscicapidae",
  "pied flycatcher", "ficedula hypoleuca", "Red", "Passeriformes", "Muscicapidae",
  "wood warbler", "phylloscopus sibilatrix", "Red", "Passeriformes", "Phylloscopidae",
  "willow tit", "poecile montanus", "Red", "Passeriformes", "Paridae",
  "marsh tit", "poecile palustris", "Red", "Passeriformes", "Paridae",
  "corn bunting", "emberiza calandra", "Red", "Passeriformes", "Emberizidae",
  "yellowhammer", "emberiza citrinella", "Red", "Passeriformes", "Emberizidae",
  "cuckoo", "cuculus canorus", "Red", "Cuculiformes", "Cuculidae",
  "swift", "apus apus", "Red", "Apodiformes", "Apodidae",
  "kestrel", "falco tinnunculus", "Red", "Falconiformes", "Falconidae",
  "merlin", "falco columbarius", "Red", "Falconiformes", "Falconidae",
  "peregrine", "falco peregrinus", "Red", "Falconiformes", "Falconidae",
  "hen harrier", "circus cyaneus", "Red", "Accipitriformes", "Accipitridae",
  "golden eagle", "aquila chrysaetos", "Red", "Accipitriformes", "Accipitridae",
  "osprey", "pandion haliaetus", "Red", "Accipitriformes", "Pandionidae",
  
  # AMBER LIST BEGIN
  "meadow pipit", "anthus pratensis", "Amber", "Passeriformes", "Motacillidae",
  "rock pipit", "anthus petrosus", "Amber", "Passeriformes", "Motacillidae",
  "grey wagtail", "motacilla cinerea", "Amber", "Passeriformes", "Motacillidae",
  "robin", "erithacus rubecula", "Amber", "Passeriformes", "Muscicapidae",
  "dunnock", "prunella modularis", "Amber", "Passeriformes", "Prunellidae",
  "blackbird", "turdus merula", "Amber", "Passeriformes", "Turdidae",
  "ring ouzel", "turdus torquatus", "Amber", "Passeriformes", "Turdidae",
  "whinchat", "saxicola rubetra", "Amber", "Passeriformes", "Muscicapidae",
  "stonechat", "saxicola rubicola", "Amber", "Passeriformes", "Muscicapidae",
  "wheatear", "oenanthe oenanthe", "Amber", "Passeriformes", "Muscicapidae",
  "long-tailed tit", "aegithalos caudatus", "Amber", "Passeriformes", "Aegithalidae",
  "great tit", "parus major", "Amber", "Passeriformes", "Paridae",
  "blue tit", "cyanistes caeruleus", "Amber", "Passeriformes", "Paridae",
  "coal tit", "periparus ater", "Amber", "Passeriformes", "Paridae",
  "goldcrest", "regulus regulus", "Amber", "Passeriformes", "Regulidae",
  "firecrest", "regulus ignicapilla", "Amber", "Passeriformes", "Regulidae",
  "bullfinch", "pyrrhula pyrrhula", "Amber", "Passeriformes", "Fringillidae",
  "linnet", "linaria cannabina", "Amber", "Passeriformes", "Fringillidae",
  "lesser redpoll", "carduelis cabaret", "Amber", "Passeriformes", "Fringillidae",
  "siskin", "spinus spinus", "Amber", "Passeriformes", "Fringillidae",
  "goldfinch", "carduelis carduelis", "Amber", "Passeriformes", "Fringillidae",
  "greenfinch", "chloris chloris", "Amber", "Passeriformes", "Fringillidae",
  "chaffinch", "fringilla coelebs", "Amber", "Passeriformes", "Fringillidae",
  "brambling", "fringilla montifringilla", "Amber", "Passeriformes", "Fringillidae",
  
  # Continue Amber species
  "jackdaw", "corvus monedula", "Amber", "Passeriformes", "Corvidae",
  "rook", "corvus frugilegus", "Amber", "Passeriformes", "Corvidae",
  "carrion crow", "corvus corone", "Amber", "Passeriformes", "Corvidae",
  "raven", "corvus corax", "Amber", "Passeriformes", "Corvidae",
  "magpie", "pica pica", "Amber", "Passeriformes", "Corvidae",
  "jay", "garrulus glandarius", "Amber", "Passeriformes", "Corvidae",
  
  "wood pigeon", "columba palumbus", "Amber", "Columbiformes", "Columbidae",
  "stock dove", "columba oenas", "Amber", "Columbiformes", "Columbidae",
  "collared dove", "streptopelia decaocto", "Amber", "Columbiformes", "Columbidae",
  
  # WADERS
  "lapwing", "vanellus vanellus", "Amber", "Charadriiformes", "Charadriidae",
  "oystercatcher", "haematopus ostralegus", "Amber", "Charadriiformes", "Haematopodidae",
  "curlew", "numenius arquata", "Amber", "Charadriiformes", "Scolopacidae",
  "snipe", "gallinago gallinago", "Amber", "Charadriiformes", "Scolopacidae",
  "woodcock", "scolopax rusticola", "Amber", "Charadriiformes", "Scolopacidae",
  
  # GULLS
  "herring gull", "larus argentatus", "Amber", "Charadriiformes", "Laridae",
  "kittiwake", "rissa tridactyla", "Amber", "Charadriiformes", "Laridae",
  "common gull", "larus canus", "Amber", "Charadriiformes", "Laridae"
)

bocc5_part2 <- tribble(
  ~common_name, ~scientific_name, ~bocc_status, ~order, ~family,
  
  # Continue Amber-listed WADERS & GULLS
  "black-headed gull", "chroicocephalus ridibundus", "Amber", "Charadriiformes", "Laridae",
  "lesser black-backed gull", "larus fuscus", "Amber", "Charadriiformes", "Laridae",
  "great black-backed gull", "larus marinus", "Amber", "Charadriiformes", "Laridae",
  "sandwich tern", "thalasseus sandvicensis", "Amber", "Charadriiformes", "Laridae",
  "common tern", "sterna hirundo", "Amber", "Charadriiformes", "Laridae",
  "arctic tern", "sterna paradisaea", "Amber", "Charadriiformes", "Laridae",
  
  # WATERFOWL
  "mallard", "anas platyrhynchos", "Amber", "Anseriformes", "Anatidae",
  "teal", "anas crecca", "Amber", "Anseriformes", "Anatidae",
  "wigeon", "mareca penelope", "Amber", "Anseriformes", "Anatidae",
  "pintail", "anas acuta", "Amber", "Anseriformes", "Anatidae",
  "shoveler", "spatula clypeata", "Amber", "Anseriformes", "Anatidae",
  "pochard", "aythya ferina", "Amber", "Anseriformes", "Anatidae",
  "tufted duck", "aythya fuligula", "Amber", "Anseriformes", "Anatidae",
  "goldeneye", "bucephala clangula", "Amber", "Anseriformes", "Anatidae",
  "goosander", "mergus merganser", "Amber", "Anseriformes", "Anatidae",
  
  # GEESE
  "pink-footed goose", "anser brachyrhynchus", "Amber", "Anseriformes", "Anatidae",
  "greylag goose", "anser anser", "Amber", "Anseriformes", "Anatidae",
  "canada goose", "branta canadensis", "Amber", "Anseriformes", "Anatidae",
  "barnacle goose", "branta leucopsis", "Amber", "Anseriformes", "Anatidae",
  
  # RAPTORS continued
  "buzzard", "buteo buteo", "Amber", "Accipitriformes", "Accipitridae",
  "red kite", "milvus milvus", "Amber", "Accipitriformes", "Accipitridae",
  
  # OWLS
  "tawny owl", "strix aluco", "Amber", "Strigiformes", "Strigidae",
  "barn owl", "tyto alba", "Amber", "Strigiformes", "Tytonidae",
  "long-eared owl", "asio otus", "Amber", "Strigiformes", "Strigidae",
  "short-eared owl", "asio flammeus", "Amber", "Strigiformes", "Strigidae",
  
  # WRENS
  "wren", "troglodytes troglodytes", "Amber", "Passeriformes", "Troglodytidae",
  
  # WOODPECKERS
  "great spotted woodpecker", "dendrocopos major", "Amber", "Piciformes", "Picidae",
  "green woodpecker", "picus viridis", "Amber", "Piciformes", "Picidae",
  "lesser spotted woodpecker", "dryobates minor", "Amber", "Piciformes", "Picidae",
  
  # GAMEBIRDS
  "red grouse", "lagopus lagopus scotica", "Amber", "Galliformes", "Phasianidae",
  "pheasant", "phasianus colchicus", "Amber", "Galliformes", "Phasianidae",
  "grey partridge", "perdix perdix", "Amber", "Galliformes", "Phasianidae",
  
  # HERONS / EGRETS
  "grey heron", "ardea cinerea", "Amber", "Pelecaniformes", "Ardeidae",
  "little egret", "egretta garzetta", "Amber", "Pelecaniformes", "Ardeidae",
  
  # GREBES
  "great crested grebe", "podiceps cristatus", "Amber", "Podicipediformes", "Podicipedidae",
  "little grebe", "tachybaptus ruficollis", "Amber", "Podicipediformes", "Podicipedidae",
  
  # SWIFTS & SWALLOWS
  "swallow", "hirundo rustica", "Amber", "Passeriformes", "Hirundinidae",
  "house martin", "delichon urbicum", "Amber", "Passeriformes", "Hirundinidae",
  "sand martin", "riparia riparia", "Amber", "Passeriformes", "Hirundinidae",
  
  # WARBLERS
  "chiffchaff", "phylloscopus collybita", "Amber", "Passeriformes", "Phylloscopidae",
  "willow warbler", "phylloscopus trochilus", "Amber", "Passeriformes", "Phylloscopidae",
  "blackcap", "sylvia atricapilla", "Amber", "Passeriformes", "Sylviidae",
  "garden warbler", "sylvia borin", "Amber", "Passeriformes", "Sylviidae",
  "whitethroat", "sylvia communis", "Amber", "Passeriformes", "Sylviidae",
  "lesser whitethroat", "curruca curruca", "Amber", "Passeriformes", "Sylviidae",
  
  # TREECREEPERS
  "treecreeper", "certhia familiaris", "Amber", "Passeriformes", "Certhiidae",
  
  # KINGFISHER
  "kingfisher", "alcedo atthis", "Amber", "Coraciiformes", "Alcedinidae",
  
  # DIPPER
  "dipper", "cinclus cinclus", "Amber", "Passeriformes", "Cinclidae",
  
  # GREEN LIST (common species)
  "robin", "erithacus rubecula", "Green", "Passeriformes", "Muscicapidae",
  "blue tit", "cyanistes caeruleus", "Green", "Passeriformes", "Paridae",
  "great tit", "parus major", "Green", "Passeriformes", "Paridae",
  "coal tit", "periparus ater", "Green", "Passeriformes", "Paridae",
  "long-tailed tit", "aegithalos caudatus", "Green", "Passeriformes", "Aegithalidae",
  "chaffinch", "fringilla coelebs", "Green", "Passeriformes", "Fringillidae",
  "goldfinch", "carduelis carduelis", "Green", "Passeriformes", "Fringillidae",
  "siskin", "spinus spinus", "Green", "Passeriformes", "Fringillidae",
  "greenfinch", "chloris chloris", "Green", "Passeriformes", "Fringillidae",
  "wren", "troglodytes troglodytes", "Green", "Passeriformes", "Troglodytidae",
  "blackbird", "turdus merula", "Green", "Passeriformes", "Turdidae",
  "dunnock", "prunella modularis", "Green", "Passeriformes", "Prunellidae",
  "wood pigeon", "columba palumbus", "Green", "Columbiformes", "Columbidae",
  "feral pigeon", "columba livia domestica", "Green", "Columbiformes", "Columbidae",
  "magpie", "pica pica", "Green", "Passeriformes", "Corvidae",
  "jackdaw", "corvus monedula", "Green", "Passeriformes", "Corvidae",
  "carrion crow", "corvus corone", "Green", "Passeriformes", "Corvidae",
  "raven", "corvus corax", "Green", "Passeriformes", "Corvidae",
  "collared dove", "streptopelia decaocto", "Green", "Columbiformes", "Columbidae",
  "great spotted woodpecker", "dendrocopos major", "Green", "Piciformes", "Picidae"
)

bocc5_part3 <- tribble(
  ~common_name, ~scientific_name, ~bocc_status, ~order, ~family,
  
  # Continue GREEN LIST
  "goldcrest", "regulus regulus", "Green", "Passeriformes", "Regulidae",
  "pheasant", "phasianus colchicus", "Green", "Galliformes", "Phasianidae",
  "moorhen", "gallinula chloropus", "Green", "Gruiformes", "Rallidae",
  "coot", "fulica atra", "Green", "Gruiformes", "Rallidae",
  "mallard", "anas platyrhynchos", "Green", "Anseriformes", "Anatidae",
  "tufted duck", "aythya fuligula", "Green", "Anseriformes", "Anatidae",
  
  # COMMON WATERBIRDS
  "mute swan", "cygnus olor", "Green", "Anseriformes", "Anatidae",
  "grey heron", "ardea cinerea", "Green", "Pelecaniformes", "Ardeidae",
  "little egret", "egretta garzetta", "Green", "Pelecaniformes", "Ardeidae",
  
  # GREBES
  "great crested grebe", "podiceps cristatus", "Green", "Podicipediformes", "Podicipedidae",
  
  # COMMON WADERS
  "redshank", "tringa totanus", "Green", "Charadriiformes", "Scolopacidae",
  "common sandpiper", "actitis hypoleucos", "Green", "Charadriiformes", "Scolopacidae",
  
  # COMMON SEA BIRDS
  "gannet", "morus bassanus", "Green", "Suliformes", "Sulidae",
  "great black-backed gull", "larus marinus", "Green", "Charadriiformes", "Laridae",
  "lesser black-backed gull", "larus fuscus", "Green", "Charadriiformes", "Laridae",
  "black-headed gull", "chroicocephalus ridibundus", "Green", "Charadriiformes", "Laridae",
  
  # COMMON COASTAL BIRDS
  "oystercatcher", "haematopus ostralegus", "Green", "Charadriiformes", "Haematopodidae",
  "curlew", "numenius arquata", "Green", "Charadriiformes", "Scolopacidae",
  
  # FAMILIAR SPECIES
  "robin", "erithacus rubecula", "Green", "Passeriformes", "Muscicapidae",
  "wren", "troglodytes troglodytes", "Green", "Passeriformes", "Troglodytidae",
  "great tit", "parus major", "Green", "Passeriformes", "Paridae",
  "chaffinch", "fringilla coelebs", "Green", "Passeriformes", "Fringillidae",
  "goldfinch", "carduelis carduelis", "Green", "Passeriformes", "Fringillidae",
  
  # WILD PREDATORS
  "kestrel", "falco tinnunculus", "Green", "Falconiformes", "Falconidae",
  "peregrine falcon", "falco peregrinus", "Green", "Falconiformes", "Falconidae",
  
  # PIGEONS & DOVES
  "wood pigeon", "columba palumbus", "Green", "Columbiformes", "Columbidae",
  "feral pigeon", "columba livia domestica", "Green", "Columbiformes", "Columbidae",
  "collared dove", "streptopelia decaocto", "Green", "Columbiformes", "Columbidae",
  
  # ADDITIONAL AMBER (not previously listed)
  "grey wagtail", "motacilla cinerea", "Amber", "Passeriformes", "Motacillidae",
  "reed bunting", "emberiza schoeniclus", "Amber", "Passeriformes", "Emberizidae",
  "sedge warbler", "acrocephalus schoenobaenus", "Amber", "Passeriformes", "Acrocephalidae",
  "grasshopper warbler", "locustella naevia", "Amber", "Passeriformes", "Locustellidae",
  
  # EXTRA RED (not previously included but common in UK lists)
  "grey partridge", "perdix perdix", "Red", "Galliformes", "Phasianidae",
  "ring ouzel", "turdus torquatus", "Red", "Passeriformes", "Turdidae",
  "yellowhammer", "emberiza citrinella", "Red", "Passeriformes", "Emberizidae",
  
  # DUPLICATES REMOVED AUTOMATICALLY LATER
  "house sparrow", "passer domesticus", "Green", "Passeriformes", "Passeridae",
  
  # END
  NA, NA, NA, NA, NA
) %>%
  filter(!is.na(common_name))

bocc5 <- bind_rows(bocc5_part1, bocc5_part2, bocc5_part3) %>%
  distinct(common_name, .keep_all = TRUE)

# Fix obvious BoCC5 mislabels (BoCC5 2021) and harmonise
bocc5 <- bocc5 %>%
  mutate(
    common_name     = tolower(common_name),
    scientific_name = tolower(scientific_name),
    bocc_status     = tolower(bocc_status)
  ) %>%
  mutate(
    bocc_status = case_when(
      # Raptors/owls (headline fixes)
      common_name %in% c("kestrel") ~ "amber",
      common_name %in% c("golden eagle") ~ "amber",
      common_name %in% c("peregrine","peregrine falcon") ~ "green",
      common_name %in% c("buzzard") ~ "green",
      common_name %in% c("red kite") ~ "green",
      common_name %in% c("tawny owl") ~ "green",
      
      # Waders/gulls/terns/waterbirds
      common_name %in% c("lapwing") ~ "red",
      common_name %in% c("curlew") ~ "red",
      common_name %in% c("redshank") ~ "amber",
      common_name %in% c("snipe") ~ "amber",
      common_name %in% c("woodcock") ~ "red",
      common_name %in% c("herring gull") ~ "red",
      common_name %in% c("kittiwake") ~ "red",
      common_name %in% c("common gull") ~ "amber",
      common_name %in% c("black-headed gull") ~ "amber",
      common_name %in% c("lesser black-backed gull") ~ "amber",
      common_name %in% c("great black-backed gull") ~ "amber",
      common_name %in% c("sandwich tern","common tern","arctic tern") ~ "amber",
      common_name %in% c("pochard") ~ "red",
      common_name %in% c("mallard","tufted duck","goldeneye","goosander",
                         "teal","wigeon","pintail","shoveler") ~ "amber",
      # Geese (season/local dependence; default split)
      common_name %in% c("pink-footed goose","barnacle goose") ~ "amber",
      common_name %in% c("greylag goose","canada goose") ~ "green",
      
      # Common passerines to Green (not usually "important" in EcIA)
      common_name %in% c("robin","blackbird","blue tit","great tit","coal tit","long-tailed tit",
                         "wren","magpie","jackdaw","carrion crow","raven","jay",
                         "wood pigeon","collared dove","great spotted woodpecker",
                         "green woodpecker","goldcrest","firecrest","siskin","goldfinch",
                         "treecreeper","stonechat","chiffchaff","blackcap","garden warbler",
                         "whitethroat","gannet","oystercatcher","grey heron","little egret",
                         "great crested grebe","little grebe","swallow","sand martin") ~ "green",
      
      # Amber/Red updates
      common_name %in% c("dunnock","wheatear","willow warbler","lesser whitethroat",
                         "kingfisher","dipper","rock pipit","meadow pipit",
                         "grey wagtail","bullfinch","linnet","brambling","stock dove") ~ "amber",
      common_name %in% c("lesser redpoll","rook","ring ouzel","house martin",
                         "greenfinch","chaffinch") ~ "red",
      
      TRUE ~ bocc_status
    )
  ) %>%
  mutate(bocc_status = tools::toTitleCase(bocc_status))

# Keep only Red/Amber birds as "important" for EcIA
bocc5_for_importance <- bocc5 %>% filter(bocc_status %in% c("Red","Amber"))

# ---- Mammals: protected species list ----

mammals_imp <- tribble(
  ~common_name,            ~scientific_name,          ~status,                    ~order,         ~family,
  # EPS (Habitats Regulations)
  "otter",                 "lutra lutra",             "EPS (fully protected)",    "Carnivora",    "Mustelidae",
  "common pipistrelle",    "pipistrellus pipistrellus","EPS (fully protected)",   "Chiroptera",   "Vespertilionidae",
  "soprano pipistrelle",   "pipistrellus pygmaeus",   "EPS (fully protected)",    "Chiroptera",   "Vespertilionidae",
  "nathusius' pipistrelle","pipistrellus nathusii",   "EPS (fully protected)",    "Chiroptera",   "Vespertilionidae",
  "lesser noctule",        "nyctalus leisleri",       "EPS (fully protected)",    "Chiroptera",   "Vespertilionidae",
  "brown long-eared bat",  "plecotus auritus",        "EPS (fully protected)",    "Chiroptera",   "Vespertilionidae",
  
  # WCA S5 / Protected (Scotland)
  "badger",                "meles meles",             "Protected (PBA 1992)",     "Carnivora",    "Mustelidae",
  "red squirrel",          "sciurus vulgaris",        "Protected (WCA S5), SBL",  "Rodentia",     "Sciuridae",
  "pine marten",           "martes martes",           "Protected (WCA S5)",       "Carnivora",    "Mustelidae",
  "water vole",            "arvicola amphibius",      "Protected (WCA S5)",       "Rodentia",     "Cricetidae",
  "wildcat",               "felis silvestris",        "EPS (fully protected)",    "Carnivora",    "Felidae",
  "mountain hare",         "lepus timidus",           "Protected (Scotland)",     "Lagomorpha",   "Leporidae",
  "beaver",                "castor fiber",            "EPS (fully protected)",    "Rodentia",     "Castoridae",
  
  # SBL Priority (often material in EcIA)
  "hedgehog",              "erinaceus europaeus",     "SBL Priority",             "Eulipotyphla", "Erinaceidae",
  "brown hare",            "lepus europaeus",         "SBL Priority",             "Lagomorpha",   "Leporidae"
) %>%
  mutate(
    common_name = tolower(common_name),
    scientific_name = tolower(scientific_name),
    bocc_status = NA_character_,
    group = "Mammal"
  )

# Optional: context mammals (not "protected species")
# other_mammals_context <- tribble(
#   ~common_name, ~scientific_name,      ~status,                              ~order,         ~family,
#   "red deer",   "cervus elaphus",      "Managed (not protected species)",    "Artiodactyla", "Cervidae",
#   "roe deer",   "capreolus capreolus", "Managed (not protected species)",    "Artiodactyla", "Cervidae",
#   "red fox",    "vulpes vulpes",       "Common (no special protection)",     "Carnivora",    "Canidae"
# ) %>% mutate(
#   common_name = tolower(common_name),
#   scientific_name = tolower(scientific_name),
#   bocc_status = NA_character_,
#   group = "Mammal"
# )

# ---- Build important-species lookup ----

lookup_all <- bind_rows(
  bocc5_for_importance %>%
    mutate(group = "Bird", status = bocc_status) %>%
    select(common_name, scientific_name, status, order, family, group),
  mammals_imp %>%
    select(common_name, scientific_name, status, order, family, group)
  # , other_mammals_context   # optionally include as non-important context
) %>%
  mutate(
    common_name     = tolower(common_name),
    scientific_name = tolower(scientific_name)
  ) %>%
  distinct(common_name, scientific_name, .keep_all = TRUE)

# ---- Normalise observed names ----

assoc <- assoc %>%
  mutate(
    common_name     = tolower(str_squish(common_name)),
    scientific_name = tolower(str_squish(scientific_name))
  )

synonyms <- tribble(
  ~common_name_obs,        ~common_name_std,
  "peregrine",             "peregrine falcon",
  "nathusius pipistrelle", "nathusius' pipistrelle"
)

assoc <- assoc %>%
  left_join(synonyms, by = c("common_name" = "common_name_obs")) %>%
  mutate(common_name = coalesce(common_name_std, common_name)) %>%
  select(-common_name_std)

# ---- Fuzzy matching ----

# Filter out NA names so fuzzyjoin doesn't choke on them
assoc_clean  <- assoc %>% dplyr::filter(!is.na(common_name))
lookup_clean <- lookup_all %>% dplyr::filter(!is.na(common_name))

matched_common <- stringdist_inner_join(
  assoc_clean,
  lookup_clean,
  by = c("common_name" = "common_name"),
  max_dist     = 2,
  distance_col = "dist_common"
)

assoc_cleans  <- assoc %>% dplyr::filter(!is.na(scientific_name))
lookup_cleans <- lookup_all %>% dplyr::filter(!is.na(scientific_name))

matched_sci <- stringdist_inner_join(
  assoc_cleans,
  lookup_cleans,
  by = c("scientific_name" = "scientific_name"),
  max_dist     = 2,
  distance_col = "dist_sci"
)

matched_all <- bind_rows(matched_common, matched_sci) %>%
  mutate(key_sci = coalesce(scientific_name.y, scientific_name.x)) %>%
  arrange(event_id, key_sci, dist_common, dist_sci) %>%
  distinct(event_id, key_sci, .keep_all = TRUE) %>%
  rename(
    common_name_obs       = common_name.x,
    scientific_name_obs   = scientific_name.x,
    common_name_ref       = common_name.y,
    scientific_name_ref   = scientific_name.y
  )

# ---- Build important species table ----

important_species_all <- matched_all %>%
  transmute(
    event_id,
    region,
    group,                                   # Bird / Mammal
    common_name     = coalesce(common_name_ref, common_name_obs),
    scientific_name = coalesce(scientific_name_ref, scientific_name_obs),
    status          = status,                # Birds: Red/Amber; Mammals: EPS/Protected/SBL
    order,
    family
  )

# Optional quick QA
# important_species_all %>% filter(group == "Bird", status == "Green") %>% distinct(common_name)  # should be empty
# important_species_all %>% filter(group == "Mammal") %>% count(status, sort = TRUE)
# assoc %>% anti_join(lookup_all, by = c("common_name","scientific_name")) %>% distinct(common_name, scientific_name)

# ---- Setup: important-only dataset for plots ----

# EcIA-ready plots use only important species (Birds: Red/Amber; Mammals: EPS/Protected/SBL)
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
}

# Exclude out-of-boundary events if event_id encodes this
dat <- important_species_all %>%
  filter(!str_detect(tolower(event_id), "\\boob\\b"))

# Consistent ordering and labels
dat <- dat %>%
  mutate(
    group = factor(group, levels = c("Bird","Mammal")),
    region = factor(region),
    species_label = str_to_title(common_name)  # nicer labels for plots
  )

# Palettes
pal_region <- c("North" = "#E76F51", "South" = "#2A9D8F") %>%
  { .[intersect(names(.), levels(dat$region))] } # keep only present regions
pal_bocc   <- c("Red" = "#C1121F", "Amber" = "#FFB000")

theme_ecia <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      legend.position = "right",
      strip.text = element_text(face = "bold")
    )
}

# ---- Plot: Richness (important only) ----

richness <- dat %>%
  distinct(region, group, scientific_name) %>%
  count(region, group, name = "n_species") %>%
  group_by(region) %>%
  mutate(pct = n_species / sum(n_species)) %>%
  ungroup()

fig_richnessimportant <- ggplot(richness, aes(x = region, y = n_species, fill = group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(aes(label = n_species),
            position = position_dodge(width = 0.75),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Bird" = "#577590", "Mammal" = "#43AA8B")) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08)),
    limits = c(0, max(richness$n_species) * 1.15)
  ) +
  labs(
    title = "Species richness by region and group (important only)",
    x = "Region",
    y = "Number of important species",
    fill = "Group"
  ) +
  theme_ecia() +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.4),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(size = 11, colour = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

fig_richnessimportant


# ---- Plot: Abundance (important only) ----

top_n_per_group <- 10

abund_species <- dat %>%
  count(group, region, species_label, name = "n_records") %>%
  group_by(group, species_label) %>%
  mutate(total = sum(n_records)) %>%
  ungroup() %>%
  group_by(group) %>%
  slice_max(order_by = total, n = top_n_per_group, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(species_label = fct_reorder(species_label, total))

fig_abundanceimportant <- ggplot(abund_species, aes(x = n_records, y = species_label, fill = region)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = n_records),
            position = position_dodge(width = 0.8), hjust = -0.2, size = 3) +
  facet_wrap(~ group, scales = "free_y") +
  scale_fill_manual(values = pal_region) +
  scale_x_continuous(expand = expansion(mult = c(0, .1))) +
  labs(
    title = "Abundance (record count) by species and region (important only)",
    x = "Number of important records",
    y = "Species",
    fill = "Region"
  ) +
  theme_ecia() +
  theme(
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

fig_abundanceimportant

# ---- Plot: Heatmap (important only) ----

# Add species-level conservation category
species_status <- dat %>%
  mutate(category = case_when(
    # Birds: BoCC
    group == "Bird" & status %in% c("Red", "Amber") ~ status,
    group == "Bird"                                 ~ "Other",
    
    # Mammals: protection categories
    group == "Mammal" & str_detect(status, regex("EPS", ignore_case = TRUE)) ~ "EPS",
    group == "Mammal" & str_detect(status, regex("Protected|WCA|PBA", ignore_case = TRUE)) ~ "Protected",
    group == "Mammal" & str_detect(status, regex("SBL", ignore_case = TRUE)) ~ "SBL Priority",
    group == "Mammal" ~ "Other",
    
    TRUE ~ "Other"
  )) %>%
  distinct(group, species_label, category)

# Build heatmap data with conservation category attached
heat <- dat %>%
  count(group, species_label, region, name = "n_records") %>%
  group_by(group, species_label) %>%
  mutate(total = sum(n_records)) %>%
  ungroup() %>%
  left_join(species_status, by = c("group", "species_label")) %>%
  mutate(
    species_label = fct_reorder(species_label, total),
    category = factor(
      category,
      levels = c("SBL Priority", "Protected", "EPS", "Amber", "Red", "Other")
    )
  )

# Unified palette
pal_comp <- c(
  "Amber"        = "#FFB000",
  "Red"          = "#C1121F",
  "EPS"          = "#2D6A4F",
  "Protected"    = "#40916C",
  "SBL Priority" = "#74C69D",
  "Other"        = "grey80"
)

fig_heatmapimportant <- ggplot(heat, aes(x = region, y = species_label)) +
  geom_tile(aes(fill = category),
            color = "grey80",  # slightly darker borders for clarity
            linewidth = 0.35) +
  geom_text(aes(label = n_records),
            color = "black",
            size = 3.6,
            fontface = "bold") +
  facet_wrap(~ group, scales = "free_y", ncol = 1, strip.position = "top") +
  scale_fill_manual(values = pal_comp, name = "Category") +
  labs(
    title = "Species–region heatmap of important records",
    subtitle = "Tile colour indicates conservation/protection category; numbers show important record counts",
    x = "Region",
    y = "Species"
  ) +
  theme_ecia() +
  theme(
    # Axis clarity
    axis.line.x = element_line(colour = "black", linewidth = 0.4),
    axis.line.y = element_line(colour = "black", linewidth = 0.4),
    axis.ticks = element_line(colour = "black"),
    
    # Label clarity
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    
    # Facet spacing and layout
    panel.spacing.y = unit(1.2, "lines"),
    strip.text = element_text(size = 11, face = "bold"),
    
    # Titles
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, margin = margin(b = 10)),
    
    # Legend
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 9),
    
    # Margins
    plot.margin = margin(10, 15, 10, 10)
  )

fig_heatmapimportant


# ---- Plot: Bird BoCC composition (important only) ----

birds_comp_counts <- dat %>%
  filter(group == "Bird", status %in% c("Red","Amber")) %>%
  mutate(status = factor(status, levels = c("Amber","Red"))) %>%  # Amber at bottom, Red on top
  count(region, status, name = "n") %>%
  group_by(region) %>%
  mutate(total = sum(n)) %>%
  ungroup()

fig_bocc_comp_counts <- ggplot(birds_comp_counts, aes(region, n, fill = status)) +
  geom_col(width = 0.7, color = "white") +
  geom_text(
    aes(label = ifelse(n > 0, n, "")),
    position = position_stack(vjust = 0.5), size = 3.2, color = "black"
  ) +
  geom_text(
    data = distinct(birds_comp_counts, region, total),
    aes(region, total, label = paste0("Total: ", total)),
    vjust = -0.5, size = 3.2, fontface = "bold", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = pal_bocc) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(
    title = "Bird BoCC status by region (important only, stacked counts)",
    x = "Region", y = "Number of important bird records", fill = "BoCC status"
  ) +
  theme_ecia()

fig_bocc_comp_counts

# ---- Plot: Mammal protection composition (important only) ----

mammals_comp_counts <- dat %>%
  filter(group == "Mammal") %>%
  mutate(protection = case_when(
    str_detect(status, regex("EPS", ignore_case = TRUE)) ~ "EPS",
    str_detect(status, regex("Protected|WCA|PBA", ignore_case = TRUE)) ~ "Protected",
    str_detect(status, regex("SBL", ignore_case = TRUE)) ~ "SBL Priority",
    TRUE ~ status
  )) %>%
  mutate(protection = factor(protection, levels = c("SBL Priority","Protected","EPS"))) %>%
  count(region, protection, name = "n") %>%
  group_by(region) %>%
  mutate(total = sum(n)) %>%
  ungroup()

pal_protect <- c("EPS" = "#2D6A4F", "Protected" = "#40916C", "SBL Priority" = "#74C69D")

fig_mammal_comp_counts <- ggplot(mammals_comp_counts, aes(region, n, fill = protection)) +
  geom_col(width = 0.7, color = "white") +
  geom_text(
    aes(label = ifelse(n > 0, n, "")),
    position = position_stack(vjust = 0.5), size = 3.2, color = "black"
  ) +
  geom_text(
    data = distinct(mammals_comp_counts, region, total),
    aes(region, total, label = paste0("Total: ", total)),
    vjust = -0.5, size = 3.2, fontface = "bold", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = pal_protect) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(
    title = "Mammal protection category by region (important only, stacked counts)",
    x = "Region", y = "Number of important mammal records", fill = "Category"
  ) +
  theme_ecia()

fig_mammal_comp_counts

# ---- Plot: Combined composition (important only) ----

# Build bird composition (Red/Amber)
birds_comp_counts <- dat %>%
  filter(group == "Bird", status %in% c("Red","Amber")) %>%
  mutate(category = factor(status, levels = c("Amber","Red"))) %>%  # Amber bottom, Red top
  count(group, region, category, name = "n")

# Build mammal composition (EPS / Protected / SBL)
mammals_comp_counts <- dat %>%
  filter(group == "Mammal") %>%
  mutate(category = case_when(
    str_detect(status, regex("EPS", ignore_case = TRUE)) ~ "EPS",
    str_detect(status, regex("Protected|WCA|PBA", ignore_case = TRUE)) ~ "Protected",
    str_detect(status, regex("SBL", ignore_case = TRUE)) ~ "SBL Priority",
    TRUE ~ status
  )) %>%
  mutate(category = factor(category, levels = c("SBL Priority","Protected","EPS"))) %>% # SBL bottom, EPS top
  count(group, region, category, name = "n")

# Combine and order categories across both groups
comp_all <- bind_rows(birds_comp_counts, mammals_comp_counts) %>%
  mutate(
    region = factor(region),
    group = factor(group, levels = c("Bird","Mammal")),
    category = factor(category, levels = c("SBL Priority","Protected","EPS","Amber","Red"))
  )

# Totals per region in each facet
comp_totals <- comp_all %>%
  group_by(group, region) %>%
  summarise(total = sum(n), .groups = "drop")

# Unified palette for all categories
pal_comp <- c(
  "Amber" = "#FFB000",
  "Red" = "#C1121F",
  "EPS" = "#2D6A4F",
  "Protected" = "#40916C",
  "SBL Priority" = "#74C69D"
)

fig_comp_combinedimportantbyregionandstatus <- ggplot(comp_all, aes(region, n, fill = category)) +
  geom_col(
    width = 0.7,
    color = "grey90",        # clearer bar edges
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = ifelse(n > 0, n, "")),
    position = position_stack(vjust = 0.5),
    size = 3.4,
    color = "black",
    fontface = "bold"
  ) +
  geom_text(
    data = comp_totals,
    aes(region, total, label = paste0("Total: ", total)),
    vjust = -0.6,
    size = 3.4,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ group, ncol = 2, scales = "fixed") +
  scale_fill_manual(values = pal_comp, drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Composition of important records by region (Bird BoCC and Mammal protection)",
    x = "Region",
    y = "Number of important records",
    fill = "Category"
  ) +
  theme_ecia() +
  theme(
    # clearer axes
    axis.line = element_line(colour = "black", linewidth = 0.4),
    axis.ticks = element_line(colour = "black"),
    
    # axis text and titles
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    
    # facet labels
    strip.text = element_text(size = 11, face = "bold"),
    
    # title
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    
    # legend
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    
    plot.margin = margin(10, 15, 10, 10)
  )

fig_comp_combinedimportantbyregionandstatus






#----plots for report----
p_rich_all
p_diverge_final
fig_abundanceimportant
fig_richnessimportant
fig_heatmapimportant
fig_comp_combinedimportantbyregionandstatus
