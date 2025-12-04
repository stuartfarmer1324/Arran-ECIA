# Arran biodiversity analysis + important species pipeline
# Uses associatedoccurences.csv as the data source
# Imputes genus-only records to synthetic taxa without altering real species
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
    basisOfRecord  = rep("HumanObservation", 8),
    individualCount = c(3, 5, 1, 2, 10, 8, 4, 6)
  )
}
assoc
# ---- Standardise fields ----

assoc <- assoc %>%
  rename(
    event_id        = eventID,
    common_name_raw = Commonme,
    scientific_raw  = scientificme
  ) %>%
  mutate(
    common_name_raw   = str_squish(common_name_raw),
    scientific_raw    = str_squish(scientific_raw),
    common_name       = str_to_lower(common_name_raw),
    scientific_name   = str_to_lower(scientific_raw),
    # NEW: title-case common name for policy-facing outputs
    common_name_title = str_to_title(common_name_raw),
    
    kingdom           = str_to_sentence(str_squish(kingdom)),
    taxonrank         = str_to_sentence(str_squish(taxonRank)),
    order_name        = str_squish(Order),
    occurrence_status = str_to_sentence(str_squish(occurrenceStatus)),
    basis_of_record   = str_to_sentence(str_squish(basisOfRecord)),
    
    # abundance: make sure it's numeric, treat missing as 0
    individualCount   = replace_na(as.numeric(individualCount), 0),
    
    # Hillside / region (same thing, N vs S)
    region = case_when(
      str_starts(event_id, "N") ~ "North",
      str_starts(event_id, "S") ~ "South",
      TRUE ~ NA_character_
    )
  )


# ---- Genus → synthetic taxon imputation (hybrid sci + common name) ----

# Helper: extract scientific genus if present
extract_genus <- function(x) {
  x <- str_trim(x)
  
  # Extract first word
  genus <- word(x, 1)
  
  # TRUE for valid genera (capitalized Latin-like)
  is_valid <- str_detect(genus, "^[A-Z][a-z]+$")
  
  # Return genus if valid, otherwise NA
  ifelse(is_valid, genus, NA_character_)
}


# Helper: detect binomial species names
is_binomial <- function(x) {
  x <- str_trim(x)
  pattern <- "^[A-Za-z]+\\s+[A-Za-z\\-]+$"
  valid_species <- str_detect(x, pattern)
  flagged_terms <- str_detect(
    x,
    regex("\\b(sp|spp|species|cf|aff|indet|unknown)\\b", ignore_case = TRUE)
  )
  valid_species & !flagged_terms
}


# Helper: simple common-name stemmer for inverts / fuzzy things
stem_common <- function(name_raw) {
  # Make sure name_raw is character
  name_raw <- tolower(name_raw) %>% str_squish()
  
  # If NA or empty -> return NA
  is_bad <- is.na(name_raw) | name_raw == ""
  name_raw[is_bad] <- NA_character_
  
  # Strip qualifiers
  name_raw <- str_replace_all(
    name_raw,
    regex("(larvae|larva|nymph|case|cased|caseless|adult|juvenile)", ignore_case = TRUE),
    ""
  )
  
  name_raw <- str_replace_all(
    name_raw,
    regex("(common|unknown|unidentified|sp\\.|species)", ignore_case = TRUE),
    ""
  )
  
  name_raw <- str_squish(name_raw)
  
  # Extract first word core
  core <- word(name_raw, 1)
  
  # Standardise a few common groups
  core <- case_when(
    core %in% c("caddis", "caddisfly") ~ "caddisfly",
    core %in% c("fly", "true")         ~ "fly",
    core %in% c("spider")              ~ "spider",
    core %in% c("tick")                ~ "tick",
    core %in% c("aphid")               ~ "aphid",
    core %in% c("ant")                 ~ "ant",
    TRUE ~ core
  )
  
  # Replace empty strings with NA
  core[core == ""] <- NA_character_
  
  core
}


assoc <- assoc %>%
  mutate(
    # Clean scientific names
    sci_clean = case_when(
      !is.na(scientific_raw) & scientific_raw != "" &
        is_binomial(str_to_title(scientific_raw)) ~ str_to_title(scientific_raw),
      TRUE ~ str_to_title(scientific_raw)
    ),
    
    sci_genus   = extract_genus(sci_clean),
    sci_epithet = word(sci_clean, 2),
    epithet_lc  = tolower(coalesce(sci_epithet, "")),
    
    missing_epithet = (
      is.na(sci_epithet) |
        epithet_lc == "" |
        str_detect(epithet_lc, regex("^(sp|spp|cf|aff|indet|unknown)$",
                                     ignore_case = TRUE))
    ),
    
    # Common-name stem for cases with no reliable scientific genus
    common_stem = stem_common(common_name_raw),
    
    # Hybrid grouping key:
    #  - if scientific genus exists → use that
    #  - else fall back to common-name stem
    #  - else "unidentified"
    genus_group = case_when(
      !is.na(sci_genus)      ~ tolower(sci_genus),
      is.na(sci_genus) &
        !is.na(common_stem)  ~ common_stem,
      TRUE                   ~ "unidentified"
    ),
    
    # Final taxon unit:
    #  - true species keep their cleaned binomial
    #  - genus-only records collapse to "<group> sp."
    taxon_unit = case_when(
      !missing_epithet ~ str_to_lower(sci_clean),           # real species
      missing_epithet  ~ paste0(genus_group, " sp."),       # genus/common-group level
      TRUE             ~ NA_character_
    )
  )

# ---- Filter occurrences for analysis ----

assoc_use <- assoc %>%
  filter(
    occurrence_status == "Present",
    !is.na(region)
  )

assoc_use <- assoc_use %>%
  mutate(
    Pollutionscore    = suppressWarnings(as.numeric(Pollutionscore)),
    event_id          = str_replace_all(event_id, " ", "_"),
    common_name_title = str_to_title(common_name_raw)
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
      
      # Bats
      str_detect(oname_lc, "chiroptera") |
        str_detect(evid, "cr") |
        str_detect(bid, "audiomoth") ~ "Bats",
      
      # Birds
      str_detect(evid, bird_event_pat) ~ "Birds",
      
      # Aquatic inverts
      str_detect(evid, "fi") |
        str_detect(oname_lc, str_c(aquatic_orders, collapse = "|")) |
        str_detect(cname_lc,
                   "caddis|mayfly|stonefly|dragonfly|damselfly|boatman|diving beetle|dytisc|corix") ~
        "Aquatic invertebrates",
      
      # Mammals (non-bats)
      str_detect(oname_lc, str_c(setdiff(mammal_orders, "chiroptera"), collapse = "|")) ~ "Mammals",
      
      # Terrestrial inverts
      str_detect(oname_lc, str_c(terrestrial_orders, collapse = "|")) |
        str_detect(evid, "tin|tip|tit") ~ "Terrestrial invertebrates",
      
      TRUE ~ "Other"
    ),
    
    hillside = region  # alias to keep old naming
  )

# Only keep rows with usable taxon_unit
oc_use <- assoc_use %>%
  filter(
    !is.na(taxon_unit),
    taxon_unit != "",
    !str_detect(taxon_unit, regex("^no records|unidentified", ignore_case = TRUE))
  )

# ---- All-taxa plots (richness / record-based, not abundance) ----

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
  scale_fill_manual(values = c("North" = "#E76F51", "South" = "#0072B2")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title    = "Richness by sampled group and hillside",
    subtitle = "Number of distinct taxa per hillside (synthetic taxa for genus-only or group-level records)",
    x = "Richness (taxa)",
    y = "Group",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "top",
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
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
    "North" = "#E76F51",
    "South" = "#0072B2",
    "Tie"   = "grey80"
  )) +
  scale_x_continuous(
    labels = abs,
    breaks = forced_breaks,
    expand = expansion(mult = c(0.12, 0.12))
  ) +
  labs(
    x = "Difference in species richness ",
    y = "Group",
    fill = "Candidate site"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
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


# 6.3 Per-group "abundance" but as record counts, not individuals
# --- Grouped abundance function (records, not individualCount) ---

plot_group_abundance <- function(df, grp_label, top_n = Inf) {
  
  df_grp <- df %>% filter(group == grp_label)
  if (nrow(df_grp) == 0) return(NULL)
  
  # Count records per species per hillside
  dat0 <- df_grp %>%
    count(hillside, taxon_unit, name = "records")
  
  # Fill missing combinations
  dat <- dat0 %>%
    tidyr::complete(taxon_unit, hillside = c("North", "South"),
                    fill = list(records = 0)) %>%
    mutate(hillside = factor(hillside, levels = c("North", "South")))
  
  # Compute total records per species for ranking
  totals <- dat %>%
    group_by(taxon_unit) %>%
    summarise(total_records = sum(records), .groups = "drop")
  
  dat <- dat %>% left_join(totals, by = "taxon_unit")
  
  # Restrict to top N species if wanted
  if (is.finite(top_n)) {
    keep <- totals %>%
      arrange(desc(total_records)) %>%
      slice_head(n = top_n) %>%
      pull(taxon_unit)
    dat <- dat %>% filter(taxon_unit %in% keep)
  }
  
  dat <- dat %>%
    mutate(taxon_unit = reorder(taxon_unit, total_records))
  
  ggplot(dat, aes(x = taxon_unit, y = records, fill = hillside)) +
    geom_col(position = position_dodge(width = 0.7), 
             width = 0.65, colour = "white") +
    geom_text(aes(label = records),
              position = position_dodge(width = 0.7),
              hjust = -0.1, size = 3.2) +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    scale_fill_manual(values = c("North" = "#E76F51", "South" = "#0072B2")) +
    labs(
      title = paste0(grp_label, " — record-based abundance (species)"),
      subtitle = "Counts of sampling-event records per species",
      x = NULL,
      y = "Record count",
      fill = "Region"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line()
    )
}

major_groups <- c(
  "Birds",
  "Mammals",
  "Bats",
  "Plants",
  "Terrestrial invertebrates",
  "Aquatic invertebrates"
)

plots_grouped_abundance <- lapply(major_groups, function(g) {
  message("Plotting: ", g)
  plot_group_abundance(oc_use, g)
})

plots_grouped_abundance


#----overall abundance plot WITH OVERALL GROUP----

# 1. Compute abundance per group x hillside
abund_group_ns <- assoc_use %>%
  filter(group != "Other") %>%                      
  group_by(group, hillside) %>%
  summarise(
    total_individuals = sum(individualCount, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    hillside = factor(hillside, levels = c("North","South"))
  )

# 2. Create "Overall" totals for each hillside
overall_ns <- assoc_use %>%
  filter(group != "Other") %>%
  group_by(hillside) %>%
  summarise(
    total_individuals = sum(individualCount, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    group = "Overall",
    hillside = factor(hillside, levels = c("North","South"))
  )

# 3. Bind rows + set factor order (groups then Overall last)
abund_group_ns2 <- bind_rows(abund_group_ns, overall_ns) %>%
  mutate(
    group = factor(group,
                   levels = c(
                     sort(unique(abund_group_ns$group)),
                     "Overall"
                   ))
  )

# 4. Log breaks
log_breaks <- scales::log_breaks()(range(abund_group_ns2$total_individuals))

# 5. Plot
overall_abundance_plot <- ggplot(abund_group_ns2,
                                 aes(x = group, y = total_individuals, fill = hillside)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65,
           colour = "white") +
  geom_text(
    aes(label = total_individuals),
    position = position_dodge(width = 0.75),
    vjust = -0.4,
    size = 3.7,
    colour = "black"
  ) +
  scale_y_log10(
    breaks = log_breaks,
    labels = scales::comma_format(),
    expand = expansion(mult = c(0, 0.15))
  ) +
  scale_fill_manual(values = c(
    "North" = "#E76F51",
    "South" = "#0072B2"
  )) +
  labs(
    x = "Group",
    y = "Total individuals (log10 scale)",
    fill = "Candidate site"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

overall_abundance_plot


#----Invertebrate environemntal quality indicators----
####Freshwater Invertebrates####
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)

# clean data 
assoc_use <- assoc_use %>%
  mutate(
    Pollutionscore = suppressWarnings(as.numeric(Pollutionscore)),
    event_id = str_replace_all(event_id, " ", "_")
  )

fresh_inverts <- assoc_use %>%
  filter(
    group == "Aquatic invertebrates",
    !is.na(Pollutionscore),
    Pollutionscore > 0
  )

# aspt calculations 
aspt_sites <- fresh_inverts %>%
  group_by(event_id, hillside) %>%
  summarise(
    BMWP = sum(Pollutionscore, na.rm = TRUE),
    n_taxa = n(),
    ASPT = BMWP / n_taxa,
    .groups = "drop"
  )

aspt_totals <- fresh_inverts %>%
  group_by(hillside) %>%
  summarise(
    BMWP = sum(Pollutionscore, na.rm = TRUE),
    n_taxa = n(),
    ASPT = BMWP / n_taxa,
    .groups = "drop"
  ) %>%
  mutate(event_id = ifelse(hillside == "North", "N_total", "S_total"))

df <- bind_rows(aspt_sites, aspt_totals)

# site labels
site_labels <- c(
  "ND4FI_A" = "Downstream NE",
  "NU4FI_A" = "Upstream NE",
  "ND4FI_B" = "Downstream NW",
  "NU4FI_B" = "Upstream NW",
  "N_total" = "North Total",
  "SD4FI_A" = "Downstream SE",
  "SU4FI_A" = "Upstream SE",
  "SD4FI_B" = "South Bog",
  "SU4FI_B" = "Upstream of South Bog",
  "S_total" = "South Total"
)

small_gap <- 0.65
group_gap <- 0.95
big_gap <- 1.3

x_ND4FI_A <- 1
x_NU4FI_A <- x_ND4FI_A + small_gap
x_ND4FI_B <- x_NU4FI_A + group_gap
x_NU4FI_B <- x_ND4FI_B + small_gap
x_N_total <- x_NU4FI_B + group_gap

x_SD4FI_A <- x_N_total + big_gap
x_SU4FI_A <- x_SD4FI_A + small_gap
x_SD4FI_B <- x_SU4FI_A + group_gap
x_SU4FI_B <- x_SD4FI_B + small_gap
x_S_total <- x_SU4FI_B + group_gap

positions <- c(
  ND4FI_A = x_ND4FI_A,
  NU4FI_A = x_NU4FI_A,
  ND4FI_B = x_ND4FI_B,
  NU4FI_B = x_NU4FI_B,
  N_total = x_N_total,
  SD4FI_A = x_SD4FI_A,
  SU4FI_A = x_SU4FI_A,
  SD4FI_B = x_SD4FI_B,
  SU4FI_B = x_SU4FI_B,
  S_total = x_S_total
)

df <- df %>%
  mutate(
    x = positions[event_id],
    label = site_labels[event_id],
    is_total = event_id %in% c("N_total", "S_total")
  )

# heatmap data – now grouped by common_name_title (policy-friendly labels)
region_heat <- assoc_use %>%
  filter(group == "Aquatic invertebrates") %>%
  group_by(hillside, common_name_title) %>%
  summarise(
    abundance = sum(individualCount, na.rm = TRUE),
    bmwpscore = mean(Pollutionscore, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    hillside         = factor(hillside, levels = c("North", "South")),
    common_name_title = factor(common_name_title, levels = sort(unique(common_name_title)))
  )

common_scale <- scale_fill_gradient(
  low = "#d6e9f9",
  high = "#4C6FA3",
  name = "BMWP Score"
)

###ASPT bar plot ##
plot_aspt <- ggplot(df, aes(x = x, y = ASPT, fill = hillside)) +
  geom_col(aes(width = ifelse(is_total, 0.8, 0.55)), colour = "white") +
  geom_text(
    aes(label = round(ASPT, 2)),
    vjust = -0.3,
    fontface = ifelse(df$is_total, "bold", "plain"),
    size = 4
  ) +
  scale_fill_manual(values = c("North" = "#E76F51", "South" = "#0072B2")) +
  scale_x_continuous(
    breaks = df$x,
    labels = df$label,
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_y_continuous(
    limits = c(0, 8.25),
    expand = c(0, 0)
  ) +
  labs(
    x = "Site",
    y = "ASPT Score",
    fill = "Candidate site"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.line.x = element_line(colour = "black", linewidth = 1),
    axis.line.y = element_line(colour = "black", linewidth = 1),
    axis.ticks = element_line(colour = "black"),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 18),
    legend.position = "right"
  )

# combined heatmap – y-axis is now common (English) names
heatmap_combined <- ggplot(region_heat,
                           aes(x = hillside,
                               y = common_name_title,
                               fill = bmwpscore)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = abundance),
            colour = "black",
            fontface = "bold",
            size = 4) +
  common_scale +
  labs(
    x = "Candidate site",
    y = "Taxon (common name)"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 18),
    legend.position = "right"
  )


# print plots
plot_aspt
heatmap_combined





####Moths, terrestrial invertebrates and Overall environmental comparison####
library(tidyverse)

# ---- Derive Terrestrial Invertebrate Traits from assoc_use ----

# expects assoc_use with:
# - group (includes "Terrestrial invertebrates")
# - hillside (North / South)
# - order_name, common_name, scientific_name
# - individualCount (abundance)

terrestrial_traits <- assoc_use %>%
  filter(
    group == "Terrestrial invertebrates",
    !is.na(individualCount),
    individualCount > 0,
    !str_detect(order_name, regex("lepidoptera", ignore_case = TRUE))
  ) %>%
  mutate(
    order_lc = tolower(order_name),
    cname_lc = tolower(coalesce(common_name, "")),
    sname_lc = tolower(coalesce(scientific_name, "")),
    
    woodland_score = case_when(
      str_detect(cname_lc, "woodlouse|woodlouse|slug|snail|millipede|centipede") ~ 2,
      str_detect(order_lc, "coleoptera|araneae|opiliones") ~ 2,
      str_detect(order_lc, "hymenoptera|hemiptera|diptera") ~ 1,
      TRUE ~ 0
    ),
    
    moisture_score = case_when(
      str_detect(cname_lc, "woodlouse|slug|snail|millipede|centipede") ~ 2,
      str_detect(order_lc, "coleoptera|araneae") ~ 1,
      TRUE ~ 0
    ),
    
    trophic_score = case_when(
      str_detect(order_lc, "araneae|opiliones") ~ 2,                                 # spiders, harvestmen
      str_detect(cname_lc, "ground beetle|carabid") ~ 2,                             # ground beetles
      str_detect(cname_lc, "woodlouse|millipede|centipede|earthworm|slug|snail") ~ 1,# detritivores
      str_detect(order_lc, "coleoptera") ~ 1.5,                                      # many beetles pred/scav
      str_detect(order_lc, "hymenoptera") ~ 1,                                       # ants, many omnivores
      TRUE ~ 0.5                                                                     # default herbivore/low-signal
    ),
    
    saproxylic_score = case_when(
      str_detect(cname_lc, "longhorn|click beetle|dor beetle|bark beetle") ~ 2,
      str_detect(order_lc, "coleoptera") ~ 1,
      TRUE ~ 0
    ),
    
    log_abund = log1p(individualCount),
    composite_score = woodland_score + moisture_score + trophic_score + saproxylic_score
  )

habitat_quality <- terrestrial_traits %>%
  group_by(hillside) %>%
  summarise(
    woodland_raw   = sum(woodland_score   * log_abund, na.rm = TRUE),
    moisture_raw   = sum(moisture_score   * log_abund, na.rm = TRUE),
    trophic_raw    = sum(trophic_score    * log_abund, na.rm = TRUE),
    saproxylic_raw = sum(saproxylic_score * log_abund, na.rm = TRUE),
    habitat_raw    = sum(composite_score  * log_abund, na.rm = TRUE),
    .groups = "drop"
  )

# ---- Terrestrial Invertebrate Metric Breakdown ----

habitat_scaled <- habitat_quality %>%
  mutate(
    woodland      = woodland_raw   / max(woodland_raw,   na.rm = TRUE),
    moisture      = moisture_raw   / max(moisture_raw,   na.rm = TRUE),
    trophic       = trophic_raw    / max(trophic_raw,    na.rm = TRUE),
    saproxylic    = saproxylic_raw / max(saproxylic_raw, na.rm = TRUE),
    Average       = habitat_raw    / max(habitat_raw,    na.rm = TRUE)   # renamed
  ) %>%
  select(hillside, woodland, moisture, trophic, saproxylic, Average)      # reordered

hab_plot <- habitat_scaled %>%
  pivot_longer(
    cols = -hillside,
    names_to = "metric",
    values_to = "score"
  ) %>%
  mutate(
    metric = factor(metric, 
                    levels = c("woodland", "moisture", "trophic", 
                               "saproxylic", "Average"))  # Average moved to end
  )

# ---- Terrestrial Invertebrate Plot with Solid Axes ----

p_terrestrial <- ggplot(hab_plot,
                        aes(x = metric, y = score, fill = hillside)) +
  geom_col(
    width = 0.6,
    position = position_dodge(width = 0.75)
  ) +
  geom_text(
    aes(label = round(score, 2)),
    position = position_dodge(width = 0.75),
    vjust = -0.2,
    size = 4
  ) +
  scale_fill_manual(values = c("North" = "#E76F51", "South" = "#0072B2")) +
  scale_y_continuous(
    limits = c(0, 1.05),
    expand = c(0, 0)
  ) +
  labs(
    x = "Metrics for environmental quality",
    y = "Scaled score (0–1)",
    fill = "Candidate site"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid = element_blank(),
    axis.line.x = element_line(colour = "black", linewidth = 1),
    axis.line.y = element_line(colour = "black", linewidth = 1)
  )

p_terrestrial

# ---- Environmental Comparison Plot with Solid Axes + Non-Overlapping Bars ----
# ---- Build Moth Air-Quality (Nitrogen Pollution) Score ----
# ---- Moth nitrogen air-quality score ----
# ---- UK Nitrogen Sensitivity Classification ----
# ---- Species name standardisation table ----
name_corrections <- tribble(
  ~raw,                   ~clean,
  "autuml rustic",        "autumn rustic",
  "autumn rustic",        "autumn rustic",
  "black rustic",         "black rustic",
  "square spot rustic",   "square spot rustic",
  "drinker moth",         "drinker moth",
  "garden dart",          "garden dart",
  "chestnut moth",        "chestnut moth",
  "dark arches",          "dark arches",
  
  # tolerant species
  "common marbled carpet", "common marbled carpet",
  "green carpet moth",     "green carpet moth",
  "marbled carpet",        "marbled carpet",
  "footman",               "footman",
  "footman moth",          "footman",
  "pug",                   "pug",
  "pug moth",              "pug",
  
  # neutral / other
  "caterpillar",           "caterpillar",
  "larva",                 "caterpillar",
  "crane fly",             "crane fly",
  "convolvulus hawk moth", "convolvulus hawk moth"
)

sensitive_species <- c(
  "drinker moth",
  "black rustic",
  "square spot rustic",
  "garden dart",
  "chestnut moth",
  "autumn rustic",
  "dark arches",
  "shuttle-shaped dart"
)

tolerant_species <- c(
  "common marbled carpet",
  "green carpet moth",
  "marbled carpet",
  "footman",
  "pug"
)

# scoring:
# 3 = sensitive (good)
# 2 = neutral
# 1 = tolerant (bad)

# ---- Robust nitrogen scoring with automatic name correction ----

moth_data <- assoc_use %>%
  filter(
    str_detect(order_name, regex("lepidoptera", ignore_case = TRUE)),
    !is.na(individualCount),
    individualCount > 0
  ) %>%
  mutate(
    cname_lc = tolower(common_name)
  ) %>%
  # apply name corrections
  left_join(name_corrections,
            by = c("cname_lc" = "raw")) %>%
  mutate(
    clean_name = coalesce(clean, cname_lc),
    
    nitrogen_score = case_when(
      clean_name %in% sensitive_species ~ 3,   # sensitive = high score
      clean_name %in% tolerant_species ~ 1,    # tolerant = low score
      TRUE ~ 2                                 # neutral
    ),
    
    air_weighted = nitrogen_score * log1p(individualCount)
  )




air_quality <- moth_data %>%
  group_by(hillside) %>%
  summarise(
    score_raw = sum(air_weighted, na.rm = TRUE),
    .groups = "drop"
  )

combined_scores_raw <- air_quality %>%
  rename(air_raw = score_raw) %>%
  left_join(
    habitat_scaled %>% select(hillside, Average),  # from your updated terrestrial plot
    by = "hillside"
  ) %>%
  rename(habitat_raw = Average) %>%
  left_join(
    aspt_totals %>%
      select(hillside, ASPT) %>%
      rename(freshwater_raw = ASPT),
    by = "hillside"
  )

# ---- Wide format for scaling ----

w <- combined_scores_raw %>%
  pivot_wider(
    names_from = hillside,
    values_from = c(air_raw, habitat_raw, freshwater_raw)
  )

# ---- Scaling helpers ----

scale_to_one <- function(xN, xS) {
  maxv <- max(xN, xS, na.rm = TRUE)
  list(N = xN / maxv, S = xS / maxv)
}

safe_extract <- function(df, colname) {
  df %>% pull({{colname}}) %>% unique() %>% max(na.rm = TRUE)
}

# ---- Apply scaling to all 3 metrics ----

air_scaled <- scale_to_one(
  safe_extract(w, air_raw_North),
  safe_extract(w, air_raw_South)
)

hab_scaled <- scale_to_one(
  safe_extract(w, habitat_raw_North),
  safe_extract(w, habitat_raw_South)
)

fresh_scaled <- scale_to_one(
  safe_extract(w, freshwater_raw_North),
  safe_extract(w, freshwater_raw_South)
)

# ---- Combined scaled values ----

scaled_scores <- tibble(
  test = c(
    "Air Pollution (Moths)",
    "Freshwater Quality (ASPT)",
    "Habitat Quality (Terrestrial Inverts)"
  ),
  North = c(air_scaled$N, fresh_scaled$N, hab_scaled$N),
  South = c(air_scaled$S, fresh_scaled$S, hab_scaled$S)
)

# ---- Overall score ----

overall_scores <- tibble(
  test = "Overall Environmental Score",
  North = mean(scaled_scores$North, na.rm = TRUE),
  South = mean(scaled_scores$South, na.rm = TRUE)
)

# ---- Long format for plotting ----

plot_scores <- scaled_scores %>%
  bind_rows(overall_scores) %>%
  pivot_longer(
    cols = c(North, South),
    names_to = "hillside",
    values_to = "scaled_score"
  )

# ---- Plot: All environmental metrics ----

p_all_env <- ggplot(
  plot_scores,
  aes(x = test, y = scaled_score, fill = hillside)
) +
  geom_col(
    width = 0.6,
    position = position_dodge(width = 0.75)
  ) +
  geom_text(
    aes(label = round(scaled_score, 2)),
    position = position_dodge(width = 0.75),
    vjust = -0.2,
    size = 4
  ) +
  scale_fill_manual(values = c("North" = "#E76F51", "South" = "#0072B2")) +
  scale_y_continuous(
    limits = c(0, 1.05),
    expand = c(0, 0)
  ) +
  labs(
    x = "Environmental quality metric",
    y = "Scaled score (0–1)",
    fill = "Candidate site"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(colour = "black", linewidth = 1),
    axis.line.y = element_line(colour = "black", linewidth = 1)
  )

p_all_env




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
  "kestrel", "falco tinnunculus", "Red", "Falconiformes", "Falconidae",
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
# ---- Normalise names & aggregate per event for matching ----

synonyms <- tribble(
  ~common_name_obs,        ~common_name_std,
  "common kestrel",        "kestrel",
  "peregrine",             "peregrine falcon",
  "nathusius pipistrelle", "nathusius' pipistrelle"
)



assoc_imp <- assoc %>%
  mutate(
    common_name     = tolower(str_squish(common_name)),
    scientific_name = tolower(str_squish(scientific_name))
  ) %>%
  left_join(synonyms, by = c("common_name" = "common_name_obs")) %>%
  mutate(common_name = coalesce(common_name_std, common_name)) %>%
  select(-common_name_std)

# Aggregate to one row per event × species (common name) first
obs_for_match <- assoc_imp %>%
  filter(!is.na(common_name)) %>%   # we match by common name for a non-specialist output
  group_by(event_id, region, common_name, scientific_name) %>%
  summarise(individualCount = sum(individualCount, na.rm = TRUE), .groups = "drop")

# ---- Exact join only (no fuzzy) to avoid over-counting ----

important_species_all <- obs_for_match %>%
  inner_join(lookup_all,
             by = c("common_name" = "common_name")) %>%
  transmute(
    event_id,
    region,
    group,                     # Bird / Mammal
    common_name,
    scientific_name = coalesce(scientific_name.y, scientific_name.x),
    status,
    order,
    family,
    individualCount
  )


# ---- Setup: important-only dataset for plots ----

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
}

dat <- important_species_all %>%
  filter(!str_detect(tolower(event_id), "\\boob\\b")) %>%
  mutate(
    group = factor(group, levels = c("Bird","Mammal")),
    region = factor(region),
    species_label = str_to_title(common_name)  # human-friendly
  )

pal_region <- c("North" = "#E76F51", "South" = "#0072B2") %>%
  { .[intersect(names(.), levels(dat$region))] }
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

# ---- Plot: Richness (important only, by common name) ----

richness <- dat %>%
  distinct(region, group, species_label) %>%
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
    title = "Important species richness by region and group",
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

# ---- Plot: Abundance (important only – TRUE abundance) ----

top_n_per_group <- 10

abund_species <- dat %>%
  group_by(group, region, species_label) %>%
  summarise(n_records = sum(individualCount, na.rm = TRUE), .groups = "drop") %>%
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
    title = "Abundance of important species by region",
    x = "Total individuals (summed individualCount)",
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

# ---- Plot: Heatmap (important only – TRUE abundance) ----

species_status <- dat %>%
  mutate(category = case_when(
    group == "Bird" & status %in% c("Red", "Amber") ~ status,
    group == "Bird"                                 ~ "Other",
    group == "Mammal" & str_detect(status, regex("EPS", ignore_case = TRUE)) ~ "EPS",
    group == "Mammal" & str_detect(status, regex("Protected|WCA|PBA", ignore_case = TRUE)) ~ "Protected",
    group == "Mammal" & str_detect(status, regex("SBL", ignore_case = TRUE)) ~ "SBL Priority",
    group == "Mammal" ~ "Other",
    TRUE ~ "Other"
  )) %>%
  distinct(group, species_label, category)

heat <- dat %>%
  group_by(group, species_label, region) %>%
  summarise(n_records = sum(individualCount, na.rm = TRUE), .groups = "drop") %>%
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
heat

# ---- Clean common names for heatmap ----
heat <- heat %>%
  mutate(
    common_name_clean = species_label %>%
      str_replace_all("_", " ") %>%   # remove underscores if present
      str_squish() %>%                # clean whitespace
      str_to_title()                  # make readable
  )


# ---- Colour palette ----
pal_comp <- c(
  "Amber"        = "#FFB000",
  "Red"          = "#C1121F",
  "EPS"          = "#1B9E77",
  "Protected"    = "#7570B3",
  "SBL Priority" = "#E7298A",
  "Other"        = "grey80"
)

# ---- Heatmap using common names ----
fig_heatmapimportant <- ggplot(heat, aes(x = region, y = common_name_clean)) +
  geom_tile(aes(fill = category),
            color = "grey",
            linewidth = 0.35) +
  geom_text(aes(label = n_records),
            color = "black",
            size = 3.6,
            fontface = "bold") +
  facet_wrap(~ group, scales = "free_y", ncol = 1, strip.position = "top") +
  scale_fill_manual(values = pal_comp, name = "Conservation Status") +
  labs(
    title = "Conservation Concern Species Heatmap",
    x = "Candidate Site",
    y = "Species (Common Names)"
  ) +
  theme_ecia() +
  theme(
    # REMOVE ALL GREY PANEL LINES
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.line.x = element_line(colour = "black", linewidth = 0.4),
    axis.line.y = element_line(colour = "black", linewidth = 0.4),
    axis.ticks = element_line(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    panel.spacing.y = unit(1.2, "lines"),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, margin = margin(b = 10)),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 9),
    plot.margin = margin(10, 15, 10, 10)
  )

fig_heatmapimportant





# ---- Plot: Bird BoCC composition (richness, not abundance) ----

birds_comp_counts <- dat %>%
  filter(group == "Bird", status %in% c("Red","Amber")) %>%
  mutate(status = factor(status, levels = c("Amber","Red"))) %>%
  distinct(region, status, species_label) %>%   # richness by common name
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
    aes(region, total, label = paste0("Total species: ", total)),
    vjust = -0.5, size = 3.2, fontface = "bold", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = pal_bocc) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(
    title = "Bird Species of Conservation Concern",
    x = "Candidate Site", y = "Number of Bird species", fill = "BoCC status"
  ) +
  theme_ecia()

fig_bocc_comp_counts

# ---- Plot: Mammal protection composition (richness, not abundance) ----

mammals_comp_counts <- dat %>%
  filter(group == "Mammal") %>%
  mutate(protection = case_when(
    str_detect(status, regex("EPS", ignore_case = TRUE)) ~ "EPS",
    str_detect(status, regex("Protected|WCA|PBA", ignore_case = TRUE)) ~ "Protected",
    str_detect(status, regex("SBL", ignore_case = TRUE)) ~ "SBL Priority",
    TRUE ~ status
  )) %>%
  mutate(protection = factor(protection, levels = c("SBL Priority","Protected","EPS"))) %>%
  distinct(region, protection, species_label) %>%
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
    aes(region, total, label = paste0("Total species: ", total)),
    vjust = -0.5, size = 3.2, fontface = "bold", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = pal_protect) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(
    title = "Important mammal species by protection category and region",
    x = "Region", y = "Number of important mammal species", fill = "Category"
  ) +
  theme_ecia()

fig_mammal_comp_counts

# ---- Combined composition (richness) ----

birds_comp_counts2 <- dat %>%
  filter(group == "Bird", status %in% c("Red","Amber")) %>%
  mutate(category = factor(status, levels = c("Amber","Red"))) %>%
  distinct(group, region, category, species_label) %>%
  count(group, region, category, name = "n")

mammals_comp_counts2 <- dat %>%
  filter(group == "Mammal") %>%
  mutate(category = case_when(
    str_detect(status, regex("EPS", ignore_case = TRUE)) ~ "EPS",
    str_detect(status, regex("Protected|WCA|PBA", ignore_case = TRUE)) ~ "Protected",
    str_detect(status, regex("SBL", ignore_case = TRUE)) ~ "SBL Priority",
    TRUE ~ status
  )) %>%
  mutate(category = factor(category, levels = c("SBL Priority","Protected","EPS"))) %>%
  distinct(group, region, category, species_label) %>%
  count(group, region, category, name = "n")

comp_all <- bind_rows(birds_comp_counts2, mammals_comp_counts2) %>%
  mutate(
    region = factor(region),
    group = factor(group, levels = c("Bird","Mammal")),
    category = factor(category, levels = c("SBL Priority","Protected","EPS","Amber","Red"))
  )

comp_totals <- comp_all %>%
  group_by(group, region) %>%
  summarise(total = sum(n), .groups = "drop")

pal_comp2 <- c(
  "Red"          = "#C1121F",   # danger / high concern
  "Amber"        = "#FFB000",   # caution / medium concern
  "EPS"          = "#1B9E77",   # teal
  "Protected"    = "#7570B3"   # purple
  
)

fig_comp_combinedimportantbyregionandstatus <- ggplot(comp_all, aes(region, n, fill = category)) +
  geom_col(
    width = 0.7,
    color = "grey90",
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
    aes(region, total, label = paste0("Total species: ", total)),
    vjust = -0.6,
    size = 3.4,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ group, ncol = 2, scales = "fixed") +
  scale_fill_manual(values = pal_comp2, drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Richness of concern species between candidate sites",
    x = "Candidate site",
    y = "Number of species",
    fill = "Conservation status"
  ) +
  theme_ecia() +
  theme(
    # REMOVE ALL GREY GRIDLINES
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.line = element_line(colour = "black", linewidth = 0.4),
    axis.ticks = element_line(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 15, 10, 10)
  )

fig_comp_combinedimportantbyregionandstatus






# ---- Potential plots for report ----
p_rich_all
p_diverge_final
overall_abundance_plot
plot_aspt
heatmap_combined
p_terrestrial
p_all_env
fig_abundanceimportant #probably wont use
fig_heatmapimportant 
fig_comp_combinedimportantbyregionandstatus

#----finalised plots for report----
setwd("C:/Users/stuar/Documents/Masters/Arran/R script/Arran field course/Final plots")

#figure 1
p_diverge_final
overall_abundance_plot

#Saving
library(ggplot2)
ggsave("p_diverge_final_A4_half.png",
       plot = p_diverge_final,
       width = 5.8, height = 8.27, units = "in", dpi = 300)

ggsave("overall_abundance_plot_A4_half.png",
       plot = overall_abundance_plot,
       width = 5.8, height = 8.27, units = "in", dpi = 300)



#Figure 2
plot_aspt
heatmap_combined
p_terrestrial
p_all_env

#Figure 3
fig_heatmapimportant 
fig_comp_combinedimportantbyregionandstatus




