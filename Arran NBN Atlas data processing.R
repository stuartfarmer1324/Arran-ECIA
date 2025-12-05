# ----------------------------------------------------------
# Packages
# ----------------------------------------------------------
library(tidyverse)
library(vroom)
library(janitor)
library(scales)

# ----------------------------------------------------------
# Paths
# ----------------------------------------------------------
out_dir <- "C:/Users/stuar/Documents/Masters/Arran/Data/NBN Atlas"

north_path      <- file.path(out_dir, "North species list.txt")
south_path      <- file.path(out_dir, "updated south candidate site.txt")
north_500_path  <- file.path(out_dir, "Terrestrial buffer North.txt")
south_500_path  <- file.path(out_dir, "TerS.txt")
north_2000_path <- file.path(out_dir, "ABUFN.txt")
south_2000_path <- file.path(out_dir, "ABUFS.txt")

# ----------------------------------------------------------
# Global colour palettes
# ----------------------------------------------------------
pal_hillside <- c(
  "North" = "#E76F51",
  "South" = "#2A9D8F"
)

pal_all <- c(
  "EPS"       = "#1B9E77",
  "Protected" = "#7570B3",
  "Priority"  = "#E7298A",
  "Red"       = "#C1121F",
  "Amber"     = "#FFB000",
  "Green"     = "#4C6FA3",
  "Other"     = "grey80"
)

# ----------------------------------------------------------
# Safe Reader
# ----------------------------------------------------------
robust_read <- function(path) {
  vroom(
    path,
    delim = ",",
    col_types = cols(.default = "c"),
    trim_ws = TRUE,
    show_col_types = FALSE,
    progress = FALSE
  )
}

# ----------------------------------------------------------
# Normalise scientific names
# ----------------------------------------------------------
clean_sci <- function(x) {
  x %>%
    tolower() %>%
    str_squish() %>%
    str_replace("\\s+agg\\.?$", "") %>%
    str_replace("\\s+subsp\\..*", "") %>%
    str_replace("\\s+var\\..*", "") %>%
    str_replace("\\s*\\(.*?\\)", "") %>%
    str_replace("=.+", "") %>%
    str_replace_all("/", " ") %>%
    str_replace_all("[^a-z\\s]", " ") %>%
    str_squish() %>%
    word(1, 2)
}

# ----------------------------------------------------------
# Read NBN files (generalised)
# ----------------------------------------------------------
read_nbn <- function(path, site_id) {
  raw <- robust_read(path)
  df  <- clean_names(raw)
  
  if (!"lsid" %in% names(df)) df$lsid <- NA
  
  lsid_parts <- str_split_fixed(df$lsid %||% "", "\\|", 5)
  colnames(lsid_parts) <- c("lsid_sci","lsid_key","lsid_vern","lsid_k","lsid_f")
  df <- bind_cols(df, as_tibble(lsid_parts))
  
  df <- df %>%
    mutate(
      site_id         = site_id,
      scientific_name = coalesce(species_name, lsid_sci) |> str_squish(),
      vernacular_name = coalesce(vernacular_name, lsid_vern),
      kingdom         = coalesce(kingdom, lsid_k),
      family          = coalesce(family,  lsid_f)
    ) %>%
    filter(!is.na(scientific_name), scientific_name != "") %>%
    distinct()
  
  df
}

# ----------------------------------------------------------
# Assign major taxonomic group
# ----------------------------------------------------------
assign_group <- function(df) {
  
  df %>%
    mutate(
      sci = tolower(scientific_name),
      fam = tolower(family),
      gen = tolower(genus),
      kng = tolower(kingdom),
      
      group = case_when(
        # Birds
        str_detect(sci, "anas|falco|turdus|corvus|larus|rissa|sterna|haematopus|numenius|podiceps|ardea|egretta|hirundo|prunella|muscicapa") ~ "Birds",
        
        # Bats
        str_detect(sci, "pipistrell|myotis|plecotus|nyctalus|eptesicus") ~ "Bats",
        
        # Mammals (non-bat)
        str_detect(sci, "sciurus|vulpes|mustela|meles|lutra|castor|lepus|erinaceus|capreolus") ~ "Mammals",
        
        # Amphibians
        str_detect(sci, "bufo|rana|lissotriton") ~ "Amphibians",
        
        # Reptiles
        str_detect(sci, "vipera|natrix|lacerta") ~ "Reptiles",
        
        # Fish
        fam %in% c("salmonidae","gobiidae","anguillidae") ~ "Fish",
        
        # Aquatic inverts
        str_detect(sci, "dytisc|corix|gammar|baet") ~ "Aquatic invertebrates",
        
        # Plants
        kng == "plantae" ~ "Plants",
        
        # Fungi
        kng == "fungi" ~ "Fungi",
        
        TRUE ~ "Terrestrial invertebrates"
      )
    )
}

# ----------------------------------------------------------
# Read candidate site datasets
# ----------------------------------------------------------
north_site <- read_nbn(north_path, "North")
south_site <- read_nbn(south_path, "South")

sites <- bind_rows(north_site, south_site) %>%
  mutate(hillside = factor(site_id, levels = c("North","South")))

site_richness <- sites %>%
  distinct(hillside, scientific_name) %>%
  count(hillside, name = "richness")

# ----------------------------------------------------------
# Read buffer datasets (AOI)
# ----------------------------------------------------------
north500_raw  <- read_nbn(north_500_path,  "North")
south500_raw  <- read_nbn(south_500_path,  "South")
north2000_raw <- read_nbn(north_2000_path, "North")
south2000_raw <- read_nbn(south_2000_path, "South")

# ----------------------------------------------------------
# Apply zone rules
# ----------------------------------------------------------

# 500 m: remove Birds, Bats, Plants
north500 <- assign_group(north500_raw) %>%
  filter(!group %in% c("Birds","Bats","Plants")) %>%
  mutate(zone="500m", hillside="North")

south500 <- assign_group(south500_raw) %>%
  filter(!group %in% c("Birds","Bats","Plants")) %>%
  mutate(zone="500m", hillside="South")

# 2000 m: only Birds + Bats
north2000 <- assign_group(north2000_raw) %>%
  filter(group %in% c("Birds","Bats")) %>%
  mutate(zone="2000m", hillside="North")

south2000 <- assign_group(south2000_raw) %>%
  filter(group %in% c("Birds","Bats")) %>%
  mutate(zone="2000m", hillside="South")

# ----------------------------------------------------------
# Combine AOI buffers
# ----------------------------------------------------------
buffers <- bind_rows(north500, south500, north2000, south2000)

# ----------------------------------------------------------
# AOI richness summaries
# ----------------------------------------------------------
aoi_rich_overall <- buffers %>%
  distinct(hillside, scientific_name) %>%
  count(hillside, name="richness")

aoi_group_rich <- buffers %>%
  distinct(hillside, group, scientific_name) %>%
  count(hillside, group, name="richness") %>%
  filter(group != "Reptiles")

# ----------------------------------------------------------
# Candidate site richness by group
# ----------------------------------------------------------
sites_group <- assign_group(sites) %>%
  distinct(hillside, group, scientific_name) %>%
  count(hillside, group, name="richness") %>%
  filter(group != "Reptiles")

# ----------------------------------------------------------
# AOI richness plot (overall)
# ----------------------------------------------------------
max_aoi_rich <- max(aoi_rich_overall$richness, na.rm=TRUE)

aoirichness <- ggplot(aoi_rich_overall,
                      aes(hillside, richness, fill=hillside)) +
  geom_col(colour="white") +
  geom_text(aes(label=richness), vjust=-0.3, size=4) +
  scale_fill_manual(values=pal_hillside) +
  scale_y_continuous(limits=c(0,max_aoi_rich*1.10), expand=c(0,0)) +
  labs(x="Area of interest", y="Species richness") +
  theme_minimal(base_size=13) +
  theme(panel.grid=element_blank(),
        axis.line=element_line(colour="black"),
        axis.ticks=element_line(colour="black"))
aoirichness



# ---- Packages ----
library(tidyverse)
library(vroom)
library(janitor)
library(scales)

# ---- Paths ----
out_dir <- "C:/Users/stuar/Documents/Masters/Arran/Data/NBN Atlas"

north_path      <- file.path(out_dir, "North species list.txt")
south_path      <- file.path(out_dir, "updated south candidate site.txt")
north_500_path  <- file.path(out_dir, "Terrestrial buffer North.txt")
south_500_path  <- file.path(out_dir, "TerS.txt")
north_2000_path <- file.path(out_dir, "ABUFN.txt")
south_2000_path <- file.path(out_dir, "ABUFS.txt")

# ---- Global colour palettes ----
pal_hillside <- c(
  "North" = "#E76F51",
  "South" = "#2A9D8F"
)

pal_all <- c(
  "EPS"       = "#2D6A4F",
  "Protected" = "#1B4332",
  "Priority"  = "#40916C",
  "Red"       = "#C1121F",
  "Amber"     = "#FFB000",
  "Green"     = "#74C69D",
  "Other"     = "grey80"
)

# ---- Safe Reader ----
robust_read <- function(path) {
  vroom(
    path,
    delim = ",",
    col_types = cols(.default = "c"),
    trim_ws = TRUE,
    show_col_types = FALSE,
    progress = FALSE
  )
}

# ---- Read & Clean NBN ----
read_nbn <- function(path, site_id) {
  
  raw <- robust_read(path)
  df  <- clean_names(raw)
  
  if (!"lsid" %in% names(df)) df$lsid <- NA
  
  lsid_parts <- str_split_fixed(df$lsid %||% "", "\\|", 5)
  colnames(lsid_parts) <- c("lsid_sci", "lsid_key", "lsid_vern", "lsid_k", "lsid_f")
  df <- bind_cols(df, as_tibble(lsid_parts))
  
  if (!"species_name"     %in% names(df)) df$species_name     <- NA
  if (!"vernacular_name"  %in% names(df)) df$vernacular_name  <- NA
  if (!"kingdom"          %in% names(df)) df$kingdom          <- NA
  if (!"family"           %in% names(df)) df$family           <- NA
  if (!"genus"            %in% names(df)) df$genus            <- NA
  
  df <- df %>%
    mutate(
      site_id         = site_id,
      scientific_name = coalesce(species_name, lsid_sci) |> str_squish(),
      vernacular_name = coalesce(vernacular_name, lsid_vern),
      kingdom         = coalesce(kingdom, lsid_k),
      family          = coalesce(family,  lsid_f)
    ) %>%
    filter(!is.na(scientific_name), scientific_name != "") %>%
    distinct()
  
  df
}

# ---- Assign Taxonomic Group ----
assign_group <- function(df) {
  
  df %>%
    mutate(
      sci = tolower(scientific_name),
      fam = tolower(family),
      gen = tolower(genus),
      kng = tolower(kingdom),
      
      group = case_when(
        # Birds
        str_detect(sci, "anas|aegithalos|accipiter|acrocephalus|actitis|alauda|anthus|emberiza|falco|muscicapa|prunella|strix|turdus|corvus|motacilla|larus|rissa|sterna|haematopus|numenius|podiceps|ardea|egretta|hirundo|delichon|riparia") ~ "Birds",
        
        # Bats
        str_detect(sci, "pipistrell|myotis|plecotus|nyctalus|eptesicus") ~ "Bats",
        
        # Mammals (non-bat)
        str_detect(sci, "sciurus|vulpes|mustela|meles|oryctolagus|capreolus|sus|martes|arvicola|lutra|castor|felis|erinaceus|lepus|cervus|phoca|halichoerus|phocoena|apodemus|myodes|sorex") ~ "Mammals",
        
        # Amphibians
        str_detect(sci, "bufo|rana|lissotriton|ichthyosaura|triturus") ~ "Amphibians",
        fam == "bufonidae" ~ "Amphibians",
        
        # Reptiles (kept in data, dropped at plotting)
        str_detect(sci, "natrix|vipera|zootoca|lacerta") ~ "Reptiles",
        
        # Fish
        fam %in% c("salmonidae","gobiidae","anguillidae","gasterosteidae","pleuronectidae") ~ "Fish",
        str_detect(sci, "anguilla|gasterosteus|salmo|platichthys|pomatoschistus") ~ "Fish",
        
        # Aquatic invertebrates
        fam %in% c("hydrachnidae","hydrachinidae","dytiscidae","corixidae","gerridae","notonectidae") ~ "Aquatic invertebrates",
        str_detect(sci, "hydra|dytisc|corix|gammar|baet|heptagen|oligochaet") ~ "Aquatic invertebrates",
        
        # Plants
        kng == "plantae" ~ "Plants",
        
        # Fungi
        kng == "fungi" ~ "Fungi",
        
        # Everything else
        TRUE ~ "Terrestrial invertebrates"
      )
    )
}

# ---- Read Candidate Site Data ----
north_site <- read_nbn(north_path, "North")
south_site <- read_nbn(south_path, "South")

sites <- bind_rows(north_site, south_site) %>%
  mutate(hillside = factor(site_id, levels = c("North","South")))

site_richness <- sites %>%
  distinct(hillside, scientific_name) %>%
  count(hillside, name = "richness")

# ---- Read Raw Buffer Data (AOI) ----
north500_raw  <- read_nbn(north_500_path,  "North")
south500_raw  <- read_nbn(south_500_path,  "South")
north2000_raw <- read_nbn(north_2000_path, "North")
south2000_raw <- read_nbn(south_2000_path, "South")

# ---- Filter Before Combining (AOI zone rules) ----
# 500 m: exclude Birds, Bats, Plants
north500 <- assign_group(north500_raw) %>%
  filter(!group %in% c("Birds", "Bats", "Plants")) %>%
  mutate(zone = "500m", hillside = "North")

south500 <- assign_group(south500_raw) %>%
  filter(!group %in% c("Birds", "Bats", "Plants")) %>%
  mutate(zone = "500m", hillside = "South")

# 2000 m: include Birds & Bats
north2000 <- assign_group(north2000_raw) %>%
  filter(group %in% c("Birds", "Bats")) %>%
  mutate(zone = "2000m", hillside = "North")

south2000 <- assign_group(south2000_raw) %>%
  filter(group %in% c("Birds", "Bats")) %>%
  mutate(zone = "2000m", hillside = "South")

# ---- Combine Cleaned Buffers (AOI dataset) ----
buffers <- bind_rows(north500, south500, north2000, south2000)

# ---- AOI Richness Summaries ----
aoi_rich_overall <- buffers %>%
  distinct(hillside, scientific_name) %>%
  count(hillside, name = "richness")

aoi_group_rich <- buffers %>%
  distinct(hillside, group, scientific_name) %>%
  count(hillside, group, name = "richness") %>%
  filter(group != "Reptiles")

# ---- Candidate Site Richness by Group ----
sites_group <- assign_group(sites) %>%
  distinct(hillside, group, scientific_name) %>%
  count(hillside, group, name = "richness") %>%
  filter(group != "Reptiles")

# ---- Candidate Site Overall Richness Plot ----
max_site_rich <- max(site_richness$richness, na.rm = TRUE)

candidaterichness <- ggplot(site_richness,
                            aes(hillside, richness, fill = hillside)) +
  geom_col(show.legend = FALSE, colour = "white") +
  geom_text(aes(label = richness),
            vjust = -0.3,
            size = 3.6) +
  scale_fill_manual(values = pal_hillside) +
  scale_y_continuous(
    limits = c(0, max_site_rich * 1.05),
    expand = c(0, 0)
  ) +
  labs(
    x = "Candidate site",
    y = "Species richness"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.ticks  = element_line(colour = "black")
  )

# ---- AOI Overall Richness Plot ----
max_aoi_rich <- max(aoi_rich_overall$richness, na.rm = TRUE)

aoirichness <- ggplot(aoi_rich_overall,
                      aes(hillside, richness, fill = hillside)) +
  geom_col(show.legend = FALSE, colour = "white") +
  geom_text(aes(label = richness),
            vjust = -0.3,
            size = 3.6) +
  scale_fill_manual(values = pal_hillside) +
  scale_y_continuous(
    limits = c(0, max_aoi_rich * 1.05),
    expand = c(0, 0)
  ) +
  labs(
    x = "Candidate site",
    y = "Species richness"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.ticks  = element_line(colour = "black")
  )

# ---- AOI Richness by Taxonomic Group (log10, with Total) ----
aoi_total <- aoi_group_rich %>%
  group_by(hillside) %>%
  summarise(group = "Total", richness = sum(richness), .groups = "drop")

aoi_group_rich_total <- bind_rows(aoi_group_rich, aoi_total) %>%
  mutate(group = factor(group, levels = c(unique(aoi_group_rich$group), "Total")))

aoigrouprichness <- ggplot(aoi_group_rich_total,
                           aes(group, richness, fill = hillside)) +
  geom_col(position = position_dodge(width = 0.75),
           width   = 0.65,
           colour  = "white") +
  geom_text(aes(label = richness),
            position = position_dodge(width = 0.75),
            vjust    = -0.3,
            size     = 3.6) +
  scale_fill_manual(values = pal_hillside) +
  scale_y_log10(
    breaks = c(1, 3, 10, 30, 100, 300, 600, 1000),
    labels = c("1","3","10","30","100","300","600","1000"),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = "Group",
    y = "Species richness (log10 scale)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.ticks  = element_line(colour = "black")
  )

aoigrouprichness
# ---- Candidate Site Richness by Group (with Total) ----
site_total <- sites_group %>%
  group_by(hillside) %>%
  summarise(group = "Total", richness = sum(richness), .groups = "drop")

sites_group_total <- bind_rows(sites_group, site_total) %>%
  mutate(group = factor(group, levels = c(unique(sites_group$group), "Total")))

max_site_group <- max(sites_group_total$richness, na.rm = TRUE)

candidatesiterichness <- ggplot(sites_group_total,
                                aes(group, richness, fill = hillside)) +
  geom_col(position = position_dodge(width = 0.75),
           width   = 0.65,
           colour  = "white") +
  geom_text(aes(label = richness),
            position = position_dodge(width = 0.75),
            vjust    = -0.25,
            size     = 3.6) +
  scale_fill_manual(values = pal_hillside) +
  scale_y_continuous(
    limits = c(0, max_site_group * 1.15),
    expand = c(0, 0)
  ) +
  labs(
    x = "Group",
    y = "Species richness"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.ticks  = element_line(colour = "black")
  )

candidatesiterichness
# ---- BoCC5 Classifications: Bird Tables ----
bocc5_part1 <- tribble(
  ~common_name, ~scientific_name, ~bocc_status, ~order, ~family,
  # RED list (subset)
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
  # AMBER (subset)
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
  "brambling", "fringilla montifringilla", "Amber", "Passeriformes", "Fringillidae",
  "jackdaw", "corvus monedula", "Amber", "Passeriformes", "Corvidae",
  "rook", "corvus frugilegus", "Amber", "Passeriformes", "Corvidae",
  "carrion crow", "corvus corone", "Amber", "Passeriformes", "Corvidae",
  "raven", "corvus corax", "Amber", "Passeriformes", "Corvidae",
  "magpie", "pica pica", "Amber", "Passeriformes", "Corvidae",
  "jay", "garrulus glandarius", "Amber", "Passeriformes", "Corvidae",
  "wood pigeon", "columba palumbus", "Amber", "Columbiformes", "Columbidae",
  "stock dove", "columba oenas", "Amber", "Columbiformes", "Columbidae",
  "collared dove", "streptopelia decaocto", "Amber", "Columbiformes", "Columbidae",
  "lapwing", "vanellus vanellus", "Amber", "Charadriiformes", "Charadriidae",
  "oystercatcher", "haematopus ostralegus", "Amber", "Charadriiformes", "Haematopodidae",
  "curlew", "numenius arquata", "Amber", "Charadriiformes", "Scolopacidae",
  "snipe", "gallinago gallinago", "Amber", "Charadriiformes", "Scolopacidae",
  "woodcock", "scolopax rusticola", "Amber", "Charadriiformes", "Scolopacidae",
  "herring gull", "larus argentatus", "Amber", "Charadriiformes", "Laridae",
  "kittiwake", "rissa tridactyla", "Amber", "Charadriiformes", "Laridae",
  "common gull", "larus canus", "Red", "Charadriiformes", "Laridae"
)

bocc5_part2 <- tribble(
  ~common_name, ~scientific_name, ~bocc_status, ~order, ~family,
  "black-headed gull", "chroicocephalus ridibundus", "Amber", "Charadriiformes", "Laridae",
  "lesser black-backed gull", "larus fuscus", "Amber", "Charadriiformes", "Laridae",
  "great black-backed gull", "larus marinus", "Red", "Charadriiformes", "Laridae",
  "sandwich tern", "thalasseus sandvicensis", "Amber", "Charadriiformes", "Laridae",
  "common tern", "sterna hirundo", "Amber", "Charadriiformes", "Laridae",
  "arctic tern", "sterna paradisaea", "Amber", "Charadriiformes", "Laridae",
  "mallard", "anas platyrhynchos", "Amber", "Anseriformes", "Anatidae",
  "teal", "anas crecca", "Amber", "Anseriformes", "Anatidae",
  "wigeon", "mareca penelope", "Amber", "Anseriformes", "Anatidae",
  "pintail", "anas acuta", "Amber", "Anseriformes", "Anatidae",
  "shoveler", "spatula clypeata", "Amber", "Anseriformes", "Anatidae",
  "pochard", "aythya ferina", "Amber", "Anseriformes", "Anatidae",
  "tufted duck", "aythya fuligula", "Amber", "Anseriformes", "Anatidae",
  "goldeneye", "bucephala clangula", "Amber", "Anseriformes", "Anatidae",
  "goosander", "mergus merganser", "Amber", "Anseriformes", "Anatidae",
  "pink-footed goose", "anser brachyrhynchus", "Amber", "Anseriformes", "Anatidae",
  "greylag goose", "anser anser", "Amber", "Anseriformes", "Anatidae",
  "canada goose", "branta canadensis", "Amber", "Anseriformes", "Anatidae",
  "barnacle goose", "branta leucopsis", "Amber", "Anseriformes", "Anatidae",
  "buzzard", "buteo buteo", "Amber", "Accipitriformes", "Accipitridae",
  "red kite", "milvus milvus", "Amber", "Accipitriformes", "Accipitridae",
  "tawny owl", "strix aluco", "Amber", "Strigiformes", "Strigidae",
  "long-eared owl", "asio otus", "Amber", "Strigiformes", "Strigidae",
  "short-eared owl", "asio flammeus", "Amber", "Strigiformes", "Strigidae",
  "wren", "troglodytes troglodytes", "Amber", "Passeriformes", "Troglodytidae",
  "great spotted woodpecker", "dendrocopos major", "Amber", "Piciformes", "Picidae",
  "green woodpecker", "picus viridis", "Amber", "Piciformes", "Picidae",
  "lesser spotted woodpecker", "dryobates minor", "Amber", "Piciformes", "Picidae",
  "grey partridge", "perdix perdix", "Amber", "Galliformes", "Phasianidae",
  "grey heron", "ardea cinerea", "Amber", "Pelecaniformes", "Ardeidae",
  "little egret", "egretta garzetta", "Amber", "Pelecaniformes", "Ardeidae",
  "great crested grebe", "podiceps cristatus", "Amber", "Podicipediformes", "Podicipedidae",
  "little grebe", "tachybaptus ruficollis", "Amber", "Podicipediformes", "Podicipedidae",
  "swallow", "hirundo rustica", "Amber", "Passeriformes", "Hirundinidae",
  "house martin", "delichon urbicum", "Amber", "Passeriformes", "Hirundinidae",
  "sand martin", "riparia riparia", "Amber", "Passeriformes", "Hirundinidae",
  "chiffchaff", "phylloscopus collybita", "Amber", "Passeriformes", "Phylloscopidae",
  "willow warbler", "phylloscopus trochilus", "Amber", "Passeriformes", "Phylloscopidae",
  "blackcap", "sylvia atricapilla", "Amber", "Passeriformes", "Sylviidae",
  "garden warbler", "sylvia borin", "Amber", "Passeriformes", "Sylviidae",
  "whitethroat", "sylvia communis", "Amber", "Passeriformes", "Sylviidae",
  "lesser whitethroat", "curruca curruca", "Amber", "Passeriformes", "Sylviidae",
  "treecreeper", "certhia familiaris", "Amber", "Passeriformes", "Certhiidae",
  "kingfisher", "alcedo atthis", "Amber", "Coraciiformes", "Alcedinidae",
  "dipper", "cinclus cinclus", "Amber", "Passeriformes", "Cinclidae",
  # Green list (subset)
  "robin", "erithacus rubecula", "Green", "Passeriformes", "Muscicapidae",
  "blue tit", "cyanistes caeruleus", "Green", "Passeriformes", "Paridae",
  "great tit", "parus major", "Green", "Passeriformes", "Paridae",
  "coal tit", "periparus ater", "Green", "Passeriformes", "Paridae",
  "long-tailed tit", "aegithalos caudatus", "Green", "Passeriformes", "Aegithalidae",
  "chaffinch", "fringilla coelebs", "Green", "Passeriformes", "Fringillidae",
  "goldfinch", "carduelis carduelis", "Green", "Passeriformes", "Fringillidae",
  "siskin", "spinus spinus", "Green", "Passeriformes", "Fringillidae",
  "greenfinch", "chloris chloris", "Green", "Passeriformes", "Fringillidae",
  "wren", "troglodytes troglotes", "Green", "Passeriformes", "Troglodytidae",
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
  "goldcrest", "regulus regulus", "Green", "Passeriformes", "Regulidae",
  "pheasant", "phasianus colchicus", "Green", "Galliformes", "Phasianidae",
  "moorhen", "gallinula chloropus", "Green", "Gruiformes", "Rallidae",
  "coot", "fulica atra", "Green", "Gruiformes", "Rallidae",
  "mallard", "anas platyrhynchos", "Green", "Anseriformes", "Anatidae",
  "tufted duck", "aythya fuligula", "Green", "Anseriformes", "Anatidae",
  "mute swan", "cygnus olor", "Green", "Anseriformes", "Anatidae",
  "grey heron", "ardea cinerea", "Green", "Pelecaniformes", "Ardeidae",
  "little egret", "egretta garzetta", "Green", "Pelecaniformes", "Ardeidae",
  "great crested grebe", "podiceps cristatus", "Green", "Podicipediformes", "Podicipedidae",
  "redshank", "tringa totanus", "Green", "Charadriiformes", "Scolopacidae",
  "common sandpiper", "actitis hypoleucos", "Green", "Charadriiformes", "Scolopacidae",
  "gannet", "morus bassanus", "Green", "Suliformes", "Sulidae",
  "great black-backed gull", "larus marinus", "Green", "Charadriiformes", "Laridae",
  "lesser black-backed gull", "larus fuscus", "Green", "Charadriiformes", "Laridae",
  "black-headed gull", "chroicocephalus ridibundus", "Green", "Charadriiformes", "Laridae",
  "oystercatcher", "haematopus ostralegus", "Green", "Charadriiformes", "Haematopodidae",
  "curlew", "numenius arquata", "Green", "Charadriiformes", "Scolopacidae",
  "robin", "erithacus rubecula", "Green", "Passeriformes", "Muscicapidae",
  "wren", "troglodytes troglotes", "Green", "Passeriformes", "Troglodytidae",
  "great tit", "parus major", "Green", "Passeriformes", "Paridae",
  "chaffinch", "fringilla coelebs", "Green", "Passeriformes", "Fringillidae",
  "goldfinch", "carduelis carduelis", "Green", "Passeriformes", "Fringillidae",
  "kestrel", "falco tinnunculus", "Green", "Falconiformes", "Falconidae",
  "peregrine falcon", "falco peregrinus", "Green", "Falconiformes", "Falconidae",
  "wood pigeon", "columba palumbus", "Green", "Columbiformes", "Columbidae",
  "feral pigeon", "columba livia domestica", "Green", "Columbiformes", "Columbidae",
  "collared dove", "streptopelia decaocto", "Green", "Columbiformes", "Columbidae",
  "grey wagtail", "motacilla cinerea", "Amber", "Passeriformes", "Motacillidae",
  "reed bunting", "emberiza schoeniclus", "Amber", "Passeriformes", "Emberizidae",
  "sedge warbler", "acrocephalus schoenobaenus", "Amber", "Passeriformes", "Acrocephalidae",
  "grasshopper warbler", "locustella naevia", "Amber", "Passeriformes", "Locustellidae",
  "grey partridge", "perdix perdix", "Red", "Galliformes", "Phasianidae",
  "ring ouzel", "turdus torquatus", "Red", "Passeriformes", "Turdidae",
  "yellowhammer", "emberiza citrinella", "Red", "Passeriformes", "Emberizidae",
  "house sparrow", "passer domesticus", "Green", "Passeriformes", "Passeridae",
  NA, NA, NA, NA, NA
) %>%
  filter(!is.na(common_name))

# ---- Compile and Harmonise BoCC5 ----
bocc5 <- bind_rows(bocc5_part1, bocc5_part2, bocc5_part3) %>%
  distinct(common_name, .keep_all = TRUE) %>%
  mutate(
    common_name     = tolower(common_name),
    scientific_name = tolower(scientific_name),
    bocc_status     = tolower(bocc_status)
  ) %>%
  mutate(
    bocc_status = case_when(
      # Raptors/owls headline fixes
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
      common_name %in% c("pink-footed goose","barnacle goose") ~ "amber",
      common_name %in% c("greylag goose","canada goose") ~ "green",
      # Common passerines to Green
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

bocc5_for_importance <- bocc5 %>%
  filter(bocc_status %in% c("Red","Amber"))

# ---- Mammals: Protected Species Lists ----
mammals_imp <- tribble(
  ~common_name,            ~scientific_name,            ~status,                    ~order,         ~family,
  "otter",                 "lutra lutra",               "EPS (fully protected)",    "Carnivora",    "Mustelidae",
  "common pipistrelle",    "pipistrellus pipistrellus", "EPS (fully protected)",    "Chiroptera",   "Vespertilionidae",
  "soprano pipistrelle",   "pipistrellus pygmaeus",     "EPS (fully protected)",    "Chiroptera",   "Vespertilionidae",
  "nathusius' pipistrelle","pipistrellus nathusii",     "EPS (fully protected)",    "Chiroptera",   "Vespertilionidae",
  "lesser noctule",        "nyctalus leisleri",         "EPS (fully protected)",    "Chiroptera",   "Vespertilionidae",
  "brown long-eared bat",  "plecotus auritus",          "EPS (fully protected)",    "Chiroptera",   "Vespertilionidae",
  "badger",                "meles meles",               "Protected (PBA 1992)",     "Carnivora",    "Mustelidae",
  "red squirrel",          "sciurus vulgaris",          "Protected (WCA S5), SBL",  "Rodentia",     "Sciuridae",
  "pine marten",           "martes martes",             "Protected (WCA S5)",       "Carnivora",    "Mustelidae",
  "water vole",            "arvicola amphibius",        "Protected (WCA S5)",       "Rodentia",     "Cricetidae",
  "wildcat",               "felis silvestris",          "EPS (fully protected)",    "Carnivora",    "Felidae",
  "mountain hare",         "lepus timidus",             "Protected (Scotland)",     "Lagomorpha",   "Leporidae",
  "beaver",                "castor fiber",              "EPS (fully protected)",    "Rodentia",     "Castoridae",
  "hedgehog",              "erinaceus europaeus",       "SBL Priority",             "Eulipotyphla", "Erinaceidae",
  "brown hare",            "lepus europaeus",           "SBL Priority",             "Lagomorpha",   "Leporidae"
)

mammals_imp_extra <- tribble(
  ~common_name,            ~scientific_name,           ~status,                       ~order,          ~family,
  "daubenton's bat",       "myotis daubentonii",       "EPS (fully protected)",        "Chiroptera",    "Vespertilionidae",
  "natterer's bat",        "myotis nattereri",         "EPS (fully protected)",        "Chiroptera",    "Vespertilionidae",
  "whiskered bat",         "myotis mystacinus",        "EPS (fully protected)",        "Chiroptera",    "Vespertilionidae",
  "brandt's bat",          "myotis brandtii",          "EPS (fully protected)",        "Chiroptera",    "Vespertilionidae",
  "red deer",              "cervus elaphus",           "None",                         "Artiodactyla",  "Cervidae",
  "roe deer",              "capreolus capreolus",      "None",                         "Artiodactyla",  "Cervidae",
  "common seal",           "phoca vitulina",           "WCA S5 (disturbance)",         "Carnivora",     "Phocidae",
  "grey seal",             "halichoerus grypus",       "WCA S5 (disturbance)",         "Carnivora",     "Phocidae",
  "harbour porpoise",      "phocoena phocoena",        "EPS (fully protected)",        "Cetacea",       "Phocoenidae",
  "weasel",                "mustela nivalis",          "None",                         "Carnivora",     "Mustelidae",
  "stoat",                 "mustela erminea",          "None",                         "Carnivora",     "Mustelidae",
  "wood mouse",            "apodemus sylvaticus",      "None",                         "Rodentia",      "Muridae",
  "bank vole",             "myodes glareolus",         "None",                         "Rodentia",      "Cricetidae",
  "common shrew",          "sorex araneus",            "None",                         "Eulipotyphla",  "Soricidae"
)

# ---- Harmonise and Combine Mammals ----
mammals_imp_extra <- mammals_imp_extra %>%
  mutate(
    common_name     = tolower(common_name),
    scientific_name = tolower(scientific_name),
    group           = "Mammal",
    bocc_status     = NA_character_
  )

mammals_imp <- mammals_imp %>%
  mutate(
    common_name     = tolower(common_name),
    scientific_name = tolower(scientific_name),
    group           = "Mammal",
    bocc_status     = NA_character_
  )

mammals_imp <- bind_rows(mammals_imp, mammals_imp_extra) %>%
  distinct(common_name, .keep_all = TRUE) %>%
  mutate(
    status = case_when(
      str_detect(status, regex("EPS", ignore_case = TRUE)) ~ "EPS",
      str_detect(status, regex("WCA|PBA|Protected \\(Scotland\\)", ignore_case = TRUE)) ~ "Protected",
      str_detect(status, regex("SBL|Priority", ignore_case = TRUE)) ~ "Priority",
      TRUE ~ NA_character_
    )
  )

# ---- Important Species Lookup (Birds + Mammals) ----
important_birds <- bocc5_for_importance %>%
  mutate(
    common_name     = tolower(str_squish(common_name)),
    scientific_name = tolower(str_squish(scientific_name)),
    status          = bocc_status,
    group           = "Bird"
  ) %>%
  select(common_name, scientific_name, status, group, order, family)

important_mammals <- mammals_imp %>%
  select(common_name, scientific_name, status, group, order, family)

lookup_important <- bind_rows(important_birds, important_mammals) %>%
  distinct(scientific_name, .keep_all = TRUE)

# ---- Normalise Scientific Names for Matching ----
clean_sci <- function(x) {
  x %>%
    tolower() %>%
    str_squish() %>%
    str_replace("\\s+agg\\.?$", "") %>%
    str_replace("\\s+subsp\\..*", "") %>%
    str_replace("\\s+var\\..*", "") %>%
    str_replace("\\s*\\(.*?\\)", "") %>%
    str_replace("=.+", "") %>%
    str_replace_all("/", " ") %>%
    str_replace_all("[^a-z\\s]", " ") %>%
    str_squish() %>%
    word(1, 2)
}

# ---- Apply Cleaning to AOI Buffers + Lookup ----
buffers_clean <- buffers %>%
  mutate(
    scientific_name_raw = scientific_name,
    scientific_name      = clean_sci(scientific_name_raw)
  )

lookup_important_clean <- lookup_important %>%
  mutate(
    scientific_name_raw = scientific_name,
    scientific_name      = clean_sci(scientific_name_raw)
  )

# ---- AOI: Important Species Join ----
important_nbn_aoi <- buffers_clean %>%
  left_join(
    lookup_important_clean,
    by = "scientific_name",
    suffix = c("_buf", "_imp")
  ) %>%
  filter(!is.na(status)) %>%
  mutate(group = group_imp) %>%
  distinct(hillside, scientific_name_raw_buf, .keep_all = TRUE) %>%
  arrange(hillside, group, scientific_name_raw_buf)

# ---- AOI: Important Species Summaries ----
important_richness_aoi <- important_nbn_aoi %>%
  count(hillside, name = "n_important_species")

important_by_status_aoi <- important_nbn_aoi %>%
  count(hillside, status, name = "n_species")

important_by_group_aoi <- important_nbn_aoi %>%
  count(hillside, group, name = "n_species")

# ---- AOI: Important Species Plots ----
max_imp_rich_aoi <- max(important_richness_aoi$n_important_species, na.rm = TRUE)

p_imp1_aoi <- ggplot(important_richness_aoi,
                     aes(hillside, n_important_species, fill = hillside)) +
  geom_col(show.legend = FALSE, colour = "white") +
  geom_text(aes(label = n_important_species),
            vjust = -0.3,
            size  = 3.6) +
  scale_fill_manual(values = pal_hillside) +
  scale_y_continuous(limits = c(0, max_imp_rich_aoi * 1.05), expand = c(0, 0)) +
  labs(
    title = "Important species richness (AOI)",
    x = "Candidate site",
    y = "Number of important species"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

max_imp_status_aoi <- max(important_by_status_aoi$n_species, na.rm = TRUE)

p_imp_status_aoi <- ggplot(important_by_status_aoi,
                           aes(status, n_species, fill = hillside)) +
  geom_col(position = position_dodge(width = 0.75),
           width   = 0.65,
           colour  = "white") +
  geom_text(aes(label = n_species),
            position = position_dodge(width = 0.75),
            vjust    = -0.3,
            size     = 3.6) +
  scale_fill_manual(values = pal_hillside) +
  scale_y_continuous(limits = c(0, max_imp_status_aoi * 1.05), expand = c(0, 0)) +
  labs(
    x = "Status",
    y = "Number of species"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank(),
    axis.line   = element_line(colour = "black"),
    axis.ticks  = element_line(colour = "black")
  )

# ---- AOI: Important Species by Group (Birds, Bats, Other mammals, Total) ----
important_clean_aoi <- important_nbn_aoi %>%
  mutate(
    sci_clean = clean_sci(scientific_name_raw_buf),
    is_bat    = str_detect(sci_clean, "pipistrell|myotis|plecotus|nyctalus|eptesicus"),
    is_mammal = group == "Mammal",
    is_bird   = group == "Bird"
  )

important_by_group_fixed_aoi <- important_clean_aoi %>%
  mutate(
    group_final = case_when(
      is_bat            ~ "Bats",
      is_mammal         ~ "Other mammals",
      is_bird           ~ "Birds",
      TRUE              ~ NA_character_
    )
  ) %>%
  filter(!is.na(group_final)) %>%
  count(hillside, group_final, name = "n_species")

important_totals_aoi <- important_by_group_fixed_aoi %>%
  group_by(hillside) %>%
  summarise(n_species = sum(n_species), .groups = "drop") %>%
  mutate(group_final = "Total")

important_by_group_all_aoi <- bind_rows(important_by_group_fixed_aoi, important_totals_aoi) %>%
  mutate(
    group_final = factor(group_final, levels = c("Birds", "Bats", "Other mammals", "Total"))
  )

max_imp_group_all_aoi <- max(important_by_group_all_aoi$n_species, na.rm = TRUE)

p_imp_group_aoi <- ggplot(
  important_by_group_all_aoi,
  aes(group_final, n_species, fill = hillside)
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width    = 0.65,
    colour   = "white"
  ) +
  geom_text(
    aes(label = n_species),
    position = position_dodge(width = 0.75),
    vjust    = -0.3,
    size     = 3.6
  ) +
  scale_fill_manual(values = pal_hillside) +
  scale_y_continuous(
    limits = c(0, max_imp_group_all_aoi * 1.12),
    expand = c(0, 0)
  ) +
  labs(
    x = "Group",
    y = "Species of concern richness"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid  = element_blank(),
    axis.line   = element_line(colour = "black"),
    axis.ticks  = element_line(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ---- AOI: Common Names for Heatmaps ----
important_clean_aoi <- important_nbn_aoi %>%
  mutate(
    common_name_buf = tolower(str_squish(common_name))
  ) %>%
  left_join(
    lookup_important_clean %>%
      select(scientific_name, common_name) %>%
      rename(common_name_imp = common_name),
    by = "scientific_name"
  ) %>%
  mutate(
    common_name_final = coalesce(common_name_imp, common_name_buf),
    common_name_final = tools::toTitleCase(common_name_final),
    sci_clean         = clean_sci(scientific_name_raw_buf),
    is_bat            = str_detect(sci_clean, "pipistrell|myotis|plecotus|nyctalus|eptesicus"),
    is_mammal         = group == "Mammal",
    is_bird           = group == "Bird"
  )

birds_classified_aoi <- important_clean_aoi %>%
  filter(is_bird) %>%
  mutate(
    classification  = status,
    species_display = common_name_final
  ) %>%
  filter(!is.na(classification))

important_clean_aoi <- important_clean_aoi %>%
  mutate(
    status_clean = case_when(
      str_detect(status, regex("EPS",        ignore_case = TRUE)) ~ "EPS",
      str_detect(status, regex("Protected",  ignore_case = TRUE)) ~ "Protected",
      str_detect(status, regex("Priority",   ignore_case = TRUE)) ~ "Priority",
      status %in% c("Red", "Amber", "Green")                      ~ status,
      TRUE ~ NA_character_
    )
  )

mammals_classified_aoi <- important_clean_aoi %>%
  filter(is_mammal | is_bat) %>%
  mutate(
    classification  = status_clean,
    species_display = common_name_final
  ) %>%
  filter(!is.na(classification))

pal_all <- c(
  "EPS"       = "#1B9E77",
  "Protected" = "#7570B3",
  "Priority"  = "#E7298A",
  "Red"       = "#C1121F",
  "Amber"     = "#FFB000",
  "Green"     = "#4C6FA3",
  "Other"     = "grey80"
)

fig_birds_heatmap_aoi <- ggplot(
  birds_classified_aoi,
  aes(
    x    = hillside,
    y    = factor(species_display),
    fill = classification
  )
) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(values = pal_all) +
  labs(
    x = "Area of interest",
    y = "Species (common names)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

fig_mammals_heatmap_aoi <- ggplot(
  mammals_classified_aoi,
  aes(
    x    = hillside,
    y    = factor(species_display),
    fill = classification
  )
) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(values = pal_all) +
  labs(
    x = "Area of interest",
    y = "Species (common names)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )
# Packages (if not already loaded)
library(tidyverse)
library(stringr)
library(ggplot2)

# ---- Candidate Site: Important Species Join ----
sites_clean <- sites %>%
  mutate(
    scientific_name_raw = scientific_name,
    scientific_name     = clean_sci(scientific_name_raw)
  )

important_site <- sites_clean %>%
  left_join(
    lookup_important_clean %>%
      rename(
        group_imp        = group,
        status_imp       = status,
        common_name_imp  = common_name
      ),
    by = "scientific_name"
  ) %>%
  filter(!is.na(status_imp)) %>%
  mutate(
    group = group_imp,
    status = status_imp,
    common_name_final = tools::toTitleCase(
      str_squish(coalesce(common_name_imp, vernacular_name, scientific_name))
    )
  ) %>%
  distinct(hillside, scientific_name, .keep_all = TRUE) %>%  # de-duplicate per hillside x species
  arrange(hillside, group, scientific_name)

# ---- Candidate Site: Important Species Summaries ----
important_richness_site <- important_site %>%
  count(hillside, name = "n_important_species")

important_by_status_site <- important_site %>%
  count(hillside, status, name = "n_species")

important_by_group_site_raw <- important_site %>%
  count(hillside, group, name = "n_species")

# ---- Candidate Site: Important Species Plots ----
max_imp_rich_site <- max(important_richness_site$n_important_species, na.rm = TRUE)

p_imp1_site <- ggplot(important_richness_site,
                      aes(hillside, n_important_species, fill = hillside)) +
  geom_col(show.legend = FALSE, colour = "white") +
  geom_text(aes(label = n_important_species),
            vjust = -0.3,
            size  = 3.6) +
  scale_fill_manual(values = pal_hillside) +
  scale_y_continuous(limits = c(0, max_imp_rich_site * 1.05), expand = c(0, 0)) +
  labs(
    x = "Candidate site",
    y = "Species of concern richness"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

max_imp_status_site <- max(important_by_status_site$n_species, na.rm = TRUE)

p_imp_status_site <- ggplot(important_by_status_site,
                            aes(status, n_species, fill = hillside)) +
  geom_col(position = position_dodge(width = 0.75),
           width   = 0.65,
           colour  = "white") +
  geom_text(aes(label = n_species),
            position = position_dodge(width = 0.75),
            vjust    = -0.3,
            size     = 3.6) +
  scale_fill_manual(values = pal_hillside) +
  scale_y_continuous(limits = c(0, max_imp_status_site * 1.05), expand = c(0, 0)) +
  labs(
    x = "Status",
    y = "Number of species"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank(),
    axis.line   = element_line(colour = "black"),
    axis.ticks  = element_line(colour = "black")
  )

# ---- Candidate Site: Important Species by Group (Birds, Bats, Other mammals, Total) ----
# FIX: use existing columns (scientific_name, scientific_name_raw, common_name_final)
# ---- Candidate Site: Important Species by Group (Birds, Bats, Other mammals, Total) ----
important_clean_site <- important_site %>%
  mutate(
    # Use the already-cleaned scientific_name; don't reference scientific_name_raw
    sci_clean = str_to_lower(scientific_name),
    is_bat    = str_detect(sci_clean, "pipistrell|myotis|plecotus|nyctalus|eptesicus"),
    is_mammal = group == "Mammal",
    is_bird   = group == "Bird",
    group_fixed = case_when(
      is_bat    ~ "Bats",
      is_mammal ~ "Other mammals",
      is_bird   ~ "Birds",
      TRUE      ~ NA_character_
    ),
    group_fixed = factor(group_fixed, levels = c("Birds", "Bats", "Other mammals")),
    # Reuse and tidy the name created earlier
    common_name_final = tools::toTitleCase(str_squish(common_name_final))
  )


important_by_group_site <- important_clean_site %>%
  filter(!is.na(group_fixed)) %>%
  count(hillside, group_fixed, name = "n_species")

important_totals_site <- important_clean_site %>%
  filter(!is.na(group_fixed)) %>%
  count(hillside, name = "n_species") %>%
  mutate(group_fixed = "Total")

important_by_group_site_all <- bind_rows(
  important_by_group_site,
  important_totals_site
) %>%
  mutate(
    group_fixed = factor(group_fixed,
                         levels = c("Birds", "Bats", "Other mammals", "Total"))
  )

max_imp_group_site <- max(important_by_group_site_all$n_species, na.rm = TRUE)

p_imp_group_site <- ggplot(
  important_by_group_site_all,
  aes(group_fixed, n_species, fill = hillside)
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width    = 0.65,
    colour   = "white"
  ) +
  geom_text(
    aes(label = n_species),
    position = position_dodge(width = 0.75),
    vjust    = -0.25,
    size     = 3.6
  ) +
  scale_fill_manual(values = pal_hillside) +
  scale_y_continuous(limits = c(0, max_imp_group_site * 1.12), expand = c(0, 0)) +
  labs(
    x = "Group",
    y = "Species of concern richness"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid  = element_blank(),
    axis.line   = element_line(colour = "black"),
    axis.ticks  = element_line(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ---- Candidate Site: Heatmaps (Common Names) ----
important_clean_site <- important_clean_site %>%
  mutate(
    status_clean = case_when(
      str_detect(status, regex("EPS",        ignore_case = TRUE)) ~ "EPS",
      str_detect(status, regex("Protected",  ignore_case = TRUE)) ~ "Protected",
      str_detect(status, regex("Priority",   ignore_case = TRUE)) ~ "Priority",
      status %in% c("Red", "Amber", "Green")                      ~ status,
      TRUE ~ NA_character_
    )
  )

birds_classified_site <- important_clean_site %>%
  filter(group_fixed == "Birds") %>%
  mutate(
    classification  = status_clean,
    species_display = common_name_final
  ) %>%
  filter(!is.na(classification))

mammals_classified_site <- important_clean_site %>%
  filter(group_fixed %in% c("Bats", "Other mammals")) %>%
  mutate(
    classification  = status_clean,
    species_display = common_name_final
  ) %>%
  filter(!is.na(classification))

# Remove legend AND remove x-axis label for birds
fig_birds_nokey <- fig_birds_heatmap_aoi +
  theme(legend.position = "none") +
  xlab(NULL)

# Mammal plot: remove legend only
fig_mammals_nokey <- fig_mammals_heatmap_aoi +
  theme(legend.position = "none")

# Combine with proportional heights
fig_birds_mammals_heatmap_aoi <-
  plot_grid(
    fig_birds_nokey,
    fig_mammals_nokey,
    ncol = 1,
    align = "v",
    axis = "lr",
    rel_heights = c(n_birds, n_mammals)
  )

# Overall richness
candidaterichness
aoirichness

# AOI important species (3 plots)
p_imp1_aoi
p_imp_status_aoi

# Candidate site important species (3 plots)
p_imp1_site
p_imp_status_site


# ---- Print key plots ----


# Richness by group
aoigrouprichness
candidatesiterichness

# AOI important species (3 plots)
p_imp_group_aoi

# Candidate site important species (3 plots)
p_imp_group_site

# AOI heatmaps
fig_birds_heatmap_aoi
fig_mammals_heatmap_aoi

# Candidate site heatmaps
p_imp_birds_candidate
p_imp_mammals_candidate

#----Heatmap panelling----
library(ggplot2)
library(cowplot)
library(dplyr)

## ---- Species counts ----
n_birds <- 47
n_mammals <- 15

## ---- Modify bird plot: remove legend, x-label, x-ticks, x-line ----
fig_birds_nokey <- fig_birds_heatmap_aoi +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),       # remove x-axis title
    axis.text.x = element_blank(),        # remove x tick labels
    axis.ticks.x = element_blank(),       # remove tick marks
    axis.line.x = element_blank()         # remove x axis line
  )

## ---- Modify mammal plot: remove legend only ----
fig_mammals_nokey <- fig_mammals_heatmap_aoi +
  theme(
    legend.position = "none"
  )

## ---- Combine with proportional heights ----
fig_birds_mammals_heatmap_aoi <- plot_grid(
  fig_birds_nokey,
  fig_mammals_nokey,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(n_birds, n_mammals)
)

fig_birds_mammals_heatmap_aoi

## ---- Save final figure ----
ggsave(
  filename = "fig_birds_mammals_heatmap_aoi.png",
  plot     = fig_birds_mammals_heatmap_aoi,
  width    = 5,
  height   = 9,
  dpi      = 1000
)



p_imp_birds_candidate
p_imp_mammals_candidate

# candidate site now
library(ggplot2)
library(cowplot)
library(dplyr)

## ---- Species counts ----
n_birds <- 47
n_mammals <- 7.5


## ---- Modify bird plot: remove legend, x-label, x-ticks, x-line ----
fig_birds_nokey <- p_imp_birds_candidate +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),       # remove x-axis title
    axis.text.x = element_blank(),        # remove x tick labels
    axis.ticks.x = element_blank(),       # remove tick marks
    axis.line.x = element_blank()         # remove x axis line
  )

## ---- Modify mammal plot: remove legend only ----
fig_mammals_nokey <- p_imp_mammals_candidate +
  theme(
    legend.position = "none"
  )

## ---- Get shared legend from one of the original plots ----
legend_candidate <- get_legend(
  p_imp_birds_candidate +
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal"
      # you can tweak legend.title/text sizes here if you like
    )
)

## ---- Combine with proportional heights (panels + legend) ----
fig_birds_mammals_heatmap_candidate <- plot_grid(
  fig_birds_nokey,
  fig_mammals_nokey,
  legend_candidate,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(n_birds, n_mammals, 5)   # small height for legend
)

fig_birds_mammals_heatmap_candidate

## ---- Save final figure ----
ggsave(
  filename = "fig_birds_mammals_heatmap_candidate_site.png",
  plot     = fig_birds_mammals_heatmap_candidate,
  width    = 5,
  height   = 6,
  dpi      = 1000
)

