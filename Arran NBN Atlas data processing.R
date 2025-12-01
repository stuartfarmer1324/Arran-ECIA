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
  
  if (!"species_name" %in% names(df)) df$species_name <- NA
  if (!"vernacular_name" %in% names(df)) df$vernacular_name <- NA
  if (!"kingdom" %in% names(df)) df$kingdom <- NA
  if (!"family" %in% names(df)) df$family <- NA
  if (!"genus" %in% names(df)) df$genus <- NA
  
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
        # Birds (broad-ish)
        str_detect(sci, "anas|aegithalos|accipiter|acrocephalus|actitis|alauda|anthus|emberiza|falco|muscicapa|prunella|strix|turdus|corvus|motacilla|larus|rissa|sterna|haematopus|numenius|podiceps|ardea|egretta|hirundo|delichon|riparia") ~ "Birds",
        
        # Bats
        str_detect(sci, "pipistrell|myotis|plecotus|nyctalus|eptesicus") ~ "Bats",
        
        # Mammals (non-bat)
        str_detect(sci, "sciurus|vulpes|mustela|meles|oryctolagus|capreolus|sus|martes|arvicola|lutra|castor|felis|erinaceus|lepus|cervus|phoca|halichoerus|phocoena|apodemus|myodes|sorex") ~ "Mammals",
        
        # Amphibians
        str_detect(sci, "bufo|rana|lissotriton|ichthyosaura|triturus") ~ "Amphibians",
        fam == "bufonidae" ~ "Amphibians",
        
        # Reptiles
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

# ---- Read Raw Buffer Data ----
north500_raw  <- read_nbn(north_500_path,  "North")
south500_raw  <- read_nbn(south_500_path,  "South")
north2000_raw <- read_nbn(north_2000_path, "North")
south2000_raw <- read_nbn(south_2000_path, "South")

# ---- Filter Before Combining (zone rules) ----
# 500 m: exclude Birds, Bats, Plants (as you originally intended)
north500 <- assign_group(north500_raw) %>%
  filter(!group %in% c("Birds", "Bats", "Plants")) %>%
  mutate(zone = "500m", hillside = "North")

south500 <- assign_group(south500_raw) %>%
  filter(!group %in% c("Birds", "Bats", "Plants")) %>%
  mutate(zone = "500m", hillside = "South")

# 2000 m: INCLUDE Birds & Bats (key fix so bats exist in AOI)
north2000 <- assign_group(north2000_raw) %>%
  filter(group %in% c("Birds", "Bats")) %>%
  mutate(zone = "2000m", hillside = "North")

south2000 <- assign_group(south2000_raw) %>%
  filter(group %in% c("Birds", "Bats")) %>%
  mutate(zone = "2000m", hillside = "South")

# ---- Combine Cleaned Buffers ----
buffers <- bind_rows(north500, south500, north2000, south2000)

# ---- AOI Richness Summaries ----
aoi_rich_overall <- buffers %>%
  distinct(hillside, scientific_name) %>%
  count(hillside, name = "richness")

aoi_group_rich <- buffers %>%
  distinct(hillside, group, scientific_name) %>%
  count(hillside, group, name = "richness")

# ---- AOI Plots ----
p1 <- ggplot(site_richness, aes(hillside, richness, fill = hillside)) +
  geom_col(show.legend = FALSE) +
  theme_minimal() +
  labs(
    title = "Candidate Site Richness",
    x = "Hillside",
    y = "Species richness"
  )

p2 <- ggplot(aoi_rich_overall, aes(hillside, richness, fill = hillside)) +
  geom_col(show.legend = FALSE) +
  theme_minimal() +
  labs(
    title = "AOI Richness (500 m + 2000 m)",
    x = "Hillside",
    y = "Species richness"
  )

p3 <- ggplot(aoi_group_rich, aes(group, richness, fill = hillside)) +
  geom_col(position = "dodge") +
  scale_y_log10(
    breaks = c(1, 3, 10, 30, 100, 300, 600, 1000),
    labels = c("1","3","10","30","100","300","600","1000")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "AOI Richness by Taxonomic Group",
    x = "Group",
    y = "Species richness (log10 scale)"
  )

print(p1)
print(p2)
print(p3)

# ---- BoCC5 Classifications ----
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
  "wren", "troglodytes troglotes", "Green", "Passeriformes", "Troglodytidae",
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
      # Geese
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

# ---- Mammals: Protected Species List ----
mammals_imp <- tribble(
  ~common_name,            ~scientific_name,           ~status,                    ~order,         ~family,
  
  # EPS (Habitats Regulations)
  "otter",                 "lutra lutra",              "EPS (fully protected)",    "Carnivora",    "Mustelidae",
  "common pipistrelle",    "pipistrellus pipistrellus","EPS (fully protected)",   "Chiroptera",   "Vespertilionidae",
  "soprano pipistrelle",   "pipistrellus pygmaeus",    "EPS (fully protected)",    "Chiroptera",   "Vespertilionidae",
  "nathusius' pipistrelle","pipistrellus nathusii",    "EPS (fully protected)",    "Chiroptera",   "Vespertilionidae",
  "lesser noctule",        "nyctalus leisleri",        "EPS (fully protected)",    "Chiroptera",   "Vespertilionidae",
  "brown long-eared bat",  "plecotus auritus",         "EPS (fully protected)",    "Chiroptera",   "Vespertilionidae",
  
  # WCA S5 / Protected (Scotland)
  "badger",                "meles meles",              "Protected (PBA 1992)",     "Carnivora",    "Mustelidae",
  "red squirrel",          "sciurus vulgaris",         "Protected (WCA S5), SBL",  "Rodentia",     "Sciuridae",
  "pine marten",           "martes martes",            "Protected (WCA S5)",       "Carnivora",    "Mustelidae",
  "water vole",            "arvicola amphibius",       "Protected (WCA S5)",       "Rodentia",     "Cricetidae",
  "wildcat",               "felis silvestris",         "EPS (fully protected)",    "Carnivora",    "Felidae",
  "mountain hare",         "lepus timidus",            "Protected (Scotland)",     "Lagomorpha",   "Leporidae",
  "beaver",                "castor fiber",             "EPS (fully protected)",    "Rodentia",     "Castoridae",
  
  # SBL Priority
  "hedgehog",              "erinaceus europaeus",      "SBL Priority",             "Eulipotyphla", "Erinaceidae",
  "brown hare",            "lepus europaeus",          "SBL Priority",             "Lagomorpha",   "Leporidae"
)

# ---- Additional mammals for Arran ----
mammals_imp_extra <- tribble(
  ~common_name,            ~scientific_name,           ~status,                       ~order,          ~family,
  
  # ---- BATS ----
  "daubenton's bat",       "myotis daubentonii",       "EPS (fully protected)",        "Chiroptera",    "Vespertilionidae",
  "natterer's bat",        "myotis nattereri",         "EPS (fully protected)",        "Chiroptera",    "Vespertilionidae",
  "whiskered bat",         "myotis mystacinus",        "EPS (fully protected)",        "Chiroptera",    "Vespertilionidae",
  "brandt's bat",          "myotis brandtii",          "EPS (fully protected)",        "Chiroptera",    "Vespertilionidae",
  
  # ---- WIDER MAMMALS ----
  "red deer",              "cervus elaphus",           "None",                         "Artiodactyla",  "Cervidae",
  "roe deer",              "capreolus capreolus",      "None",                         "Artiodactyla",  "Cervidae",
  "common seal",           "phoca vitulina",           "WCA S5 (disturbance)",         "Carnivora",     "Phocidae",
  "grey seal",             "halichoerus grypus",       "WCA S5 (disturbance)",         "Carnivora",     "Phocidae",
  "harbour porpoise",      "phocoena phocoena",        "EPS (fully protected)",        "Cetacea",       "Phocoenidae",
  
  # Mustelids
  "weasel",                "mustela nivalis",          "None",                         "Carnivora",     "Mustelidae",
  "stoat",                 "mustela erminea",          "None",                         "Carnivora",     "Mustelidae",
  
  # Small mammals
  "wood mouse",            "apodemus sylvaticus",      "None",                         "Rodentia",      "Muridae",
  "bank vole",             "myodes glareolus",         "None",                         "Rodentia",      "Cricetidae",
  "common shrew",          "sorex araneus",            "None",                         "Eulipotyphla",  "Soricidae"
)

# ---- Harmonise and combine mammals ----
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
    # Normalise and collapse to EPS / Protected / Priority / NA (non-important)
    status = case_when(
      str_detect(status, regex("EPS", ignore_case = TRUE)) ~ "EPS",
      str_detect(status, regex("WCA|PBA|Protected \\(Scotland\\)", ignore_case = TRUE)) ~ "Protected",
      str_detect(status, regex("SBL|Priority", ignore_case = TRUE)) ~ "Priority",
      TRUE ~ NA_character_
    )
  )

# ---- Important Species Lookup ----
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
    # Remove 'agg.' etc.
    str_replace("\\s+agg\\.?$", "") %>%
    # Remove subsp/var/nothospecies
    str_replace("\\s+subsp\\..*", "") %>%
    str_replace("\\s+var\\..*", "") %>%
    # Remove parentheses and contents
    str_replace("\\s*\\(.*?\\)", "") %>%
    # Remove equals synonyms (=xxx)
    str_replace("=.+", "") %>%
    # Replace "/" with space
    str_replace_all("/", " ") %>%
    # Remove remaining punctuation
    str_replace_all("[^a-z\\s]", " ") %>%
    str_squish() %>%
    # Keep only genus + species
    word(1, 2)
}

# Apply to NBN buffers
buffers_clean <- buffers %>%
  mutate(
    scientific_name_raw = scientific_name,
    scientific_name      = clean_sci(scientific_name_raw)
  )

# Apply to lookup table
lookup_important_clean <- lookup_important %>%
  mutate(
    scientific_name_raw = scientific_name,
    scientific_name      = clean_sci(scientific_name_raw)
  )

# ---- Important NBN Join ----
important_nbn <- buffers_clean %>%
  left_join(
    lookup_important_clean,
    by = "scientific_name",
    suffix = c("_buf", "_imp")
  ) %>%
  # Only keep species with a non-NA status (Red/Amber birds, EPS/Protected/Priority mammals)
  filter(!is.na(status)) %>%
  mutate(group = group_imp) %>%
  distinct(hillside, scientific_name_raw_buf, .keep_all = TRUE) %>%
  arrange(hillside, group, scientific_name_raw_buf)

print(
  important_nbn %>%
    select(hillside, zone, scientific_name_raw_buf, common_name, group, status)
)

# ---- Important Species Summaries ----
important_richness <- important_nbn %>%
  count(hillside, name = "n_important_species")

important_by_status <- important_nbn %>%
  count(hillside, status, name = "n_species")

p_imp1 <- ggplot(important_richness, aes(hillside, n_important_species, fill = hillside)) +
  geom_col(show.legend = FALSE) +
  theme_minimal() +
  labs(
    title = "Important Species Richness (NBN AOI)",
    x = "Hillside",
    y = "Number of important species"
  )

p_imp2 <- ggplot(important_by_status, aes(status, n_species, fill = hillside)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Conservation Category of Important NBN Species",
    x = "Status",
    y = "Number of species"
  )

print(p_imp1)
print(p_imp2)

head(important_nbn)

# ---- FILTER ONLY IMPORTANT / PROTECTED MAMMALS & BATS ----
important_clean <- important_nbn %>%
  mutate(
    sci_clean = clean_sci(scientific_name_raw_buf),
    is_bat    = str_detect(sci_clean, "pipistrell|myotis|plecotus|nyctalus|eptesicus"),
    is_bird   = group == "Bird",
    is_mammal = group == "Mammal" & !is_bat
  )

# ---- 1. BIRDS ----
birds_classified <- important_clean %>%
  filter(is_bird) %>%
  left_join(
    bocc5 %>%
      mutate(scientific_name = clean_sci(scientific_name)) %>%
      select(scientific_name, bocc_status),
    by = c("sci_clean" = "scientific_name")
  ) %>%
  mutate(
    classification = case_when(
      bocc_status %in% c("Red","Amber","Green") ~ bocc_status,
      TRUE ~ "Other"
    ),
    species = scientific_name_raw_buf
  )

# ---- 2 + 3. ALL IMPORTANT / PROTECTED MAMMALS (including bats) ----
mammals_classified <- important_clean %>%
  filter(group == "Mammal") %>%   # <--- includes bats now too
  mutate(
    classification = case_when(
      status == "EPS"       ~ "EPS",
      status == "Protected" ~ "Protected",
      status == "Priority"  ~ "Priority",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(classification)) %>%  # keep only important mammals
  mutate(species = scientific_name_raw_buf)

# ---- COLOUR PALETTE ----
pal_all <- c(
  "EPS"       = "#2D6A4F",
  "Protected" = "#1B4332",
  "Priority"  = "#40916C",
  "Red"       = "#C1121F",
  "Amber"     = "#FFB000",
  "Green"     = "#74C69D",
  "Other"     = "grey80"
)

# ---- BIRD HEATMAP ----
fig_birds_heatmap <- ggplot(
  birds_classified,
  aes(
    x   = hillside,
    y   = factor(species),
    fill = classification
  )
) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(values = pal_all) +
  labs(
    title = "Bird Species by Hillside (Important)",
    x = "Hillside",
    y = "Bird Species"
  ) +
  theme_minimal()

# ---- MAMMAL HEATMAP (BATS INCLUDED) ----
fig_mammals_heatmap <- ggplot(
  mammals_classified,
  aes(
    x   = hillside,
    y   = factor(species),
    fill = classification
  )
) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(values = pal_all) +
  labs(
    title = "Protected & Priority Mammals (Including Bats)",
    x = "Hillside",
    y = "Mammal Species"
  ) +
  theme_minimal()


fig_birds_heatmap
fig_mammals_heatmap
#end of script