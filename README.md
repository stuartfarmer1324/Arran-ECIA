# Arran-ECIA

This repository contains the R script used to generate the Arran biodiversity baseline plots (`Baseline graphs for Arran.R`). Any edits you make to that file directly change the analysis workflow: when you save the file and re-run it in R, the updated logic and annotations are executed. There is no intermediate build step—Git tracks the plain R code itself—so commits reflect exactly what will run.

## Running the script
1. Open `Baseline graphs for Arran.R` in your R session or IDE.
2. Ensure the required packages are installed (`tidyverse`, `stringr`, `fuzzyjoin`, `stringdist`, `forcats`, `scales`). The `patchwork` package is optional and only used if you want multi-panel layouts.
3. Set the `OCCURS_PATH` environment variable to point at your `associatedoccurences.csv` file. If the file is missing, the script automatically uses the built-in demo tibble so you can test the pipeline.
4. Source the script (e.g., `source("Baseline graphs for Arran.R")`) to produce the plots and summaries.

## Editing guidance
- Comments in the script describe why each processing step exists (e.g., genus imputation, grouping rules). Removing or changing those lines changes the underlying logic; keep them if you want the documented behavior.
- Because the script is linear, edits apply in order—changing an early transformation (such as field standardisation or grouping) will affect all downstream plots and tables.
- When you commit changes, Git records the exact lines added or removed. Anyone pulling the repository will run the updated code without extra configuration.
