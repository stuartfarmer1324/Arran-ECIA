library(ggplot2)
library(dplyr)
library(gridExtra)

# -----------------------------
# DATA
# -----------------------------
north <- data.frame(
  class = c("A2.1","A3.1","B1.2","C1.1","D1.1"),
  percent = c(1.3, 0.8, 57.8, 12.8, 27.2)
)

south <- data.frame(
  class = c("A2.2","A3.1","B1.2","B4","B5","C2","C1.2","D5","E1.7","J1.4","A1.1.1"),
  percent = c(7.7,0.8,8.0,0.8,2.9,6.4,56.9,1.8,2.4,1.4,11.0)
)

# -----------------------------
# COLOURS
# -----------------------------
palette <- c(
  "A1.1.1"="#24591C","A2.1"="#4E8C31","A2.2"="#7CCF60","A3.1"="#2F6B1A",
  "B1.2"="#E8D163","B4"="#E69F00","B5"="#A3C2FF",
  "C1.1"="#7A5C2E","C1.2"="#C8A96A","C2"="#A2A743",
  "D1.1"="#C77C6D","D5"="#A95E55",
  "E1.7"="#67A9CF","J1.4"="#984EA3"
)

# -----------------------------
# PIE FUNCTION (NO LABELS)
# -----------------------------
make_pie <- function(df, title){
  
  df <- df %>%
    mutate(
      fraction = percent / sum(percent),
      legend_label = paste0(class, " – ", percent, "%")
    )
  
  ggplot(df, aes(x="", y=fraction, fill=class)) +
    geom_col(color="white") +
    coord_polar(theta="y") +
    scale_fill_manual(
      values = palette,
      breaks = df$class,
      labels = df$legend_label
    ) +
    theme_void() +
    ggtitle(title) +
    theme(
      plot.title = element_text(size=22, face="bold", hjust=0.5),
      legend.title = element_blank(),
      legend.text = element_text(size=12)
    )
}

# -----------------------------
# Generate each pie
# -----------------------------
p_north <- make_pie(north, "North slope")
p_south <- make_pie(south, "South slope")

# -----------------------------
# Final layout:
# Pie A: Legend LEFT   | Pie
# Pie B: Pie           | Legend RIGHT
# -----------------------------
grid.arrange(
  # Panel A — legend left, pie right
  arrangeGrob(p_north, ncol=1),
  
  # Panel B — pie left, legend right
  arrangeGrob(p_south, ncol=1),
  
  ncol = 2,
  widths = c(1,1)
)

