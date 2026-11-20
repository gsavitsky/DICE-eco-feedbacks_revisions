# Waterfall of the SCC by carbon-cycle process (Shapley values; Lenton variant iv; welfare to 2200).
# Input: vis_corrected/feedback_factorial/shapley_values.csv (ch2_scc_shapley.py). Output: fig_scc_waterfall.png
suppressPackageStartupMessages({library(readr); library(dplyr); library(ggplot2); library(patchwork)})
args <- commandArgs(trailingOnly = TRUE)
INFILE <- if (length(args) >= 1) args[1] else "shapley_values.csv"
YRS <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else c(2020, 2100)
OUTFILE <- if (length(args) >= 3) args[3] else "fig_scc_waterfall.png"
S <- read_csv(file.path("vis_corrected/feedback_factorial", INFILE), show_col_types = FALSE)
labs_p <- c(ocean = "Ocean buffer\nsaturation", oceanT = "Ocean\nwarming", fert = "CO2\nfertilisation",
            photoT = "Photosynthesis\nwarming response", vegT = "Plant\nrespiration", soilT = "Soil\nrespiration",
            perma = "Permafrost")
wf <- function(metric, title) {
  d <- S %>% filter(metric == !!metric)
  off <- d$all_off[1]; on <- d$all_on[1]
  p <- d %>% mutate(lab = labs_p[process]) %>% select(lab, value = shapley)
  steps <- bind_rows(tibble(lab = "All feedbacks\nat pre-industrial", value = off, type = "total"),
                     p %>% mutate(type = ifelse(value > 0, "reinforcing", "balancing")),
                     tibble(lab = "All feedbacks\non (DICE-CP)", value = on, type = "total")) %>%
    mutate(i = row_number(), end = ifelse(type == "total", value, off + cumsum(ifelse(type == "total", 0, value))),
           start = ifelse(type == "total", 0, end - value), lab = factor(lab, levels = lab))
  ggplot(steps) +
    geom_rect(aes(xmin = i - 0.4, xmax = i + 0.4, ymin = start, ymax = end, fill = type)) +
    geom_segment(data = steps %>% filter(i < max(i)), aes(x = i + 0.4, xend = i + 0.6, y = end, yend = end), color = "grey40") +
    geom_text(aes(x = i, y = pmax(start, end), label = ifelse(type == "total", sprintf("$%.0f", value), sprintf("%+.1f", value))),
              vjust = -0.4, size = 4.2, fontface = "bold") +
    scale_x_continuous(breaks = steps$i, labels = steps$lab) +
    scale_fill_manual(values = c(total = "grey45", reinforcing = "#d62728", balancing = "#1f77b4"), name = NULL) +
    coord_cartesian(ylim = c(min(steps$start, steps$end) * 0.8, max(steps$start, steps$end) * 1.08)) +
    labs(x = NULL, y = "SCC (2019$ per tCO2)", title = title) +
    theme_bw(base_size = 14) + theme(axis.text.x = element_text(size = 10.5), legend.position = "bottom",
                                     panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
                                     plot.title = element_text(face = "bold"))
}
fig <- wf(paste0("SCC", YRS[1]), paste0("(a) Emission in ", YRS[1])) / wf(paste0("SCC", YRS[2]), paste0("(b) Emission in ", YRS[2])) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path("vis_corrected/feedback_factorial", OUTFILE), fig, width = 12, height = 10, dpi = 300)
