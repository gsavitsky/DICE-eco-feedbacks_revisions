# Meeting / paper figure: SCC (2020 pulse, welfare to 2200) decomposed by carbon-cycle feedback (Shapley values),
# Lenton (7 processes) vs DFAIR (2 aggregate terms + permafrost), same axis.
# Inputs: shapley_values.csv (ch2_scc_shapley.py), dfair_shapley_values.csv (ch2_dfair_shapley.py)
suppressPackageStartupMessages({library(readr); library(dplyr); library(ggplot2); library(patchwork)})
d <- "vis_corrected/feedback_factorial"
lab <- c(ocean = "Ocean buffer\nsaturation", oceanT = "Ocean\nwarming", fert = "CO2\nfertilisation",
         photoT = "Photosynthesis\nwarming response", vegT = "Plant\nrespiration", soilT = "Soil\nrespiration",
         perma = "Permafrost", irC = "Uptake\nsaturation (r_C)", irT = "Temperature\nterm (r_T)")
wf <- function(file, metric, title, off_lab, on_lab, ylim) {
  S <- read_csv(file.path(d, file), show_col_types = FALSE) %>% filter(metric == !!metric)
  off <- S$all_off[1]; on <- S$all_on[1]
  st <- bind_rows(tibble(lab = off_lab, value = off, type = "total"),
                  S %>% transmute(lab = lab[process], value = shapley, type = ifelse(shapley > 0, "reinforcing", "balancing")),
                  tibble(lab = on_lab, value = on, type = "total")) %>%
    mutate(i = row_number(), end = ifelse(type == "total", value, off + cumsum(ifelse(type == "total", 0, value))),
           start = ifelse(type == "total", 0, end - value))
  ggplot(st) +
    geom_rect(data = tibble(type = c("reinforcing", "balancing", "total")), aes(xmin = 1, xmax = 1, ymin = 0, ymax = 0, fill = type), alpha = 0) +
    geom_rect(aes(xmin = i - 0.4, xmax = i + 0.4, ymin = start, ymax = end, fill = type)) +
    geom_segment(data = st %>% filter(i < max(i)), aes(x = i + 0.4, xend = i + 0.6, y = end, yend = end), color = "grey40") +
    geom_text(aes(x = i, y = pmax(start, end), label = ifelse(type == "total", sprintf("$%.0f", value), sprintf("%+.1f", value))),
              vjust = -0.4, size = 4, fontface = "bold") +
    scale_x_continuous(breaks = st$i, labels = st$lab) +
    scale_fill_manual(values = c(total = "grey55", reinforcing = "#d62728", balancing = "#1f77b4"),
                      breaks = c("reinforcing", "balancing", "total"), limits = c("reinforcing", "balancing", "total"), drop = FALSE,
                      labels = c("Reinforcing (raises SCC)", "Balancing (lowers SCC)", "Total"), name = NULL) +
    coord_cartesian(ylim = ylim) +
    labs(x = NULL, y = "SCC 2020 (2019$ per tCO2)", title = title) +
    theme_bw(base_size = 13) + theme(axis.text.x = element_text(size = 9.5), panel.grid.minor = element_blank(),
                                     panel.grid.major.x = element_blank(), plot.title = element_text(face = "bold"))
}
yl <- c(0, 63)
a <- wf("shapley_values.csv", "SCC2020", "(a) Lenton carbon cycle (DICE-CP): seven processes",
        "All feedbacks at\npre-industrial", "All on\n(DICE-CP)", yl)
b <- wf("dfair_shapley_values.csv", "SCC2020", "(b) DFAIR (DICE 2023 + permafrost): two aggregate terms",
        "Feedback\nterms off", "All on (DICE 2023\n+ permafrost)", yl)
fig <- (a + theme(legend.position = "none")) /
  (b + theme(legend.position = "bottom") + guides(fill = guide_legend(override.aes = list(alpha = 1))))
ggsave(file.path(d, "fig_scc_decomposition_lenton_vs_dfair.png"), fig, width = 12, height = 10, dpi = 300)
