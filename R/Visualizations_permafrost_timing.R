# Permafrost warming over time (2020-2200): total permafrost warming and the part from carbon already thawed by 2020
# ("committed": same run with no new thaw after 2020), for three permafrost strengths, Lenton (DICE-CP) and DFAIR.
# Input: vis_corrected/ch2_permafrost_warming_timeseries.csv (ch2_warming_time.py). Output: fig_permafrost_timing.png
suppressPackageStartupMessages({library(readr); library(dplyr); library(ggplot2)})
d <- read_csv("vis_corrected/ch2_permafrost_warming_timeseries.csv", show_col_types = FALSE) %>%
  filter(m %in% c(0.5, 1, 2)) %>%
  mutate(strength = factor(case_when(m == 0.5 ~ "Half Hector (~22 PgC/°C)", m == 1 ~ "Hector (~38 PgC/°C)",
                                     TRUE ~ "Double Hector (~60 PgC/°C)"),
                           levels = c("Half Hector (~22 PgC/°C)", "Hector (~38 PgC/°C)", "Double Hector (~60 PgC/°C)")),
         carbon = factor(ifelse(carbon == "Lenton", "Lenton carbon cycle (DICE-CP)", "DFAIR (DICE 2023 + permafrost)"),
                         levels = c("Lenton carbon cycle (DICE-CP)", "DFAIR (DICE 2023 + permafrost)")))
p <- ggplot(d, aes(year)) +
  geom_area(aes(y = pmin(dT_committed, dT_total), fill = "Carbon thawed before 2020"), alpha = 0.85) +
  geom_ribbon(aes(ymin = pmin(dT_committed, dT_total), ymax = dT_total, fill = "New thaw after 2020"), alpha = 0.85) +
  geom_line(aes(y = dT_total), linewidth = 0.7) +
  facet_grid(carbon ~ strength) +
  scale_fill_manual(values = c("Carbon thawed before 2020" = "#8c564b", "New thaw after 2020" = "#e7a37a"), name = NULL) +
  scale_x_continuous(breaks = c(2050, 2100, 2150, 2200)) +
  labs(x = NULL, y = "Warming added by permafrost (°C)") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
                                   strip.text = element_text(face = "bold", size = 11))
ggsave("vis_corrected/figures/fig_permafrost_timing.png", p, width = 12, height = 7, dpi = 300)
