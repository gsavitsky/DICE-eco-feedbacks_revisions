# Figures for the two Chapter 2 hooks (corrected models, 2026-11-16)
#   fig_scc_permafrost_bounds.png : bounding the SCC of permafrost thaw across the literature range
#   fig_permafrost_welfare_range.png : welfare and discounting across plausible permafrost outcomes
# Input: vis_corrected/ch2_permafrost_sweep.csv and ch2_runs_summary.csv
#        (from Optimizing_DICE-CP/permafrost_revision_test/ch2_final_runs.py)

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

dir.create("vis_corrected/figures", showWarnings = FALSE, recursive = TRUE)

base_theme <- theme_bw(base_size = 17) +
  theme(
    text = element_text(color = "black"),
    axis.title = element_text(color = "black", size = 20, face = "bold"),
    axis.text = element_text(color = "black", size = 19, face = "bold"),
    legend.text = element_text(color = "black", size = 19),
    legend.title = element_text(color = "black", size = 20),
    plot.tag = element_text(color = "black", size = 21, face = "bold"),
    strip.text = element_text(color = "black", size = 19, face = "bold"),
    strip.background = element_rect(fill = "grey92")
  )

sweep <- read_csv("vis_corrected/ch2_permafrost_sweep.csv", show_col_types = FALSE) %>%
  filter(m <= 8) %>%   # above m = 8 the spin-up depletes the 2020 thawed pool, so release no longer increases with m
  mutate(
    damage = factor(damage, levels = c("standard", "Howard-Sterner"),
                    labels = c("Standard damage", "Howard-Sterner damage")),
    discounting = factor(discounting, levels = c("standard", "zero time preference"),
                         labels = c("Standard discounting", "Zero time preference")),
    carbon = factor(carbon, levels = c("Lenton", "DFAIR"),
                    labels = c("Lenton carbon cycle (DICE-CP)", "DFAIR (DICE 2023 + permafrost)")),
    cell = interaction(damage, discounting, sep = ", ")
  )

# literature range of the permafrost feedback (PgC released per degree C by 2100)
lit <- tibble(
  xmin = 3.1, xmax = 41, central = 18, high = 63
)
lit_layers_y <- function(ymin = -Inf, ymax = Inf) list(
  annotate("rect", xmin = lit$xmin, xmax = lit$xmax, ymin = ymin, ymax = ymax, fill = "grey80", alpha = 0.45),
  geom_vline(xintercept = lit$central, linetype = "dashed", linewidth = 0.6),
  geom_vline(xintercept = lit$high, linetype = "dotted", linewidth = 0.8)
)
lit_layers <- list(
  annotate("rect", xmin = lit$xmin, xmax = lit$xmax, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.45),
  geom_vline(xintercept = lit$central, linetype = "dashed", linewidth = 0.6),
  geom_vline(xintercept = lit$high, linetype = "dotted", linewidth = 0.8)
)
lit_note <- "Shading: AR6 range (3.1-41 PgC/°C); dashed: AR6 central (18); dotted: Schädel et al. 2026 (~63)"

cell_colors <- c(
  "Standard damage, Standard discounting" = "#1f77b4",
  "Howard-Sterner damage, Standard discounting" = "#d62728",
  "Standard damage, Zero time preference" = "#2ca02c",
  "Howard-Sterner damage, Zero time preference" = "#9467bd"
)

# ---- Hook 1: SCC bounds ----
p_pct <- ggplot(sweep, aes(PgC_per_C_2100, SCC_pct, color = cell, linetype = carbon)) +
  lit_layers +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_color_manual(values = cell_colors, name = NULL, guide = guide_legend(ncol = 2)) +
  scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  labs(x = "Permafrost feedback (PgC per °C, 2100)", y = "Change in 2020 SCC (%)", tag = "a") +
  base_theme

p_abs <- ggplot(sweep %>% filter(carbon == levels(carbon)[1]),
                aes(PgC_per_C_2100, SCC, color = cell)) +
  lit_layers_y(30, 15000) +
  geom_hline(aes(yintercept = SCC_noperma, color = cell), linewidth = 0.4, linetype = "longdash", show.legend = FALSE) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_color_manual(values = cell_colors, name = NULL, guide = guide_legend(ncol = 2)) +
  scale_y_log10() +
  labs(x = "Permafrost feedback (PgC per °C, 2100)", y = "2020 SCC (2019$ per tCO₂, log scale)", tag = "b") +
  base_theme

fig1 <- (p_pct | p_abs) + plot_layout(guides = "collect") +
  plot_annotation(caption = lit_note) & theme(legend.position = "bottom", legend.box = "vertical",
                                               plot.caption = element_text(size = 15))
ggsave("vis_corrected/figures/fig_scc_permafrost_bounds.png", fig1, width = 20, height = 9, dpi = 300)

# ---- Hook 2: welfare and discounting across the range of permafrost outcomes ----
p_w <- ggplot(sweep, aes(PgC_per_C_2100, welfare_pct, color = cell, linetype = carbon)) +
  lit_layers +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_color_manual(values = cell_colors, name = NULL, guide = guide_legend(ncol = 2)) +
  scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  labs(x = "Permafrost feedback (PgC per °C, 2100)", y = "Change in welfare, U_acc (%)", tag = "a") +
  base_theme

p_t <- ggplot(sweep %>% filter(damage == "Standard damage", discounting == "Standard discounting"),
              aes(PgC_per_C_2100, dT_peak, linetype = carbon)) +
  lit_layers +
  geom_line(linewidth = 1, color = "black") + geom_point(size = 2, color = "black") +
  scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  labs(x = "Permafrost feedback (PgC per °C, 2100)", y = "Added peak warming (°C)", tag = "b") +
  base_theme

fig2 <- (p_w | p_t) + plot_layout(widths = c(1.4, 1), guides = "collect") +
  plot_annotation(caption = lit_note) & theme(legend.position = "bottom", legend.box = "vertical",
                                               plot.caption = element_text(size = 15))
ggsave("vis_corrected/figures/fig_permafrost_welfare_range.png", fig2, width = 20, height = 9, dpi = 300)

cat("Hook figures written to vis_corrected/figures\n")

# ---- Emissions under the DICE 2023 cost-benefit optimal controls (Reviewer 1) ----
em_file <- "vis_corrected/ch2_emissions_std.csv"
if (file.exists(em_file)) {
  em <- read_csv(em_file, show_col_types = FALSE) %>%
    filter(Years <= 2200) %>%
    pivot_longer(-Years, names_to = "Config", values_to = "E") %>%
    mutate(Config = factor(Config, levels = c("DICE 2023", "DFAIR + permafrost", "DICE-C", "DICE-CP")))
  p_em <- ggplot(em, aes(Years, E, color = Config, linetype = Config)) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("DICE 2023" = "#1f77b4", "DFAIR + permafrost" = "#9467bd",
                                  "DICE-C" = "#ff7f0e", "DICE-CP" = "#d62728"), name = NULL) +
    scale_linetype_manual(values = c("solid", "longdash", "dashed", "dotdash"), name = NULL) +
    labs(x = "Year", y = "Total CO₂ emissions (GtCO₂ per year)") +
    base_theme + theme(legend.position = c(0.75, 0.8), legend.background = element_blank())
  ggsave("vis_corrected/figures/fig_emissions.png", p_em, width = 11, height = 7, dpi = 300)
}
