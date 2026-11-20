# Before vs. now: DICE-C / DICE-CP in the original submission (STELLA output, vis/) and after the corrections
# (corrected Python model, vis_corrected/). Standard damages and discounting only (runs 1, 9, 13). Window 2020-2200.
# Output: vis_corrected/figures/before_after/

suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)})

out_dir <- "vis_corrected/figures/before_after"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
runs <- c(Run_1 = "DICE 2023", Run_5 = "DICE 2023 + permafrost", Run_9 = "DICE-C", Run_13 = "DICE-CP")
vers <- c(before = "Original submission", now = "Corrected")
series_colors <- c("DICE 2023" = "#1f77b4", "DICE 2023 + permafrost" = "#9467bd", "DICE-C" = "#ff7f0e", "DICE-CP" = "#d62728")
END <- 2200

read_var <- function(dir, file, version) {
  read_csv(file.path(dir, "by_variable", file), show_col_types = FALSE) %>%
    select(Years, all_of(names(runs))) %>%
    pivot_longer(-Years, names_to = "run", values_to = "value") %>%
    mutate(config = factor(runs[run], levels = runs), version = factor(vers[version], levels = vers))
}
both <- function(file) bind_rows(read_var("vis/Welfare sensitivity", file, "before") %>% filter(config != "DICE 2023 + permafrost"),
                                 read_var("vis_corrected/Welfare sensitivity", file, "now"))

theme_ba <- theme_bw(base_size = 15) +
  theme(text = element_text(color = "black"), axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 15), strip.background = element_rect(fill = "grey92"),
        legend.position = "bottom", panel.grid.minor = element_blank())

# ---- 1. Temperature ----
tatm <- both("Climate_Atmospheric_temperature_anomaly_TATM_.csv") %>% filter(Years <= END)
p1 <- ggplot(tatm, aes(Years, value, color = config)) +
  geom_line(linewidth = 1.1) + facet_wrap(~version) +
  scale_color_manual(values = series_colors, name = NULL) +
  labs(x = NULL, y = "Warming above pre-industrial (°C)") + theme_ba
ggsave(file.path(out_dir, "ba1_temperature.png"), p1, width = 12, height = 5.5, dpi = 300)

# ---- 2. Atmospheric CO2 ----
co2 <- both("Climate_atmospheric_carbon.csv") %>% filter(Years <= END) %>% mutate(value = value / 2.124)
p2 <- ggplot(co2, aes(Years, value, color = config)) +
  geom_line(linewidth = 1.1) + facet_wrap(~version) +
  scale_color_manual(values = series_colors, name = NULL) +
  labs(x = NULL, y = expression(bold("Atmospheric CO"[2] * " (ppm)"))) + theme_ba
ggsave(file.path(out_dir, "ba2_co2.png"), p2, width = 12, height = 5.5, dpi = 300)

# ---- 3. Permafrost carbon: cumulative release to the atmosphere and frozen fraction ----
# Pools from export_permafrost_before_after.py (original: Python replication of the STELLA module, which did not save
# its permafrost pools; matches the STELLA DICE-CP temperature to within 0.01 C).
pf <- read_csv("vis_corrected/permafrost_before_after.csv", show_col_types = FALSE) %>%
  mutate(series = paste0(config, " (", version, ")"))
pf_cols <- c("DICE-CP (Original submission)" = "grey55", "DICE-CP (Corrected)" = "#d62728",
             "DICE 2023 + permafrost (Corrected)" = "#9467bd")
pf_lty <- c("DICE-CP (Original submission)" = "dashed", "DICE-CP (Corrected)" = "solid",
            "DICE 2023 + permafrost (Corrected)" = "solid")
pf_long <- pf %>% transmute(year, series, `Permafrost carbon released to the atmosphere since 2020 (PgC)` = released,
                            `Share of pre-industrial permafrost carbon still frozen` = frozen_frac) %>%
  pivot_longer(-c(year, series), names_to = "panel", values_to = "value")
p3 <- ggplot(pf_long, aes(year, value, color = series, linetype = series)) +
  geom_line(linewidth = 1.1) + facet_wrap(~panel, scales = "free_y") +
  scale_color_manual(values = pf_cols, name = NULL) + scale_linetype_manual(values = pf_lty, name = NULL) +
  labs(x = NULL, y = NULL) + theme_ba + theme(strip.text = element_text(size = 13))
ggsave(file.path(out_dir, "ba3_permafrost_carbon.png"), p3, width = 13, height = 5.8, dpi = 300)

# ---- 4. Welfare relative to DICE 2023 (accumulated utility to 2200 in both versions) ----
uacc <- both("Economy_Accumulated_Utility_Uacc_.csv") %>% filter(Years <= END) %>% group_by(version, config) %>%
  filter(Years == max(Years)) %>% ungroup()
horizon <- uacc %>% group_by(version) %>% summarise(end = max(Years))
wel <- uacc %>% group_by(version) %>% mutate(pct = 100 * (value / value[config == "DICE 2023"] - 1)) %>%
  ungroup() %>% filter(config != "DICE 2023") %>% left_join(horizon, by = "version") %>%
  mutate(version_lab = factor(paste0(version, "\n(to ", end, ")"),
                              levels = paste0(vers, "\n(to ", horizon$end[match(vers, horizon$version)], ")")))
p4 <- ggplot(wel, aes(config, pct, fill = version_lab)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%+.2f%%", pct), vjust = ifelse(pct < 0, 1.4, -0.5)),
            position = position_dodge(width = 0.75), size = 5, fontface = "bold") +
  scale_fill_manual(values = c("grey70", "grey20"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  labs(x = NULL, y = "Welfare change vs. DICE 2023 (%)") + theme_ba
ggsave(file.path(out_dir, "ba4_welfare.png"), p4, width = 8, height = 6, dpi = 300)

# ---- summary numbers for the meeting ----
pk <- tatm %>% group_by(version, config) %>% summarise(peak_T = max(value), peak_year = Years[which.max(value)],
                                                     T2100 = value[Years == 2100], T2150 = value[Years == 2150],
                                                     T2200 = value[Years == 2200], .groups = "drop")
write_csv(left_join(pk, wel %>% select(version, config, welfare_pct = pct), by = c("version", "config")),
          file.path(out_dir, "before_after_summary.csv"))
print(pk %>% mutate(across(where(is.numeric), ~ round(.x, 2))))
print(wel %>% select(version, config, end, pct))
