# Chapter 3 figures, corrected physics (generated from ch4_figures.R, 2026-11-17).
# Reads figdata_corrected/ (ch3_outcomes_corrected.py); writes figures_ch3_corrected/. Reporting window 2020-2300.
# =============================================================================
# Chapter 4 (re-optimization paper) figures
#
# Reads ONLY figdata/*.csv, written by chapter4_outcomes_table.py from the final
# saved trajectories (DICE 2023 dt=1; DICE-C and DICE-CP at 5 substeps/yr).
# Every annotated number is computed from those files, not typed in, so the
# figures always agree with chapter4_outcomes_table.md.
#
# Regenerate the data first if any trajectory changes:
#   ~/anaconda3/envs/dice/bin/python chapter4_outcomes_table.py
# Then, from this folder:
#   Rscript ch4_figures.R
# Output: figures_ch4/*.pdf, *.png and caption_notes.md (caption numbers, computed from the data)
# Figure order: 1 policy, 2 warming, 3 permafrost, 4 atmospheric carbon, 5 welfare differences,
# 6 DICE-CP temperature under GAMS vs re-optimized policy (or 5+6 combined), A1 resolution, A2 multiple starts.
#
# Style follows Chapter 3 (DICE_eco_feedbacks/DICE2023 Rproj/Visualizations.R):
# same configuration colours/linetypes and the same base theme.
# =============================================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

# Run from the script's own folder (works with Rscript and RStudio "Source")
this_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NULL)
if (is.null(this_dir)) {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) this_dir <- dirname(normalizePath(f))
}
if (!is.null(this_dir)) setwd(this_dir)

out_dir <- Sys.getenv("CH4_FIG_DIR", "figures_ch3_corrected")
dir.create(out_dir, showWarnings = FALSE)

# Set TRUE to put atmospheric carbon in Figure 2 (warming) as panel (c) instead of a separate Figure 4
ATM_CARBON_AS_FIG2D <- FALSE
# Industrial emissions (Fig 1c) in GtCO2/yr, as DICE reports them; TRUE plots GtC/yr instead
EMISSIONS_IN_GTC <- FALSE
# Appendix figures A1 (resolution) and A2 (multiple starts)
MAKE_APPENDIX <- TRUE
# TRUE: welfare differences and the DICE-CP temperature comparison as one two-panel Figure 5
COMBINE_FIG5_6 <- FALSE
# End year of the policy figure (Fig 1). All other figures run 2020-2300.
FIG1_END <- as.numeric(Sys.getenv("CH4_FIG1_END", "2150"))

# Standard x-axis ticks: every 100 years on the full 2020-2300 horizon,
# every 20 years (40 in the small panels) on a short window.
x_breaks <- function(end = 2300, small = FALSE)
  if (end > 2200) c(2020, seq(2050, 2300, 50)) else seq(2020, end, if (small) 40 else 20)
x_full <- scale_x_continuous(limits = c(2020, 2300), breaks = x_breaks())

# ---- data ----------------------------------------------------------------------
cfg_levels <- c("DICE 2023", "DICE-C", "DICE-CP")
ts   <- read_csv("figdata_corrected/ch4_timeseries.csv", show_col_types = FALSE) %>%
  mutate(config = factor(config, levels = cfg_levels))
ch34 <- read_csv("figdata_corrected/ch4_ch3_vs_ch4_dicecp.csv", show_col_types = FALSE)
runs <- read_csv("figdata_corrected/ch4_runs_resolution_multistart.csv", show_col_types = FALSE)
uopt <- read_csv("figdata_corrected/ch4_uacc_opt_vs_gams.csv", show_col_types = FALSE)
topt <- read_csv("figdata_corrected/ch4_tatm_opt_vs_gams.csv", show_col_types = FALSE)

GTCO2_PER_GTC <- 3.666   # same factor the models use

# ---- style (Chapter 3) ----------------------------------------------------------
base_theme <- theme_bw(base_size = 17) +
  theme(
    text = element_text(color = "black"),
    axis.title = element_text(color = "black", size = 20, face = "bold"),
    axis.text = element_text(color = "black", size = 19, face = "bold"),
    legend.text = element_text(color = "black", size = 19),
    legend.title = element_text(color = "black", size = 20),
    plot.tag = element_text(color = "black", size = 21, face = "bold"),
    plot.margin = margin(t = 20, r = 12, b = 10, l = 10)
  )
series_colors    <- c("DICE 2023" = "#1f77b4", "DICE-C" = "#ff7f0e", "DICE-CP" = "#d62728")
series_linetypes <- c("DICE 2023" = "solid",   "DICE-C" = "dashed",  "DICE-CP" = "dotdash")
lw <- 1.3              # line width for model series
ann_size <- 5.5        # annotation text size (mm)

scale_cfg <- list(
  scale_color_manual(values = series_colors, name = NULL),
  scale_linetype_manual(values = series_linetypes, name = NULL)
)

# PDF device that can draw μ, °, → and en dashes: Cairo if it loads (needs XQuartz on
# macOS), otherwise the native macOS Quartz PDF device, otherwise plain pdf().
pdf_device <- local({
  ok <- function(dev) tryCatch({ f <- tempfile(fileext = ".pdf"); dev(f); plot.new(); text(0.5, 0.5, "μ °C →"); dev.off(); TRUE },
                               error = function(e) FALSE, warning = function(w) FALSE)
  if (ok(function(f) cairo_pdf(f))) cairo_pdf
  else if (.Platform$OS.type == "unix" && capabilities("aqua") &&
           ok(function(f) quartz(type = "pdf", file = f)))
    function(filename, width, height, ...) quartz(type = "pdf", file = filename, width = width, height = height)
  else pdf
})

save_fig <- function(p, name, width, height) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = height, device = pdf_device)
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  message("saved ", name)
}

# ---- derived values (all from the CSVs) -----------------------------------------
first_year <- function(df, cond) { y <- df$year[cond]; if (length(y)) min(y) else NA }

mu1 <- ts %>% filter(!is.na(miu)) %>% group_by(config) %>%
  summarise(year = min(year[miu >= 1 - 1e-6]), .groups = "drop")   # mu first reaches 1.0

cum_ind <- ts %>% filter(year >= 2020, year <= 2099) %>% group_by(config) %>%
  summarise(cum_gtc = sum(eind_gtco2) / GTCO2_PER_GTC, .groups = "drop")         # 2020-2100, GtC

cp  <- ts %>% filter(config == "DICE-CP")
c0  <- cp$Cperm[cp$year == 2020]
loss_total <- c0 - cp$Cperm[cp$year == 2300]
y50 <- first_year(cp, c0 - cp$Cperm >= 0.5 * loss_total)
y90 <- first_year(cp, c0 - cp$Cperm >= 0.9 * loss_total)
thaw_pk <- cp %>% slice_max(thawedC, n = 1, with_ties = FALSE)

loss_2100   <- c0 - cp$Cperm[cp$year == 2100]
offset_2100 <- cum_ind$cum_gtc[cum_ind$config == "DICE-C"] - cum_ind$cum_gtc[cum_ind$config == "DICE-CP"]

ca_pk <- ts %>% group_by(config) %>% slice_max(Ca, n = 1, with_ties = FALSE) %>% ungroup()

wide_T <- ts %>% select(config, year, TATM) %>% pivot_wider(names_from = config, values_from = TATM)
yrs_above <- ts %>% filter(year <= 2419) %>% group_by(config) %>%
  summarise(`1.5 °C` = sum(TATM > 1.5), `2 °C` = sum(TATM > 2), `2.5 °C` = sum(TATM > 2.5), .groups = "drop") %>%
  pivot_longer(-config, names_to = "threshold", values_to = "years") %>%
  mutate(threshold = factor(threshold, levels = c("1.5 °C", "2 °C", "2.5 °C")))


# Caption notes: numbers for the captions are collected here and written to
# figures_ch4/caption_notes.md at the end, so captions never need retyping.
notes <- character()
note <- function(...) notes <<- c(notes, paste0(...))
legend_top <- theme(legend.position = "top", legend.key.width = unit(2.5, "cm"))
cp_col <- series_colors[["DICE-CP"]]

# =============================================================================
# Figure 1. Optimal policy: a planner facing permafrost abates sooner
# =============================================================================
pol <- ts %>% filter(year <= FIG1_END, !is.na(miu))

p1a <- ggplot(pol, aes(year, miu)) +
  geom_line(aes(colour = config, linetype = config), linewidth = lw) +
  scale_cfg + scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1.2, 0.2)) +
  scale_x_continuous(limits = c(2020, FIG1_END), breaks = x_breaks(FIG1_END)) +
  labs(x = "Year", y = "Emissions control rate (μ)", tag = "(a)") +
  base_theme + legend_top

p1b <- ggplot(ts %>% filter(year <= FIG1_END, !is.na(s)), aes(year, s)) +
  geom_line(aes(colour = config, linetype = config), linewidth = lw) +
  scale_cfg + scale_x_continuous(limits = c(2020, FIG1_END), breaks = x_breaks(FIG1_END, small = TRUE)) +
  labs(x = "Year", y = "Savings rate (s)", tag = "(b)") +
  base_theme + theme(legend.position = "none")

# Emissions units: DICE reports emissions in GtCO2/yr; set EMISSIONS_IN_GTC <- TRUE (top of
# script) to plot GtC/yr instead, matching the GtC units of the carbon stocks in Figs 3-4.
emis <- ts %>% filter(year <= FIG1_END, !is.na(eind_gtco2)) %>%
  mutate(e = if (EMISSIONS_IN_GTC) eind_gtco2 / GTCO2_PER_GTC else eind_gtco2)
emis_ylab <- if (EMISSIONS_IN_GTC) expression(atop(bold("Industrial CO"[2] * " emissions"), bold("(GtC yr"^-1 * ")"))) else
  expression(atop(bold("Industrial CO"[2] * " emissions"), bold("(GtCO"[2] * " yr"^-1 * ")")))
p1c <- ggplot(emis, aes(year, e)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.5) +
  geom_line(aes(colour = config, linetype = config), linewidth = lw) +
  scale_cfg + scale_x_continuous(limits = c(2020, FIG1_END), breaks = x_breaks(FIG1_END, small = TRUE)) +
  labs(x = "Year", y = emis_ylab, tag = "(c)") +
  base_theme + theme(legend.position = "none")

fig1 <- p1a / (p1b | p1c) + plot_layout(heights = c(1.3, 1))
save_fig(fig1, "fig1_optimal_policy", 16, 14)

note("## Figure 1 (optimal policy)\n",
     "- (a) μ first reaches 1.0 in ", paste0(mu1$year, " (", mu1$config, ")", collapse = ", "),
     ". The GAMS μ ceiling binds in 2030–34 (0.24) for all three, and in 2035–39 (0.36) for DICE-CP only; ",
     "it rises to 1.1 in 2120, allowing net-negative industrial CO₂.\n",
     "- (b) The short savings spike around 2120 coincides with the μ ceiling rising from 1.0 to 1.1 ",
     "(inherited GAMS schedule; anticipation is observed, not derived). From 2205",
     if (FIG1_END < 2205) " (outside this window)" else "", " s is fixed at 0.28, as in GAMS.\n",
     "- (c) Cumulative industrial CO₂ 2020–2100: ",
     paste0(cum_ind$config, " ", round(cum_ind$cum_gtc), " GtC", collapse = ", "), ".\n")

# =============================================================================
# Figure 2. Warming that optimal policy does not offset (+ optional (c) carbon)
# =============================================================================
p2a <- ggplot(ts, aes(year, TATM)) +
  geom_line(aes(colour = config, linetype = config), linewidth = lw) +
  scale_cfg + x_full +
  labs(x = "Year", y = "Atmospheric temperature\nrise (°C)", tag = "(a)") +
  base_theme + legend_top

p2b <- ggplot(yrs_above, aes(threshold, years, fill = config)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = years), position = position_dodge(width = 0.8), vjust = -0.4, size = ann_size - 0.5) +
  scale_fill_manual(values = series_colors, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Threshold", y = "Years above\nthreshold (2020–2419)", tag = "(b)") +
  base_theme + theme(legend.position = "none")

# atmospheric carbon: Figure 4, or Figure 2(c) if ATM_CARBON_AS_FIG2D
p_ca <- ggplot(ts, aes(year, Ca)) +
  geom_line(aes(colour = config, linetype = config), linewidth = lw) +
  scale_cfg + x_full +
  scale_y_continuous(labels = comma) +
  labs(x = "Year", y = "Atmospheric carbon (GtC)") +
  base_theme

if (ATM_CARBON_AS_FIG2D) {
  fig2 <- p2a / (p2b | (p_ca + labs(tag = "(c)") + theme(legend.position = "none"))) +
    plot_layout(heights = c(1.2, 1))
  save_fig(fig2, "fig2_warming", 16, 13)
} else {
  fig2 <- p2a / p2b + plot_layout(heights = c(1.2, 1))
  save_fig(fig2, "fig2_warming", 16, 13)
  save_fig(p_ca + legend_top, "fig4_atmospheric_carbon", 14, 7.5)
}

onfloor <- ts %>% filter(config == "DICE-C") %>% mutate(f = TATM < 0.5 + 1e-3)
floor_first <- min(onfloor$year[onfloor$f])
floor_from  <- max(onfloor$year[!onfloor$f]) + 1          # start of the final stretch on the floor
ca_early <- ts %>% filter(config %in% c("DICE-C", "DICE-CP"), year <= 2030) %>% group_by(config) %>%
  summarise(c0 = Ca[year == 2020], cmin = min(Ca), ymin = year[which.min(Ca)], .groups = "drop")
note("## Figure 2 (warming)\n",
     "- (a) DICE-C first reaches 0.5 °C in ", floor_first, " and stays there from ", floor_from,
     ": this is the lower bound on temperature inherited from the GAMS code (TATM.lo = 0.5, ",
     "imposed by Nordhaus for numerical stability), not a physical result.\n",
     "- (b) Colours as in (a); no separate legend.\n",
     if (ATM_CARBON_AS_FIG2D) "- (c) See Figure 4 notes.\n" else "")
note("## Figure 4 (atmospheric carbon)\n",
     "- Peaks: ", paste0(sprintf("%s %s GtC in %d", ca_pk$config, comma(round(ca_pk$Ca)), ca_pk$year), collapse = "; "), ".\n",
     "- DICE-C and DICE-CP start from the 2020 state of a historical spin-up of the Lenton carbon cycle (419 ppm), ",
     "so there is no initial transient; DICE 2023 starts from its own 2020 DFAIR state.\n")

# =============================================================================
# Figure 3. Permafrost carbon in DICE-CP: GAMS DICE 2023 policy vs re-optimized policy
#   "GAMS policy" = the Chapter 4 DICE-CP model run with the GAMS DICE 2023 controls
#   (the Chapter 3 policy); it reproduces Chapter 3 to within 0.01 °C.
# =============================================================================
pol_lab <- c(gams = "GAMS DICE 2023 policy", opt = "Re-optimized policy")
perm <- ch34 %>%
  select(year, Cperm_gams = Cperm_ch4_model_gams_controls, Cperm_opt = Cperm_ch4_optimized,
         thawedC_gams = thawedC_ch4_model_gams_controls, thawedC_opt = thawedC_ch4_optimized) %>%
  pivot_longer(-year, names_to = c(".value", "policy"), names_sep = "_") %>%
  mutate(policy = factor(pol_lab[policy], levels = pol_lab))
pol_scale <- scale_linetype_manual(values = setNames(c("dotdash", "solid"), pol_lab), name = NULL)

p3a <- ggplot(perm, aes(year, Cperm, linetype = policy)) +
  geom_line(colour = cp_col, linewidth = lw) + pol_scale +
  x_full +
  labs(x = "Year", y = "Frozen permafrost carbon (GtC)", tag = "(a)") + base_theme

p3b <- ggplot(perm, aes(year, thawedC, linetype = policy)) +
  geom_line(colour = cp_col, linewidth = lw) + pol_scale +
  x_full +
  labs(x = "Year", y = "Thawed permafrost carbon (GtC)", tag = "(b)") + base_theme

fig3 <- (p3a | p3b) + plot_layout(guides = "collect") &
  theme(legend.position = "top", legend.key.width = unit(2.5, "cm"), legend.box.margin = margin(b = 10))
fig3 <- fig3 + plot_annotation(theme = theme(plot.margin = margin(t = 5)))
save_fig(fig3, "fig3_permafrost_carbon", 16, 8)

pw <- perm %>% select(year, policy, Cperm, thawedC)
lossf <- function(pl, y) c0 - pw$Cperm[pw$policy == pl & pw$year == y]
tpk <- pw %>% group_by(policy) %>% slice_max(thawedC, n = 1, with_ties = FALSE) %>% ungroup()
note("## Figure 3 (permafrost carbon, DICE-CP)\n",
     "- Both lines are DICE-CP. The GAMS-policy line is the same model run with the GAMS DICE 2023 ",
     "controls (the Chapter 2 runs).\n",
     sprintf("- (a) Frozen carbon lost by 2100: %.0f GtC (GAMS policy) vs %.0f GtC (re-optimized); by 2300: %.0f vs %.0f GtC.\n",
             lossf(pol_lab[["gams"]], 2100), lossf(pol_lab[["opt"]], 2100),
             lossf(pol_lab[["gams"]], 2300), lossf(pol_lab[["opt"]], 2300)),
     "- (b) Thawed pool peaks: ", paste0(sprintf("%.1f GtC in %d (%s)", tpk$thawedC, tpk$year, tpk$policy), collapse = "; "), ".\n")

# =============================================================================
# Figure 5. Welfare: U_acc differences (raw units)
# Figure 6. DICE-CP temperature, GAMS DICE 2023 policy vs re-optimized policy
#   With COMBINE_FIG5_6 = TRUE both are written as one figure (fig5_welfare_temperature).
#   "GAMS policy" = Chapter 4 DICE-CP model run with the GAMS controls, as in Figure 3.
# =============================================================================
pct <- function(a, b) 100 * (a / b - 1)
# Delta U_acc relative to optimized DICE 2023; colour = model, line type = policy
du <- uopt %>%
  transmute(year,
            `DICE-C_opt` = dicec_opt - dice2023_opt, `DICE-C_gams` = dicec_gams - dice2023_opt,
            `DICE-CP_opt` = dicecp_opt - dice2023_opt, `DICE-CP_gams` = dicecp_gams - dice2023_opt) %>%
  pivot_longer(-year, names_to = c("config", "policy"), names_sep = "_", values_to = "d") %>%
  mutate(config = factor(config, levels = cfg_levels),
         policy = factor(pol_lab[policy], levels = pol_lab))

p5 <- ggplot(du, aes(year, d, colour = config, linetype = policy)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = lw) +
  scale_color_manual(values = series_colors, name = NULL, drop = TRUE) +
  pol_scale + x_full +
  labs(x = "Year", y = expression(bold(Delta * U[acc] * " from DICE 2023"))) +
  base_theme + legend_top + theme(legend.box = "vertical")

c34 <- ch34 %>% select(year, gams = TATM_ch4_model_gams_controls, opt = TATM_ch4_optimized)
long34 <- c34 %>% pivot_longer(-year, names_to = "policy", values_to = "TATM") %>%
  mutate(policy = factor(pol_lab[policy], levels = pol_lab))
p6 <- ggplot() +
  geom_ribbon(data = c34, aes(year, ymin = opt, ymax = gams), fill = cp_col, alpha = 0.15) +
  geom_line(data = long34, aes(year, TATM, linetype = policy), colour = cp_col, linewidth = lw) +
  pol_scale + x_full +
  labs(x = "Year", y = "DICE-CP atmospheric\ntemperature rise (°C)") +
  base_theme + legend_top

if (COMBINE_FIG5_6) {
  save_fig((p5 + labs(tag = "(a)")) / (p6 + labs(tag = "(b)")), "fig5_welfare_temperature", 14, 14)
} else {
  save_fig(p5, "fig5_welfare_cost", 14, 8)
  save_fig(p6, "fig6_ch3_vs_ch4", 14, 7.5)
}

fig5_name <- if (COMBINE_FIG5_6) "Figure 5(a)" else "Figure 5"
fig6_name <- if (COMBINE_FIG5_6) "Figure 5(b)" else "Figure 6"
u_at <- function(col, y) uopt[[col]][uopt$year == y]
note("## ", fig5_name, " (welfare differences)\n",
     "- ΔU_acc relative to the optimized DICE 2023 run (model units; no scale_2 constant). Solid = re-optimized policy, ",
     "dot-dash = GAMS DICE 2023 policy (the Chapter 2 runs: the same models driven by the GAMS controls).\n",
     sprintf("- U_acc(2300): DICE 2023 %.1f; DICE-C %.1f re-optimized / %.1f GAMS policy; DICE-CP %.1f / %.1f.\n",
             u_at("dice2023_opt", 2300), u_at("dicec_opt", 2300), u_at("dicec_gams", 2300),
             u_at("dicecp_opt", 2300), u_at("dicecp_gams", 2300)),
     sprintf("- ΔU_acc(2300) from DICE 2023: DICE-C %+.1f / %+.1f; DICE-CP %+.1f / %+.1f (re-optimized / GAMS policy).\n",
             u_at("dicec_opt", 2300) - u_at("dice2023_opt", 2300), u_at("dicec_gams", 2300) - u_at("dice2023_opt", 2300),
             u_at("dicecp_opt", 2300) - u_at("dice2023_opt", 2300), u_at("dicecp_gams", 2300) - u_at("dice2023_opt", 2300)),
     sprintf("- Gain from re-optimizing: DICE-C %+.1f (%+.3f%%), DICE-CP %+.1f (%+.3f%%).\n",
             u_at("dicec_opt", 2300) - u_at("dicec_gams", 2300), pct(u_at("dicec_opt", 2300), u_at("dicec_gams", 2300)),
             u_at("dicecp_opt", 2300) - u_at("dicecp_gams", 2300), pct(u_at("dicecp_opt", 2300), u_at("dicecp_gams", 2300))))
pkg <- c34 %>% slice_max(gams, n = 1); pko <- c34 %>% slice_max(opt, n = 1)
dmod <- max(abs(ch34$TATM_ch4_model_gams_controls - ch34$TATM_ch3))
note("## ", fig6_name, " (DICE-CP temperature, GAMS policy vs re-optimized)\n",
     sprintf("- Re-optimizing changes DICE-CP peak warming from %.2f °C (%d) to %.2f °C (%d), %+.2f °C; ",
             pkg$gams, pkg$year, pko$opt, pko$year, pko$opt - pkg$gams),
     sprintf("in 2300 from %.2f to %.2f °C (%+.2f °C).\n",
             c34$gams[c34$year == 2300], c34$opt[c34$year == 2300], c34$opt[c34$year == 2300] - c34$gams[c34$year == 2300]),
     "- Both lines use the same model; the difference is due to the policy alone.\n")


# =============================================================================
# Figure 6 (five configurations). Temperature: DICE 2023 (optimized), and DICE-C and
#   DICE-CP under the GAMS DICE 2023 policy and the re-optimized policy.
#   Same layout as Figure 5: colour = model, line type = policy, legends on top.
# =============================================================================
t5 <- topt %>%
  transmute(year, `DICE 2023_opt` = dice2023_opt, `DICE-C_gams` = dicec_gams, `DICE-C_opt` = dicec_opt,
            `DICE-CP_gams` = dicecp_gams, `DICE-CP_opt` = dicecp_opt) %>%
  pivot_longer(-year, names_to = c("config", "policy"), names_sep = "_", values_to = "TATM") %>%
  mutate(config = factor(config, levels = cfg_levels), policy = factor(pol_lab[policy], levels = pol_lab))
p6five <- ggplot(t5, aes(year, TATM, colour = config, linetype = policy)) +
  geom_line(linewidth = lw) +
  scale_color_manual(values = series_colors, name = NULL) +
  pol_scale + x_full +
  labs(x = "Year", y = "Atmospheric temperature\nrise (°C)") +
  base_theme + legend_top + theme(legend.box = "vertical")
save_fig(p6five, "fig6_temp_five_configs", 14, 8)
pk5 <- t5 %>% group_by(config, policy) %>% slice_max(TATM, n = 1, with_ties = FALSE) %>% ungroup()
note("## Figure 6, five configurations (temperature)\n",
     "- DICE 2023 is shown at its own optimum (solid). DICE-C and DICE-CP: dot-dash = GAMS DICE 2023 policy ",
     "(Chapter 3 models driven by the GAMS controls), solid = re-optimized.\n",
     "- Peaks: ", paste0(sprintf("%s, %s: %.2f °C (%d)", pk5$config, pk5$policy, pk5$TATM, pk5$year), collapse = "; "), ".\n",
     "- DICE-CP's GAMS-policy peak is flat to ±0.0001 °C over 2124–2125.\n")

if (MAKE_APPENDIX) {
# =============================================================================
# Appendix A1. Time-step resolution: 1 vs 2 vs 5 substeps per year
# =============================================================================
keep_a1 <- tribble(~config, ~substeps, ~start,                     # one run per config and resolution
  "DICE-C", 1, "GAMS-controls start",
  "DICE-C", 2, "conservative start",
  "DICE-C", 5, "warm start from the 2-substep solution",
  "DICE-CP", 1, "warm start (earlier solution)",
  "DICE-CP", 2, "warm start from the annual-grid solution",
  "DICE-CP", 5, "warm start from the 2-substep solution")
res <- runs %>% filter(year <= 2150) %>% semi_join(keep_a1, by = c("config", "substeps", "start")) %>%
  mutate(res = factor(paste0(substeps, ifelse(substeps == 1, " substep/yr (annual)", " substeps/yr")),
                      levels = c("1 substep/yr (annual)", "2 substeps/yr", "5 substeps/yr")))
pA1 <- ggplot(res, aes(year, miu, colour = res, linetype = res, linewidth = res)) +
  geom_line() +
  facet_wrap(~config) +
  scale_color_manual(values = c("grey55", "#4a4a4a", "black"), name = NULL) +
  scale_linetype_manual(values = c("solid", "solid", "dashed"), name = NULL) +
  scale_linewidth_manual(values = c(1.3, 2.2, 1.1), name = NULL) +
  scale_x_continuous(breaks = seq(2020, 2150, 40)) +
  labs(x = "Year", y = "Emissions control rate (μ)") +
  base_theme + theme(legend.position = "top", legend.key.width = unit(2.5, "cm"),
                     strip.text = element_text(size = 19, face = "bold"))
save_fig(pA1, "figA1_resolution", 16, 7.5)
note("## Figure A1 (time-step resolution)\n",
     "- The annual grid (1 substep per year) is numerically unstable for the corrected surface ocean; ",
     "solutions at 2 and 5 substeps per year coincide, and results are reported at 5.\n")

# =============================================================================
# Appendix A2. DICE-CP from multiple starting points
#   (a) reporting resolution (2 and 5 substeps/yr): all starts reach the same optimum
#   (b) annual grid: the cold start ends in a worse local optimum
# =============================================================================
ms <- runs %>% filter(config == "DICE-CP")
pal_a <- c("warm start from the 2-substep solution" = cp_col,
           "aggressive start" = "#7f3b08", "warm start from the annual-grid solution" = "#f4a3a3")
ms_a <- ms %>% filter(substeps >= 2) %>%
  mutate(lab = factor(paste0(substeps, " substeps/yr, ", start),
                      levels = c("5 substeps/yr, warm start from the 2-substep solution",
                                 "2 substeps/yr, aggressive start",
                                 "2 substeps/yr, warm start from the annual-grid solution")))
pA2a <- ggplot(ms_a, aes(year, miu, colour = lab, linetype = lab, linewidth = lab)) +
  geom_line() +
  scale_color_manual(values = unname(pal_a), name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = NULL) +
  scale_linewidth_manual(values = c(2.6, 1.2, 1.4), name = NULL) +
  x_full +
  coord_cartesian(ylim = c(0, 1.15)) +
  guides(colour = guide_legend(ncol = 1)) +
  labs(x = "Year", y = "DICE-CP emissions\ncontrol rate (μ)", tag = "(a)",
       subtitle = "Reporting resolution (2 and 5 substeps per year)") +
  base_theme + theme(legend.position = "top", legend.key.width = unit(2.2, "cm"),
                     legend.text = element_text(size = 16), plot.subtitle = element_text(size = 18, face = "bold"))

ms_b <- ms %>% filter(substeps == 1) %>%
  mutate(lab = factor(start, levels = c("warm start (earlier solution)", "conservative start",
                                        "aggressive start", "cold start (GAMS controls)")))
pA2b <- ggplot(ms_b, aes(year, miu, colour = lab, linetype = lab, linewidth = lab)) +
  geom_line() +
  scale_color_manual(values = c("grey35", "#b8860b", "#2a9d8f", "#5b7c99"), name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "solid"), name = NULL) +
  scale_linewidth_manual(values = c(2.2, 1.4, 1.4, 1.2), name = NULL) +
  x_full +
  coord_cartesian(ylim = c(0, 1.15)) +
  guides(colour = guide_legend(ncol = 2)) +
  labs(x = "Year", y = "DICE-CP emissions\ncontrol rate (μ)", tag = "(b)",
       subtitle = "Annual grid (1 substep per year)") +
  base_theme + theme(legend.position = "top", legend.key.width = unit(2.2, "cm"),
                     legend.text = element_text(size = 16), plot.subtitle = element_text(size = 18, face = "bold"))

save_fig(pA2a / pA2b, "figA2_multistart_dicecp", 16, 14)

ms_tab <- ms %>% distinct(substeps, start, Uacc2420) %>% arrange(desc(substeps), desc(Uacc2420))
cold_holiday <- ms %>% filter(start == "cold start (GAMS controls)", year > 2200, miu < 0.5) %>% pull(year)
note("## Figure A2 (DICE-CP multiple starts)\n",
     "- (a) At the reporting resolution every start that converged reaches the same optimum.\n",
     "- (b) On the annual grid, the cold start (GAMS controls) ends in a worse local optimum: μ collapses ",
     "to near zero around ", min(cold_holiday), "–", max(cold_holiday), ".\n",
     "- U_acc(2300) by run:\n\n| Substeps/yr | Start | U_acc(2300) |\n|---|---|---|\n",
     paste0("| ", ms_tab$substeps, " | ", ms_tab$start, " | ",
            ifelse(is.na(ms_tab$Uacc2420), "trajectory not saved", sprintf("%.1f", ms_tab$Uacc2420)), " |", collapse = "\n"),
     "\n\n- A third 2-substep start (conservative) did not converge and is not shown.\n")

}  # MAKE_APPENDIX

writeLines(c("# Caption notes (generated by ch3_figures_corrected.R; numbers computed from figdata_corrected/)\n", notes),
           file.path(out_dir, "caption_notes.md"))
message("All figures written to ", normalizePath(out_dir))
