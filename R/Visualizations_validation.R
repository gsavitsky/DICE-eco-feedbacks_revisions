# Model validation figure: (a-c) C4MIP carbon-cycle feedbacks (1%/yr to 4xCO2) of Lenton (variant iv) and DFAIR against
# the 11 CMIP6 Earth-system models of Arora et al. (2020, Table A1); (d) emissions-driven SSP check, CO2 rise 2020-2100,
# against RCMIP/MAGICC (ssp_check_v2.py).
suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)})
A <- read_csv("vis_corrected/arora2020_cmip6_4xco2.csv", show_col_types = FALSE) %>%
  mutate(beta_total = beta_land + beta_ocean, gamma_total = gamma_land + gamma_ocean)
C <- read_csv("vis_corrected/feedback_factorial/c4mip_compare.csv", show_col_types = FALSE)
len <- C %>% filter(grepl("^Lenton", model)); dfa <- C %>% filter(grepl("^DFAIR", model))
pal <- c("CMIP6 ESMs" = "grey55", "CMIP6 mean" = "black", "Lenton (iv)" = "#d62728", "DFAIR (DICE 2023)" = "#1f77b4")
sc <- function(xv, yv, xl, yl, title, extra = NULL) {
  m <- tibble(x = mean(A[[xv]]), y = mean(A[[yv]]))
  p <- ggplot() +
    geom_point(data = A, aes(.data[[xv]], .data[[yv]], colour = "CMIP6 ESMs"), size = 2.4) +
    geom_point(data = m, aes(x, y, colour = "CMIP6 mean"), shape = 4, size = 4, stroke = 1.4) +
    geom_point(data = len, aes(.data[[xv]], .data[[yv]], colour = "Lenton (iv)"), size = 4, shape = 17)
  if (!is.null(extra)) p <- p + geom_point(data = extra, aes(.data[[xv]], .data[[yv]], colour = "DFAIR (DICE 2023)"), size = 4, shape = 15)
  p + scale_colour_manual(values = pal, name = NULL) + labs(x = xl, y = yl, title = title) +
    theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 12), panel.grid.minor = element_blank())
}
a <- sc("beta_land", "gamma_land", "Land β (PgC/ppm)", "Land γ (PgC/°C)", "(a) Land")
b <- sc("beta_ocean", "gamma_ocean", "Ocean β (PgC/ppm)", "Ocean γ (PgC/°C)", "(b) Ocean")
c <- sc("beta_total", "gamma_total", "Total β (PgC/ppm)", "Total γ (PgC/°C)", "(c) Land + ocean", extra = dfa)
S <- read_csv("vis_corrected/ssp_check_v2.csv", show_col_types = FALSE) %>%   # copied from permafrost_revision_test
  mutate(model = recode(model, "RCMIP/MAGICC concentrations" = "RCMIP / MAGICC", "DFAIR (DICE 2023)" = "DFAIR",
                        "Lenton [Lenton iv]" = "Lenton (iv)", "Lenton + permafrost [Lenton iv]" = "Lenton (iv) + permafrost"),
         model = factor(model, levels = c("RCMIP / MAGICC", "DFAIR", "DFAIR + permafrost", "Lenton (iv)", "Lenton (iv) + permafrost")),
         scenario = toupper(sub("ssp(\\d)(\\d)(\\d)", "SSP\\1-\\2.\\3", scenario)))
d <- ggplot(S, aes(scenario, d2100, fill = model)) + geom_col(position = position_dodge(0.85), width = 0.8) +
  geom_hline(yintercept = 0, colour = "grey40") +
  scale_fill_manual(values = c("RCMIP / MAGICC" = "grey25", "DFAIR" = "#9ecae1", "DFAIR + permafrost" = "#1f77b4",
                               "Lenton (iv)" = "#fcae91", "Lenton (iv) + permafrost" = "#d62728"), name = NULL) +
  labs(x = NULL, y = "CO2 rise 2020-2100 (ppm)", title = "(d) Emissions-driven SSPs (RCMIP v5.1.0)") +
  theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 12), legend.position = "bottom",
                                   panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
fig <- ((a + theme(legend.position = "none")) | (b + theme(legend.position = "none")) | (c + theme(legend.position = "right"))) /
  d + plot_layout(heights = c(1, 1.1))
ggsave("vis_corrected/figures/fig_model_validation.png", fig, width = 14, height = 9.5, dpi = 300)
