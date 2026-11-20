# Pairwise Shapley interaction indices for the 2020 SCC (Lenton variant iv, welfare to 2200), shown as a heatmap,
# with each process's own Shapley value on the diagonal. Input: shapley_interactions.csv, shapley_values.csv.
suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr); library(ggplot2)})
d <- "vis_corrected/feedback_factorial"
I <- read_csv(file.path(d, "shapley_interactions.csv"), show_col_types = FALSE) %>% filter(metric == "SCC2020")
S <- read_csv(file.path(d, "shapley_values.csv"), show_col_types = FALSE) %>% filter(metric == "SCC2020")
lv <- c("ocean", "oceanT", "fert", "photoT", "vegT", "soilT", "perma")
nm <- c(ocean = "Ocean buffer", oceanT = "Ocean warming", fert = "CO2 fertilisation", photoT = "Photosynthesis (T)",
        vegT = "Plant respiration", soilT = "Soil respiration", perma = "Permafrost")
M <- bind_rows(I %>% select(i, j, value = interaction, sign_consistent), I %>% select(i = j, j = i, value = interaction, sign_consistent),
               S %>% transmute(i = process, j = process, value = shapley, sign_consistent = TRUE)) %>%
  mutate(i = factor(nm[i], levels = nm[lv]), j = factor(nm[j], levels = rev(nm[lv])),
         diag = as.character(i) == as.character(j))
p <- ggplot(M, aes(i, j, fill = value)) + geom_tile(color = "white") +
  geom_tile(data = M %>% filter(diag), color = "black", linewidth = 0.6) +
  geom_text(aes(label = ifelse(sign_consistent, sprintf("%+.1f", value), sprintf("%+.1f*", value)),
                fontface = ifelse(sign_consistent, "bold", "plain"), colour = sign_consistent), size = 3.8) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey45"), guide = "none") +
  scale_fill_gradient2(low = "#1f77b4", mid = "white", high = "#d62728", midpoint = 0, name = "$ per tCO2") +
  labs(x = NULL, y = NULL, title = "SCC 2020: own contribution (diagonal) and pairwise interactions",
       caption = "* sign changes across the 32 settings of the other five processes (bold: same sign in all 32)") +
  theme_minimal(base_size = 13) + theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
ggsave(file.path(d, "fig_scc_interactions.png"), p, width = 9, height = 7.5, dpi = 300)
